import os
import glob
import boto3
import csv
from botocore.exceptions import NoCredentialsError, ClientError
from dotenv import load_dotenv # pyright: ignore[reportMissingImports]

load_dotenv()

R2_ACCOUNT_ID = os.getenv("ACCOUNT_ID")
R2_ACCESS_KEY_ID = os.getenv("R2_ACCESS_KEY_ID")
R2_SECRET_ACCESS_KEY = os.getenv("R2_SECRET_ACCESS_KEY")
R2_BUCKET_NAME = os.getenv("R2_BUCKET_NAME")

OUTPUTS_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "../Optimize_3D_ARNStructure/outputs")

def get_r2_client():
    """Initializes and returns the boto3 client for Cloudflare R2."""
    if not all([R2_ACCOUNT_ID, R2_ACCESS_KEY_ID, R2_SECRET_ACCESS_KEY, R2_BUCKET_NAME]):
        print("Error: Missing R2 credentials in .env file.")
        print("Please ensure R2_ACCOUNT_ID, R2_ACCESS_KEY_ID, R2_SECRET_ACCESS_KEY, and R2_BUCKET_NAME are set.")
        return None

    return boto3.client(
        's3',
        endpoint_url=f"https://{R2_ACCOUNT_ID}.r2.cloudflarestorage.com",
        aws_access_key_id=R2_ACCESS_KEY_ID,
        aws_secret_access_key=R2_SECRET_ACCESS_KEY,
        region_name="auto"
    )

def upload_file_to_r2(s3_client, bucket_name, local_path, object_key, content_type):
    """Uploads a file to R2 with a specific content type."""
    try:
        print(f"Uploading {local_path} to {object_key}...")
        s3_client.upload_file(
            local_path,
            bucket_name,
            object_key,
            ExtraArgs={'ContentType': content_type}
        )
        print(f"  -> Successfully uploaded {object_key}")
        return True
    except FileNotFoundError:
        print(f"  -> Error: File {local_path} not found.")
    except ClientError as e:
        print(f"  -> Error uploading {object_key}: {e}")
    return False

def get_pdb_filename(local_vis_dir):
    """Reads the last line of folding_vis.csv to find the final PDB filename."""
    csv_path = os.path.join(local_vis_dir, "folding_vis.csv")
    if os.path.exists(csv_path):
        try:
            with open(csv_path, mode='r', encoding='utf-8') as f:
                reader = csv.reader(f)
                rows = list(reader)
                if len(rows) > 1:
                    header = [h.strip().lower() for h in rows[0]]
                    
                    # 1. Try exact known names
                    pdb_idx = -1
                    for name in ["path_fichier_pdb", "pdb_path"]:
                        if name in header:
                            pdb_idx = header.index(name)
                            break
                    
                    # 2. Try searching for "pdb" or "path" in headers
                    if pdb_idx == -1:
                        for i, h in enumerate(header):
                            if "pdb" in h or "path" in h:
                                pdb_idx = i
                                break
                    
                    # 3. Fallback to the last column
                    if pdb_idx == -1:
                        pdb_idx = len(header) - 1

                    last_row = rows[-1]
                    if last_row and len(last_row) > pdb_idx:
                        pdb_path = last_row[pdb_idx].strip()
                        if pdb_path:
                            return os.path.basename(pdb_path)
        except Exception as e:
            print(f"  -> Warning: failed to parse {csv_path}: {e}")

    # Fallbacks: check if either common name exists locally
    for fallback in ["best_structure_full_atom.pdb", "vis_best_full_atom.pdb"]:
        if os.path.exists(os.path.join(local_vis_dir, fallback)):
            return fallback
            
    # Default fallback
    return "vis_best_full_atom.pdb"


def upload_directory_to_r2(directory_path, bucket_name):
    """Reads subdirectories from metrics.csv, checks their presence in R2, and conditionally uploads files."""
    s3_client = get_r2_client()
    if not s3_client:
        return

    # metrics.csv is located in the parent directory of outputs
    csv_file_path = os.path.abspath(os.path.join(directory_path, "..", "metrics.csv"))
    if not os.path.exists(csv_file_path):
        print(f"Error: metrics.csv not found at {csv_file_path}")
        return

    print(f"Reading sequence directories from {csv_file_path}...")
    vis_dirs = []
    with open(csv_file_path, mode='r', encoding='utf-8') as f:
        reader = csv.reader(f)
        header = next(reader, None)
        if header:
            header = [h.strip() for h in header]
            try:
                vis_dir_idx = header.index("Vis_Dir")
            except ValueError:
                vis_dir_idx = -1

            if vis_dir_idx != -1:
                for row in reader:
                    if not row or len(row) <= vis_dir_idx:
                        continue
                    vis_dir = row[vis_dir_idx].strip()
                    # Skip empty paths, or PDB file paths
                    if vis_dir and not vis_dir.endswith(".pdb") and not vis_dir.endswith(".cif"):
                        vis_dirs.append(vis_dir)

    # Filter out duplicate directories to save API calls
    unique_vis_dirs = list(dict.fromkeys(vis_dirs))
    print(f"Found {len(unique_vis_dirs)} unique sequence directories to process.")

    for vis_dir in unique_vis_dirs:
        # Resolve the local path of the visualization directory
        local_vis_dir = os.path.abspath(os.path.join(directory_path, "..", vis_dir))
        
        if not os.path.exists(local_vis_dir) or not os.path.isdir(local_vis_dir):
            print(f"Skipping {vis_dir}: Local directory does not exist.")
            continue

        # Get local paths
        pdb_filename = get_pdb_filename(local_vis_dir)
        local_pdb_path = os.path.join(local_vis_dir, pdb_filename)
        local_video_path = os.path.join(local_vis_dir, "folding_animation.mp4")

        if not os.path.exists(local_pdb_path) and not os.path.exists(local_video_path):
            print(f"Skipping {vis_dir}: Neither pdb ({pdb_filename}) nor video exists locally.")
            continue

        # Determine relative path from directory_path (which is outputs)
        sub_rel_path = os.path.relpath(local_vis_dir, directory_path).replace(os.sep, '/')
        prefix = f"PDB/{sub_rel_path}/"
        
        pdb_key = f"PDB/{sub_rel_path}/{pdb_filename}"
        video_key = f"PDB/{sub_rel_path}/folding_animation.mp4"

        try:
            # Query R2 for files starting with this prefix
            response = s3_client.list_objects_v2(Bucket=bucket_name, Prefix=prefix)
            contents = response.get('Contents', [])
            existing_keys = {obj['Key'] for obj in contents}
        except ClientError as e:
            print(f"Error querying R2 for prefix {prefix}: {e}")
            continue

        subdir_exists_on_r2 = len(existing_keys) > 0

        if not subdir_exists_on_r2:
            print(f"Distant subdirectory '{prefix}' does not exist on R2. Uploading PDB and video...")
            if os.path.exists(local_pdb_path):
                upload_file_to_r2(s3_client, bucket_name, local_pdb_path, pdb_key, "chemical/x-pdb")
            if os.path.exists(local_video_path):
                upload_file_to_r2(s3_client, bucket_name, local_video_path, video_key, "video/mp4")
        else:
            # Distant subdirectory exists, check individual files
            if video_key not in existing_keys:
                print(f"Video '{video_key}' is missing on R2. Uploading...")
                if os.path.exists(local_video_path):
                    upload_file_to_r2(s3_client, bucket_name, local_video_path, video_key, "video/mp4")
                else:
                    print(f"  -> Warning: Local video file not found at {local_video_path}")
            
            if pdb_key not in existing_keys:
                print(f"PDB '{pdb_key}' is missing on R2. Uploading...")
                if os.path.exists(local_pdb_path):
                    upload_file_to_r2(s3_client, bucket_name, local_pdb_path, pdb_key, "chemical/x-pdb")
                else:
                    print(f"  -> Warning: Local PDB file not found at {local_pdb_path}")

    print("Upload process finished.")

if __name__ == "__main__":
    upload_directory_to_r2(OUTPUTS_DIR, R2_BUCKET_NAME)
