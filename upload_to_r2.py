import os
import glob
import boto3
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

def upload_directory_to_r2(directory_path, bucket_name):
    """Uploads all .pdb and .cif files from the given directory to R2."""
    s3_client = get_r2_client()
    if not s3_client:
        return

    # Only upload files from directories starting with "vis"
    search_pattern_pdb = os.path.join(directory_path, "vis*", "**", "*.pdb")
    search_pattern_cif = os.path.join(directory_path, "vis*", "**", "*.cif")
    
    files_to_upload = glob.glob(search_pattern_pdb, recursive=True) + glob.glob(search_pattern_cif, recursive=True)

    if not files_to_upload:
        print(f"No .pdb or .cif files found in {directory_path}")
        return

    print(f"Found {len(files_to_upload)} files to upload to R2 bucket '{bucket_name}'...")

    for file_path in files_to_upload:
        relative_path = os.path.relpath(file_path, directory_path)
        
        # Place everything under the PDB/ folder in the bucket
        object_key = f"PDB/{relative_path.replace(os.sep, '/')}"

        try:
            print(f"Uploading {object_key}...")
            content_type = "text/plain" 
            if file_path.endswith(".pdb"):
                content_type = "chemical/x-pdb"
            elif file_path.endswith(".cif"):
                content_type = "chemical/x-cif"

            s3_client.upload_file(
                file_path, 
                bucket_name, 
                object_key,
                ExtraArgs={'ContentType': content_type}
            )
            print(f"  -> Successfully uploaded {object_key}")
        except FileNotFoundError:
            print(f"  -> Error: File {file_path} not found.")
        except NoCredentialsError:
            print("  -> Error: Credentials not available. Check your .env file.")
            break
        except ClientError as e:
            print(f"  -> Error uploading {object_key}: {e}")

    print("Upload process finished.")

if __name__ == "__main__":
    upload_directory_to_r2(OUTPUTS_DIR, R2_BUCKET_NAME)
