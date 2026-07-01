from arango import ArangoClient  # pyright: ignore[reportMissingImports]
import pandas as pd  # pyright: ignore[reportMissingImports]
import json
import os
from dotenv import load_dotenv  # pyright: ignore[reportMissingImports]
from pprint import pprint


'''
This .py file will be called by the main task, it will be given a path to the specific vis_  dir.
In this dir are found :
    folding_vis.csv
    metrics.csv
From those we can get all needed infos. 
Those infos will be inserted into an ArangoDB database.
'''


def read_folding_vis(path):
    csv_path = path + '/folding_vis.csv'
    df = pd.read_csv(csv_path, skipinitialspace=True)
    if len(df.columns) > 2:
        # Rebuild a clean DataFrame with known columns to avoid duplicate/extra column issues
        clean_df = pd.DataFrame({
            "phase": df.iloc[:, 0],
            "epoch": df.iloc[:, 1],
            "score": pd.to_numeric(df.iloc[:, 2], errors='coerce'),
            "pdb_path": df.iloc[:, 3] if len(df.columns) > 3 else ""
        })
        df = clean_df.dropna(subset=["score"])
    return df

def read_metric(path):
    csv_path = path + '/metrics.csv'
    df = pd.read_csv(csv_path, skipinitialspace=True)
    df.columns = df.columns.str.strip()
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)
    return df

def parse_vis(vis_df):
    frames = []
    import math
    for _, row in vis_df.iterrows():
        score_val = float(row.iloc[2])
        if pd.isna(score_val) or not math.isfinite(score_val):
            score_val = 0.0
        interframe = {
            "phase": str(row.iloc[0]),
            "epoch": str(row.iloc[1]),
            "score": score_val,
            "pdb_path": str(row.iloc[3])
        }
        frames.append(interframe)
    return frames 

def get_std_mean(vis_df):
    import math
    mean = vis_df["score"].mean()
    std = vis_df["score"].std()
    
    if pd.isna(mean) or not math.isfinite(mean):
        mean = 0.0
    if pd.isna(std) or not math.isfinite(std):
        std = 0.0
        
    return std, mean
    

def safe_float(val, default=0.0):
    if pd.isna(val):
        return default
    try:
        fval = float(val)
        import math
        if not math.isfinite(fval):
            return default
        return fval
    except (ValueError, TypeError):
        return default

def safe_int(val, default=0):
    if pd.isna(val):
        return default
    try:
        return int(float(val))  # handle cases like "1.0" being converted to int
    except (ValueError, TypeError):
        return default

def parse_metrics(metrics_df):
    row = metrics_df.iloc[0]

    #verfication que la sequence contient que des lettres ARN valides 
    sequence = str(row.get('Sequence', ''))
    if sequence and not all(c in 'augc' for c in sequence.lower()):
        print(f"Warning: Sequence '{sequence}' contains non-RNA letters.")
    
    document = {
        "methods": str(row.get('Method', '')),
        "optimization_mode": str(row.get('Optimization_Mode','')),
        "score_function": str(row.get('Score_Function', '')),
        "score_weights": str(row.get('Score_Weights', '')),
        "length": safe_int(row.get('Sequence_Length', 0)),
        "bead_atom": str(row.get('Bead_Atom', '')),
        "chain": str(row.get('Chain', '')),
        "time": safe_float(row.get('Wall_Time_s', 0.0)),
        "gpu_time": safe_float(row.get('GPU_Time_s', 0.0)),
        "video_path": "folding_animation.mp4",
        "final_score": safe_float(row.get('Final_Score', 0.0)),
        "best_score_step": safe_int(row.get('Best_Score_Step', 0)),
        "molecule": str(row.get('Molecule', '')),
        "local_filepath": str(row.get('Out_Name', '')),
        "potential": safe_float(row.get('Potential', 0.0)),
        "bond": safe_float(row.get('Bond', 0.0)),
        "wca": safe_float(row.get('WCA', 0.0)),
        "rmsd": safe_float(row.get('RMSD', 0.0)),
        "rmsd_bead": str(row.get('RMSD_bead', '')),
        "vis_dir": str(row.get('Vis_Dir', '')),
        "type": str(row.get('Type', ''))
    }
    
    if document["type"].lower() == "switch":
        document["pdb_initial"] = str(row.get('pdb_initial', ''))
        document["cible"] = str(row.get('pdb_cible', ''))
    
    
    top_level_info = {
        "sequence": str(row.get('Sequence', '')),
        "name": str(row.get('Name_Seq', '')),
        "organism": str(row.get('Organism', '')),
    }
    
    return document, top_level_info

def get_date(path):
    filename = os.path.basename(path)
    parts= filename.split("_")
    for part in parts:
        if len(part)==8 and part.isdigit():
            return f"{part[:4]}-{part[4:6]}-{part[6:]}"
    return ""


def check_and_get_version(sequence, bead_atom, chain, new_final_score, db, all_documents):
    # AQL query to retrieve documents matching sequence, bead_atom, and chain
    query = """
    FOR doc IN sequences
        FILTER doc.sequence == @sequence 
          AND doc.metrics.bead_atom == @bead_atom 
          AND doc.metrics.chain == @chain
        RETURN doc
    """
    cursor = db.aql.execute(query, bind_vars={
        "sequence": sequence,
        "bead_atom": bead_atom,
        "chain": chain
    })
    db_documents = list(cursor)
    
    local_documents = [doc for doc in all_documents if doc.get("sequence") == sequence and doc.get("metrics", {}).get("bead_atom") == bead_atom and doc.get("metrics", {}).get("chain") == chain]
    
    all_matching_docs = db_documents + local_documents
    
    if len(all_matching_docs) == 0:
        return "1.0"
        
    best_existing_score = float('inf')
    max_version = 0.0
    
    for doc in all_matching_docs:
        vers_str = doc.get("vers")
        if not vers_str:
            vers_str = "1.0"
            
        try:
            vers = float(vers_str)
        except ValueError:
            vers = 1.0
            
        if vers > max_version:
            max_version = vers
            
        metrics = doc.get("metrics", {})
        score_str = metrics.get("final_score")
        if score_str:
            try:
                score = float(score_str)
                if score < best_existing_score:
                    best_existing_score = score
            except ValueError:
                pass
                            
    if new_final_score < best_existing_score:
        print("yipeeeee")
        return str(round(max_version + 0.1, 1))
    else:
        return None
    
    
def prepare_send_to_arango(metrics, top_level_info, frames, avg, std, db, all_documents):
    sequence = top_level_info.get("sequence", "")
    bead_atom = metrics.get("bead_atom", "")
    chain = metrics.get("chain", "")
    try:
        new_final_score = float(metrics.get("final_score", float('inf')))
    except ValueError:
        new_final_score = float('inf')
        
    version = check_and_get_version(sequence, bead_atom, chain, new_final_score, db, all_documents)
    
    if version is None:
        return None
        
    last_pdb = frames[-1]["pdb_path"] if frames else ""
    document = {
        "sequence": sequence,
        "name": top_level_info.get("name", ""),
        "organism": top_level_info.get("organism", ""),
        "date": get_date(metrics.get("local_filepath", "")),
        "vers": version,
        "file": last_pdb,
        "metrics": metrics,
        "RMSD_avg": avg,
        "RMSD_std": std,
        "interframes": frames
    }
    return document


def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence", type=str, required=False, default="", help="Sequence string passed by API")
    args = parser.parse_args()
    
    load_dotenv()
    ARANGO_URL = os.getenv("ARANGO_URL", "http://localhost:8529")
    
    url_stripped = ARANGO_URL.strip().rstrip('/')
    if url_stripped in ("https:/", "https:", "https://", "http:/", "http:", "http://") or not url_stripped:
        print(f"Error: Invalid ARANGO_URL '{ARANGO_URL}' configured in .env.")
        return

    ARANGO_USER = os.getenv("ARANGO_USER", "root")
    ARANGO_PASSWORD = os.getenv("ARANGO_PASSWORD", "")
    ARANGO_DB = os.getenv("ARANGO_DB", "ARN")

    client = ArangoClient(hosts=ARANGO_URL)
    
    try:
        db = client.db(ARANGO_DB, username=ARANGO_USER, password=ARANGO_PASSWORD)
        if not db.has_collection("sequences"):
            collection = db.create_collection("sequences")
        else:
            collection = db.collection("sequences")
    except Exception as direct_error:
        try:
            sys_db = client.db('_system', username=ARANGO_USER, password=ARANGO_PASSWORD)
            if not sys_db.has_database(ARANGO_DB):
                sys_db.create_database(ARANGO_DB)
            db = client.db(ARANGO_DB, username=ARANGO_USER, password=ARANGO_PASSWORD)
            if not db.has_collection("sequences"):
                collection = db.create_collection("sequences")
            else:
                collection = db.collection("sequences")
        except Exception as sys_error:
            print(f"Connection Error: {direct_error}")
            raise direct_error
       
    origin_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "Optimize_3D_ARNStructure"))
    csv_file = "metrics.csv"
    source_file = os.path.abspath(os.path.join(origin_path, csv_file))

    all_documents = []
    
    if os.path.exists(source_file):
        df = pd.read_csv(source_file, skipinitialspace=True)
        df.columns = df.columns.str.strip()
        df = df.map(lambda x: x.strip() if isinstance(x, str) else x)
        
        for i, row in df.iterrows():
            vis_dir_val = str(row.get('Vis_Dir', '')).strip()
            if not vis_dir_val or vis_dir_val.lower() == 'nan':
                print(f"Row {i} skipped: Vis_Dir is empty or nan.")
                continue
                
            vis_dir_path = os.path.join(origin_path, vis_dir_val)
            if not os.path.isdir(vis_dir_path):
                print(f"Row {i} skipped: Directory {vis_dir_path} does not exist.")
                continue
                
            try:
                vis_df = read_folding_vis(vis_dir_path)
                frames = parse_vis(vis_df)
                std, mean = get_std_mean(vis_df)
                
                # Read local metrics.csv for full prediction-specific details (Sequence, Chain, Organism, Name_Seq, PDB info)
                local_metrics_df = read_metric(vis_dir_path)
                metrics_dict, top_level_info = parse_metrics(local_metrics_df)
                
                if metrics_dict is None:
                    print(f"Warning: Row {i} local metrics_dict is None, skipping.")
                    continue
                
                # Overwrite/enrich with global row keys that may not exist locally (Score_Weights, WCA, RMSD, etc.)
                for key in ["Score_Weights", "Optimization_Mode", "WCA", "RMSD", "RMSD_bead", "Vis_Dir"]:
                    val = row.get(key)
                    if pd.notna(val):
                        field_name = key.lower()
                        if field_name == "score_weights":
                            metrics_dict["score_weights"] = str(val)
                        elif field_name == "optimization_mode":
                            metrics_dict["optimization_mode"] = str(val)
                        elif field_name == "wca":
                            metrics_dict["wca"] = safe_float(val)
                        elif field_name == "rmsd":
                            metrics_dict["rmsd"] = safe_float(val)
                        elif field_name == "rmsd_bead":
                            metrics_dict["rmsd_bead"] = str(val)
                        elif field_name == "vis_dir":
                            metrics_dict["vis_dir"] = str(val)
                    
                if args.sequence:
                    top_level_info["sequence"] = args.sequence
                    
                seq_val = top_level_info.get("sequence", "")
                if not seq_val or str(seq_val).lower() == 'nan':
                    top_level_info["sequence"] = "UNKNOWN"
                    
                if not metrics_dict.get("vis_dir") or str(metrics_dict.get("vis_dir")).lower() == 'nan':
                    metrics_dict["vis_dir"] = vis_dir_val
                    
                final_document = prepare_send_to_arango(metrics_dict, top_level_info, frames, mean, std, db, all_documents)
                if final_document is not None:
                    all_documents.append(final_document)
                    
            except Exception as e:
                print(f"Error processing row {i} with vis_dir {vis_dir_path}: {e}")
                continue

    output_file = "arango_insert.json"
    with open(output_file, 'w') as f:
        json.dump(all_documents, f, indent=4)
        
    if all_documents:
        batch_size = 10
        inserted_count = 0
        for i in range(0, len(all_documents), batch_size):
            batch = all_documents[i:i + batch_size]
            try:
                collection.insert_many(batch)
                inserted_count += len(batch)
            except Exception as e:
                print(f"Batch insert failed: {e}. Retrying one-by-one...")
                for doc in batch:
                    try:
                        collection.insert(doc)
                        inserted_count += 1
                    except Exception as single_e:
                        pass
        print(f"{inserted_count} documents insérés dans ArangoDB")
    else:
        print("No new documents to insert into ArangoDB (all were discarded as identical or worse).")
        
    client.close()

if __name__ == '__main__':
    main()
