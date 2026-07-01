from botocore import parsers  # pyright: ignore[reportMissingImports]
from annotated_types import doc  # pyright: ignore[reportMissingImports]
import pymongo  # pyright: ignore[reportMissingImports]
import pandas as pd  # pyright: ignore[reportMissingImports]
import json
import os
from dotenv import load_dotenv  # pyright: ignore[reportMissingImports]
from pymongo import MongoClient  # pyright: ignore[reportMissingImports]
from pymongo import AsyncMongoClient  # pyright: ignore[reportMissingImports]
from pprint import pprint


'''
This .py file will be called by the main task, it will be given a path to the specific vis_  dir.
In this dir are found :
    folding_vis.csv
    metrics.csv
From those we can get all needed infos. 
Those infos will be inserted into a mongoDB database.
'''


def read_folding_vis(path):
    csv_path = path + '/folding_vis.csv'
    df = pd.read_csv(csv_path, skipinitialspace=True)
    if len(df.columns) > 0:
        # Strip and lowercase columns for uniform access
        df.columns = df.columns.str.strip().str.lower()
        
        # 1. phase
        phase_col = "phase" if "phase" in df.columns else df.columns[0]
        # 2. epoch
        epoch_col = "epoch" if "epoch" in df.columns else df.columns[1]
        
        # 3. score (evolution value)
        score_col = None
        for name in ["score", "rmsd", "value"]:
            if name in df.columns:
                score_col = name
                break
        if score_col is None:
            # Fallback to index 5 (score) if it exists, otherwise index 2
            score_col = df.columns[5] if len(df.columns) > 5 else df.columns[2]
            
        # 4. pdb_path
        pdb_col = None
        for name in ["path_fichier_pdb", "pdb_path", "filepath", "path"]:
            if name in df.columns:
                pdb_col = name
                break
        if pdb_col is None:
            # Fallback to last column
            pdb_col = df.columns[-1]
            
        # Rebuild a clean DataFrame
        clean_df = pd.DataFrame({
            "phase": df[phase_col],
            "epoch": df[epoch_col],
            "score": pd.to_numeric(df[score_col], errors='coerce'),
            "pdb_path": df[pdb_col]
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
        return int(float(val))
    except (ValueError, TypeError):
        return default

def parse_metrics(metrics_df):
    row = metrics_df.iloc[0]

    #verfication que la sequence contient que des lettres ARN valides 
    sequence = str(row.get('Sequence', ''))
    if not sequence or not all(c in 'augc' for c in sequence.lower()):
        print(f"Warning: Sequence '{sequence}' contains illegal letters or is empty. Skipping.")
        return None, None
    
    document = {
        "methods": str(row.get('Method', '')),
        "score_function": str(row.get('Score_Function', '')),
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
        "bond": safe_float(row.get('Bond', 0.0))
    }
    
    
    top_level_info = {
        "sequence": str(row.get('Sequence', '')),
        "name": str(row.get('Name_Seq', '')),
        "organism": str(row.get('Organism', '')),
    }
    
    return document, top_level_info

def get_date(path):
    import re
    # 1. Search for a standalone 8-digit date string anywhere in path
    m = re.search(r'(?:^|[\/_])(\d{8})(?:$|[\/_])', path)
    if m:
        part = m.group(1)
        return f"{part[:4]}-{part[4:6]}-{part[6:]}"
        
    # 2. Fallback to parsing folder components
    parts = path.replace('\\', '/').split('/')
    for p in parts:
        subparts = p.split('_')
        for part in subparts:
            if len(part) == 8 and part.isdigit():
                return f"{part[:4]}-{part[4:6]}-{part[6:]}"
                
    # 3. Last fallback: original basename split logic
    filename = os.path.basename(path)
    parts = filename.split("_")
    for part in parts:
        if len(part) == 8 and part.isdigit():
            return f"{part[:4]}-{part[4:6]}-{part[6:]}"
    return ""


def check_and_get_version(sequence, bead_atom,chain, new_final_score, collection, all_documents):
    db_documents = list(collection.find({"sequence": sequence, "metrics.bead_atom": bead_atom, "metrics.chain": chain}))
    
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
    
    
def prepare_send_to_mongo(metrics, top_level_info, frames, avg, std, collection, all_documents):
    sequence = top_level_info.get("sequence", "")
    bead_atom = metrics.get("bead_atom", "")
    chain = metrics.get("chain", "")
    try:
        new_final_score = float(metrics.get("final_score", float('inf')))
    except ValueError:
        new_final_score = float('inf')
        
    version = check_and_get_version(sequence, bead_atom,chain, new_final_score, collection, all_documents)
    
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
        "RMSD_avg":avg,
        "RMSD_std":std,
        "interframes": frames
    }
    return document


def main():
    load_dotenv()
    MONGO_URI = os.getenv("API_USER")
    client = MongoClient(MONGO_URI, tls=False)
    collection = client["arn"]["sequences"]
    
    
    origin_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "Optimize_3D_ARNStructure"))
       
    csv_file = "metrics.csv"
    source_file = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "Optimize_3D_ARNStructure", csv_file))


    vis_dirs = []
    if os.path.exists(source_file):
        import csv
        with open(source_file, 'r', encoding='utf-8') as f:
            reader = csv.reader(f)
            header = next(reader, None)
            if header:
                for row in reader:
                    if not row:
                        continue
                    last_val = row[-1].strip()
                    if last_val and ("vis_" in last_val or "outputs/" in last_val) and not last_val.endswith(".pdb") and not last_val.endswith(".cif"):
                        vis_dirs.append({"Vis_Dir": last_val})
    source = pd.DataFrame(vis_dirs)
    
    all_documents = []
    
    for i,row in source.iterrows():
        path = os.path.join(origin_path,row["Vis_Dir"])
        try:
            vis_df = read_folding_vis(path)
            metrics_df = read_metric(path)

            frames = parse_vis(vis_df)
            std,mean = get_std_mean(vis_df)
            metrics_dict, top_level_info = parse_metrics(metrics_df)
            
            if metrics_dict is None:
                continue
            
            final_document = prepare_send_to_mongo(metrics_dict, top_level_info, frames, mean, std, collection, all_documents)
            if final_document is not None:
                all_documents.append(final_document)
        except Exception as e:
            print(f"Error processing {path}: {e}")
            continue

        #print(f"{json.dumps(final_document, indent=4)}")

    output_file = "mongo_insert.json"
    with open(output_file, 'w') as f:
        json.dump(all_documents, f, indent=4)
        
    if all_documents:
        collection.insert_many(all_documents)
        print(f"{len(all_documents)} documents insérés dans MongoDB")
    else:
        print("No new documents to insert into MongoDB (all were discarded as identical or worse).")
        
    client.close()

    #print(f"Successfully wrote {len(all_documents)} documents to {output_file}")

if __name__ == '__main__':
    main()