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
    return df

def read_metric(path):
    csv_path = path + '/metrics.csv'
    df = pd.read_csv(csv_path, skipinitialspace=True)
    df.columns = df.columns.str.strip()
    df = df.map(lambda x: x.strip() if isinstance(x, str) else x)
    return df

def parse_vis(vis_df):
    frames = []
    for _, row in vis_df.iterrows():
        interframe = {
            "phase": str(row.iloc[0]),
            "epoch": str(row.iloc[1]),
            "score": float(row.iloc[2]),
            "pdb_path": str(row.iloc[3])
        }
        frames.append(interframe)
    return frames 

def get_std_mean(vis_df):
    mean = vis_df["score"].mean()
    std = vis_df["score"].std()
    return std, mean
    

def parse_metrics(metrics_df):
    row = metrics_df.iloc[0]

    #verfication que la sequence contient que des lettres ARN valides 
    sequence = str(row.get('Sequence', ''))
    if not sequence or not all(c in 'augc' for c in sequence.lower()):
        print(f"Warning: Sequence '{sequence}' contains illegal letters or is empty. Skipping.")
        return None, None
    
    document = {
        "methods": str(row.get('Method', '')),
        "optimization_mode": str(row.get('Optimization_Mode','')),
        "score_function": str(row.get('Score_Function', '')),
        "score_weights": str(row.get('Score_Weights', '')),
        "length": int(row.get('Sequence_Length', 0) or 0),
        "bead_atom": str(row.get('Bead_Atom', '')),
        "chain": str(row.get('Chain', '')),
        "time": float(row.get('Wall_Time_s', 0.0) or 0.0),
        "gpu_time": float(row.get('GPU_Time_s', 0.0) or 0.0),
        "video_path": "folding_animation.mp4",
        "final_score": float(row.get('Final_Score', 0.0) or 0.0),
        "best_score_step": int(row.get('Best_Score_Step', 0) or 0),
        "molecule": str(row.get('Molecule', '')),
        "local_filepath": str(row.get('Out_Name', '')),
        "potential": float(row.get('Potential', 0.0) or 0.0),
        "bond": float(row.get('Bond', 0.0) or 0.0),
        "wca": float(row.get('WCA', 0.0) or 0.0),
        "rmsd": float(row.get('RMSD', 0.0) or 0.0),
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
            out_name = str(row.get('Out_Name', ''))
            if not out_name or not out_name.endswith('.pdb'):
                continue
                
            vis_dir_val = str(row.get('Vis_Dir', '')).strip()
            vis_dir_path = None
            if vis_dir_val and vis_dir_val.lower() != 'nan':
                vis_dir_path = os.path.join(origin_path, vis_dir_val)
                if not os.path.isdir(vis_dir_path):
                    vis_dir_path = None
                    
            if not vis_dir_path:
                # e.g. out_name: outputs/opt_bs_C4'_cgRNASP_Seq_1_161850_793965_228916.pdb
                # we need to extract Seq_1_161850_793965_228916
                import re
                filename = os.path.basename(out_name)
                # Match anything that starts with Seq_ followed by digits/underscores
                m = re.search(r'(Seq_[0-9]+_[0-9_]+)', filename)
                if not m:
                    continue
                
                seq_part = m.group(1)
                # clean up trailing _full_atom if matched
                seq_part = seq_part.replace('_full_atom', '').strip('_')
                
                bead_atom = str(row.get('Bead_Atom', "C4'"))
                vis_dir_name = f"vis_{seq_part}_{bead_atom}"
                
                out_dir = os.path.dirname(out_name)
                vis_dir_path = os.path.join(origin_path, out_dir, vis_dir_name)
                
                if not os.path.isdir(vis_dir_path):
                    # Try fallback without bead_atom if not found
                    vis_dir_name = f"vis_{seq_part}"
                    vis_dir_path = os.path.join(origin_path, out_dir, vis_dir_name)
                    if not os.path.isdir(vis_dir_path):
                        continue
            
            try:
                vis_df = read_folding_vis(vis_dir_path)
                frames = parse_vis(vis_df)
                std, mean = get_std_mean(vis_df)
                
                seq_val = args.sequence if args.sequence else str(row.get('Sequence', ''))
                if not seq_val:
                    # If we really don't have a sequence from args or csv, skip
                    continue
                    
                metrics_dict = {
                    "methods": str(row.get('Method', '')),
                    "score_function": str(row.get('Score_Function', '')),
                    "score_weights": str(row.get('Score_Weights', '')),
                    "optimization_mode": str(row.get('Optimization_Mode', '')),
                    "length": int(row.get('Sequence_Length', 0) or 0),
                    "bead_atom": str(row.get('Bead_Atom', '')),
                    "chain": str(row.get('Chain', 'A')),
                    "time": float(row.get('Wall_Time_s', 0.0) or 0.0),
                    "gpu_time": float(row.get('GPU_Time_s', 0.0) or 0.0),
                    "video_path": "folding_animation.mp4",
                    "final_score": float(row.get('Final_Score', 0.0) or 0.0),
                    "best_score_step": int(row.get('Best_Score_Step', 0) or 0),
                    "molecule": str(row.get('Molecule', 'RNA')),
                    "local_filepath": out_name,
                    "potential": float(row.get('Potential', 0.0) or 0.0),
                    "bond": float(row.get('Bond', 0.0) or 0.0),
                    "wca": float(row.get('WCA', 0.0) or 0.0),
                    "rmsd": float(row.get('RMSD', 0.0) or 0.0),
                    "rmsd_bead": str(row.get('RMSD_bead', '')),
                    "vis_dir": str(row.get('Vis_Dir', '')),
                    "type": str(row.get('Type', ''))
                }
                
                if metrics_dict["type"].lower() == "switch":
                    metrics_dict["pdb_initial"] = str(row.get('pdb_initial', ''))
                    metrics_dict["cible"] = str(row.get('pdb_cible', ''))
                
                top_level_info = {
                    "sequence": seq_val,
                    "name": str(row.get('Name_Seq', 'N/A')),
                    "organism": str(row.get('Organism', 'N/A')),
                }
                
                final_document = prepare_send_to_arango(metrics_dict, top_level_info, frames, mean, std, db, all_documents)
                if final_document is not None:
                    all_documents.append(final_document)
                    
            except Exception as e:
                print(f"Error processing {vis_dir_path}: {e}")
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
