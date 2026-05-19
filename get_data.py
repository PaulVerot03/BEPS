import pandas as pd
import json

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
            "score": str(row.iloc[2]),
            "pdb_path": str(row.iloc[3])
        }
        frames.append(interframe)
    return frames

def parse_metrics(metrics_df):
    row = metrics_df.iloc[0]
    
    document = {
        "methods": str(row.get('Method', '')),
        "score_function": str(row.get('Score_Function', '')),
        "length": str(row.get('Sequence_Length', '')),
        "bead_atom": str(row.get('Bead_Atom', '')),
        "chain": str(row.get('Chain', '')),
        "time": str(row.get('Wall_Time_s', '')),
        "gpu_time": str(row.get('GPU_Time_s', '')),
        "video_path": "folding_video.mp4",
        "final_score": str(row.get('Final_Score', '')),
        "best_score_step": str(row.get('Best_Score_Step', '')),
        "molecule": str(row.get('Molecule', '')),
        "local_filepath": str(row.get('Out_Name', '')),
        "potential": str(row.get('Potential', '')),
        "bond": str(row.get('Bond', ''))
    }
    
    # TODO : query la base pour connaitre la version de l'enregistrement et la metre a jour
    top_level_info = {
        "sequence": str(row.get('Sequence', '')),
        "name": str(row.get('Name_Seq', '')),
        "organism": str(row.get('Organism', '')),
    }
    
    return document, top_level_info

def prepare_send_to_mongo(metrics, top_level_info, frames):
    last_pdb = frames[-1]["pdb_path"] if frames else ""
    
    document = {
        "sequence": top_level_info.get("sequence", ""),
        "name": top_level_info.get("name", ""),
        "organism": top_level_info.get("organism", ""),
        "date": "",
        "vers": "",
        "file": last_pdb,
        "metrics": metrics,
        "RMSD": "",
        "interframes": frames
    }
    return document

def main():
    path = "/home/paul/Documents/Evry/Stage/Optimize_3D_ARNStructure/outputs/vis_20260519_114243"
    
    # Read CSVs using pandas
    vis_df = read_folding_vis(path)
    metrics_df = read_metric(path)

    # Parse to dicts
    frames = parse_vis(vis_df)
    metrics_dict, top_level_info = parse_metrics(metrics_df)
    
    # Prepare final document
    final_document = prepare_send_to_mongo(metrics_dict, top_level_info, frames)
    
    print(f"{json.dumps(final_document, indent=4)}")

if __name__ == '__main__':
    main()