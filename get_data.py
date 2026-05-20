import pandas as pd
import json
import os

from pymongo import MongoClient 
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
def get_std_mean(vis_df):
    mean = vis_df["score"].mean()
    std = vis_df["score"].std()
    return std, mean
    

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

#date

def get_date(path):
    filename = os.path.basename(path)
    parts= filename.split("_")
    for part in parts:
        if len(part)==8 and part.isdigit():
            return f"{part[:4]}-{part[4:6]}-{part[6:]}"
    return ""


#versionning
def get_version(sequence,collection):
    documents= list(collection.find({"sequence":sequence}))
    if len(documents) == 0:
        return "1.0"
    else:
        versions= [float(doc["vers"]) for doc in documents if doc.get("vers")]
        derniere_version = max(versions)
        nouvelle_version = round (derniere_version + 0.1, 1)
        return str(nouvelle_version)
    
    
def prepare_send_to_mongo(metrics, top_level_info, frames, avg, std, collection):
    
    last_pdb = frames[-1]["pdb_path"] if frames else ""
    document = {
        "sequence": top_level_info.get("sequence", ""),
        "name": top_level_info.get("name", ""),
        "organism": top_level_info.get("organism", ""),
        "date": get_date(metrics.get("local_filepath", "")),
        "vers": get_version(top_level_info.get("sequence", ""), collection),
        "file": last_pdb,
        "metrics": metrics,
        "RMSD_avg":avg,
        "RMSD_std":std,
        "interframes": frames
    }
    return document


def main():
    
    client = MongoClient(host="localhost", port=27017)
    collection = client["rna_optimizer"]["structures"]

    origin_path = "/home/paul/Documents/Evry/Stage/Optimize_3D_ARNStructure/"
    
    # Method,Score_Function,Sequence_Length,Bead_Atom,Wall_Time_s,GPU_Time_s,Final_Score,Best_Score_Step,Molecule,Out_Name,Potential,Bond,Vis_Dir
    source_file = "/home/paul/Documents/Evry/Stage/Optimize_3D_ARNStructure/metrics.csv"
    source = pd.read_csv(source_file, skipinitialspace=True)
    
    all_documents = []
    
    for i,row in source.iterrows():
        path = os.path.join(origin_path,row["Vis_Dir"])
        #print(row["Vis_Dir"])
        vis_df = read_folding_vis(path)
        metrics_df = read_metric(path)

        frames = parse_vis(vis_df)
        std,mean = get_std_mean(vis_df)
        metrics_dict, top_level_info = parse_metrics(metrics_df)
        
        final_document = prepare_send_to_mongo(metrics_dict, top_level_info, frames, mean, std, collection)
        all_documents.append(final_document)
        #print(f"{json.dumps(final_document, indent=4)}")

    output_file = "mongo_insert.json"
    with open(output_file, 'w') as f:
        json.dump(all_documents, f, indent=4)
    
    collection.insert_many(all_documents)   
    client.close()

    print(f"{len(all_documents)} documents inseres dans Mongodb")

    print(f"Successfully wrote {len(all_documents)} documents to {output_file}")

if __name__ == '__main__':
    main()