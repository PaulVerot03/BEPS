import csv
import json

'''
This .py file will be called by the main task, it will be given a path to the specific vis_  dir.
In this dir are found :
    folding_vis.csv
    metrics.csv
From those we can get all needed infos. 
Those infos will be inserted into a mongoDB database.

"sequence":"",
    "name":"",
    "organism":"",
    "date":"",
    "vers":"",
    "file":"",
    "metrics":{
        "methods":"",
        "score_function":"",
        "length":"",
        "bead_atom":"",
        "time":"",
        "gpu_time":"",
        "final_score":"",
        "best_score_step":"",
        "molecule":"",
        "local_filepath":"",
        "potential":"",
        "bond":""
    },
    "RMSD":"",
    "interframes": [
        {
            "phase":"",
            "epoch":"",
            "score":"",
            "video_path":"",
            "pdb_path":""
        }
    ]

'''

def read_folding_vis(path):
    csv_path = path + '/folding_vis.csv'
    data = [] # phase,epoch,score,path_fichier_png,path_fichier_pdb
    with open(csv_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            data.append(row)
    return data

def read_metric(path):
    csv_path = path + '/metrics.csv'
    data = [] 
    # Method,Score_Function,Sequence_Length,Bead_Atom,Wall_Time_s,GPU_Time_s,Final_Score,Best_Score_Step,Molecule,Out_Name,Metric,logs
    with open(csv_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            data.append(row)
    return data

def parse_vis(vis):
    frames = []
    for i in range(len(vis)):
        parts = vis[i].split()
        # parts = phase,epoch,score,path_fichier_png,path_fichier_pdb
        phase = parts[0]
        epoch = parts[1]
        score = parts[2]
        pdb_path = parts[4]
        if i != len(vis)-1:
            interframe = {
                "phase": phase,
                "epoch": epoch,
                "score": score,
                "pdb_path": pdb_path
            },
        if i == len(vis)-1: # avoid trailing comma
            interframe = {
                "phase": phase,
                "epoch": epoch,
                "score": score,
                "pdb_path": pdb_path
            }
            
        frames.append(interframe)
    return(frames)

def parse_metrics(metrics):
    document = {}
    parts = metrics.split()
    method = parts[0]
    score_function = parts[1]
    sequence_length = parts[2]
    bead_atom = parts[3]
    wall_time = parts[4]
    gpu_time = parts[5]
    final_score = parts[6]
    best_score_step = parts[7]
    molecule = parts[8]
    out_name = parts[9]
    potential = parts[10]
    bond = parts[11]
    
    
def prepare_send_to_mongo(metrics, frames):
    last_pdb = frames[len(frames-1)].pdb_path
    document = {
    "sequence":"",
    "name":"",
    "organism":"",
    "date":"",
    "vers":"",
    "file":last_pdb,
    "metrics":{
        metrics
    },
    "RMSD":"",
    "video_path":"",
    "interframes": [
        frames 
    ]
    }
    
def main():
    path = "/home/paul/Documents/Evry/Stage/Optimize_3D_ARNStructure/outputs/vis_20260518_140042"
    vis = read_folding_vis(path)
    metrics = read_metric(path)

    frames = parse_vis(vis)
    infos = parse_metrics(metrics)
    prepare_send_to_mongo(frames,infos)


if __name__ == '__main__':
    main()