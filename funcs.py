from urllib.request import urlretrieve
import zipfile as zf
import numpy as np
import math
import json
import os
from StatTools.analysis.dfa import dfa
from StatTools.analysis.utils import analyse_zero_cross_ff
import matplotlib.pyplot as plt
from scipy import stats
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans

def one_dimension(sequence):
    rw = {"AG-CT": [], "AC-GT": [], "AT-CG": []}
    # cs = {"AG-CT": [], "AC-GT": [], "AT-CG": []}
    tracker = 0
    for base in sequence:
        if base in ["A", "G"]:
            tracker += 1
            # cs["AG-CT"].append(1)
        else:
            tracker -= 1
            # cs["AG-CT"].append(-1)
        rw["AG-CT"].append(tracker)

    tracker = 0
    for base in sequence:
        if base in ["A", "C"]:
            tracker += 1
            # cs["AC-GT"].append(1)
        else:
            tracker -= 1
            # cs["AC-GT"].append(-1)
        rw["AC-GT"].append(tracker)
    tracker = 0
    for base in sequence:
        if base in ["A", "T"]:
            tracker += 1
            # cs["AT-CG"].append(1)
        else:
            tracker -= 1
            # cs["AT-CG"].append(-1)
        rw["AT-CG"].append(tracker)
    return rw#, cs

def two_dimension(sequence):
    rw = {}
    rw["AT-CG"] = [[],[]]
    rw["AG-CT"] = [[], []]
    rw["AC-GT"] = [[], []]
    atcg_x_tracker = 0
    atcg_y_tracker = 0
    acgt_x_tracker = 0
    acgt_y_tracker = 0
    agct_x_tracker = 0
    agct_y_tracker = 0
    for base in sequence:
        match base:
            case "A":
                atcg_x_tracker += 1
                acgt_x_tracker += 1
                agct_x_tracker += 1
            case "G":
                atcg_y_tracker -= 1
                acgt_y_tracker += 1
                agct_x_tracker -= 1
            case "C":
                atcg_y_tracker += 1
                acgt_x_tracker -= 1
                agct_y_tracker += 1
            case "T":
                atcg_x_tracker -= 1
                acgt_y_tracker -= 1
                agct_y_tracker -= 1
        rw["AT-CG"][0].append(atcg_x_tracker)
        rw["AT-CG"][1].append(atcg_y_tracker)
        rw["AG-CT"][0].append(agct_x_tracker)
        rw["AG-CT"][1].append(agct_y_tracker)
        rw["AC-GT"][0].append(acgt_x_tracker)
        rw["AC-GT"][1].append(acgt_y_tracker)
    return rw

def three_dimension(sequence):
    rw = [[], [], []]
    x_tracker = 0
    y_tracker = 0
    z_tracker = 0
    for base in sequence:
        match base:
            case "A":
                x_tracker += 0.75
                y_tracker -= 0.25
                z_tracker -= 0.25
            case "C":
                x_tracker -= 0.25
                y_tracker += 0.75
                z_tracker -= 0.25
            case "G":
                x_tracker -= 0.25
                y_tracker -= 0.25
                z_tracker += 0.75
            case "T":
                pass
        rw[0].append(x_tracker)
        rw[1].append(y_tracker)
        rw[2].append(z_tracker)
    return rw

def mean_pos(series, dimension):
    match dimension:
        case 1:
            total_pos = 0
            for element in series:
                total_pos += element
            mean = total_pos/len(series)
            return mean
        case 2:
            total_pos_x = 0
            total_pos_y = 0
            for i in range(len(series[0])):
                total_pos_x += series[0][i]
                total_pos_y += series[1][i]
            mean_x = total_pos_x/len(series[0])
            mean_y = total_pos_y/len(series[0])
            return mean_x, mean_y
        case 3:
            total_pos_x = 0
            total_pos_y = 0
            total_pos_z = 0
            for i in range(len(series[0])):
                total_pos_x += series[0][i]
                total_pos_y += series[1][i]
                total_pos_z += series[2][i]
            mean_x = total_pos_x/len(series[0])
            mean_y = total_pos_y/len(series[0])
            mean_z = total_pos_z/len(series[0])
            return mean_x, mean_y, mean_z
    print("invalid dimension")
    return 0

def katz(sequence: list, dimension: int):
    match dimension:
        case 1:
            n = len(sequence)-1
            y0 = sequence[0]
            max_d = float('-inf')
            l = 0
            x = 0
            for base in sequence:
                y = base - y0
                d = math.sqrt(x**2+y**2)
                max_d = max(max_d, d)
                x+=1
            for i in range(len(sequence)-1):
                y = sequence[i]-sequence[i+1]
                x = 1
                d = math.sqrt(x**2+y**2)
                l+=d
        case 2:
            n = len(sequence[0])-1
            x0 = sequence[0][0]
            y0 = sequence[1][0]
            max_d = float("-inf")
            l = 0
            for i in range(len(sequence[0])):
                x = sequence[0][i] - x0
                y = sequence[1][i] - y0
                d = math.sqrt(x**2+y**2)
                max_d = max(max_d, d)
            for i in range(len(sequence[0])-1):
                x = sequence[0][i] - sequence[0][i+1]
                y = sequence[1][i] - sequence[1][i+1]
                d = math.sqrt(x ** 2 + y ** 2)
                l += d
        case 3:
            n = len(sequence[0]) - 1
            x0 = sequence[0][0]
            y0 = sequence[1][0]
            z0 = sequence[2][0]
            max_d = float("-inf")
            l = 0
            for i in range(len(sequence[0])):
                x = sequence[0][i]-x0
                y = sequence[1][i] - y0
                z = sequence[2][i] - z0
                d = math.sqrt(z**2 + x**2 + y**2)
                max_d = max(max_d, d)
            for i in range(len(sequence[0])-1):
                x = sequence[0][i] - sequence[0][i + 1]
                y = sequence[1][i] - sequence[1][i + 1]
                z = sequence[2][i] - sequence[2][i + 1]
                d = math.sqrt(x**2 + y**2 + z**2)
                l += d
    k = math.log10(n)/(math.log10(n) + math.log10(max_d/l))
    return k

def dfa_hurst(s):
    s = np.array(s)
    min_window = 4
    max_window = len(s)//4
    demeaned_s = s-np.mean(s)
    cumsum_s = []
    tracker = 0
    for element in demeaned_s:
        tracker += element
        cumsum_s.append(tracker)
    seg_size = np.unique(np.logspace(np.log10(min_window), np.log10(max_window), 20, dtype = int))
    fn = []
    for size in seg_size:
        segments = []
        for i in range(0, len(cumsum_s), size):
            segments.append(cumsum_s[i:i+size])
        rms_vals = []
        for segment in segments:
            if len(segment) == 1:
                continue
            x = np.arange(len(segment))
            fit = np.polyfit(x, segment, 1)
            rms = np.sqrt(np.mean((segment - np.polyval(fit, x))**2))
            rms_vals.append(rms)
        #fn.append(np.sqrt(np.mean(np.square(rms_vals))))
        fn.append(np.mean(rms_vals))
    x = np.log10(seg_size)
    y = np.log10(fn)
    slope, intercept, r, p, stderr = stats.linregress(x, y)
    return slope

def file_import(gen_id):
    if gen_id[:3] == "GCA":
        url = (
            f"https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{gen_id}/download?include_annotation_type=GENOME_FASTA"
        )
        filename = f"static/{gen_id}.zip"
        urlretrieve(url, filename)

        with zf.ZipFile(f"static/{gen_id}.zip", "r") as seq:
            seq.extract("ncbi_dataset/data/assembly_data_report.jsonl", f"sequences/{gen_id}")
            with open(f"sequences/{gen_id}/ncbi_dataset/data/assembly_data_report.jsonl", "r") as f:
                assembly = json.load(f)
                assembly_name = assembly["assemblyInfo"]["assemblyName"]
            seq.extract(f"ncbi_dataset/data/{gen_id}/{gen_id}_{assembly_name}_genomic.fna",
                        path=f"sequences/{gen_id}")

        os.rename(f"sequences/{gen_id}/ncbi_dataset/data/assembly_data_report.jsonl",
                  f"sequences/{gen_id}/data_report.jsonl")
        os.rename(f"sequences/{gen_id}/ncbi_dataset/data/{gen_id}/{gen_id}_{assembly_name}_genomic.fna",
                  f"sequences/{gen_id}/rna.fna")
        os.rmdir(f"sequences/{gen_id}/ncbi_dataset/data/{gen_id}")
        os.rmdir(f"sequences/{gen_id}/ncbi_dataset/data")
        os.rmdir(f"sequences/{gen_id}/ncbi_dataset")

    else:
        url = (
            f"https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{gen_id}/download?include_annotation_type=FASTA_RNA"
        )
        filename = f"static/{gen_id}.zip"
        urlretrieve(url, filename)

        with zf.ZipFile(f"static/{gen_id}.zip", "r") as seq:
            seq.extract("ncbi_dataset/data/rna.fna",
                        f"sequences/{gen_id}")
            seq.extract("ncbi_dataset/data/data_report.jsonl",
                        f"sequences/{gen_id}")

        os.rename(f"sequences/{gen_id}/ncbi_dataset/data/data_report.jsonl",
                  f"sequences/{gen_id}/data_report.jsonl")
        os.rename(f"sequences/{gen_id}/ncbi_dataset/data/rna.fna",
                  f"sequences/{gen_id}/rna.fna")
        os.rmdir(f"sequences/{gen_id}/ncbi_dataset/data")
        os.rmdir(f"sequences/{gen_id}/ncbi_dataset")

def seq_extract(gen_id):
    with open(f"sequences/{gen_id}/rna.fna", "r") as rna:
        raw_seq = rna.read().splitlines()
        sequence = ""
        first = False
        for line in raw_seq:
            if line[0] == ">":
                if first:
                    break
                else:
                    first = True
            else:
                sequence = sequence + line
        return sequence

def get_name(gen_id):
    if gen_id[:3] == "GCA":
        with open(f"sequences/{gen_id}/data_report.jsonl", "r") as n:
            report = json.load(n)
            org_name = report["organism"]["organismName"]
            return org_name+" full genome"
    else:
        with open(f"sequences/{gen_id}/data_report.jsonl", "r") as n:
            report = json.load(n)
            org_name = report["commonName"]
            gene_name = report["description"]
    return org_name + " " + gene_name

def k_means(results, test = False):
    #if test == False:
    labels = [get_name(gen_id) for gen_id in results.keys()]
    all_params = []
    for gen_id, dimensions in results.items():
        params_ext = []
        for dimension, mappings in dimensions.items():
            for mapping, params in mappings.items():
                for param in params.values():
                   params_ext.append(param)
        all_params.append(params_ext)
    print(all_params)
    scaler = StandardScaler()
    all_params_scaled = scaler.fit_transform(all_params)
    print(all_params_scaled)
    n_components = min(len(all_params), len(all_params[0]))
    pca = PCA(n_components)
    all_params_scaled_pca = pca.fit_transform(all_params_scaled)
    print(all_params_scaled_pca)
    max_silhouette = float('-inf')
    opt_n = 0
    opt_clusters = None
    for n in range(2, min(25, len(all_params))):
        kmeans_clusters = KMeans(n_clusters=n, random_state=10).fit_predict(all_params_scaled_pca)
        print(kmeans_clusters)
        silhouette = silhouette_score(all_params_scaled_pca, kmeans_clusters)
        max_silhouette = max(max_silhouette, silhouette)
        if max_silhouette==silhouette:
            opt_n = n
            opt_clusters = kmeans_clusters
    return opt_n, opt_clusters, max_silhouette