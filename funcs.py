from urllib.request import urlretrieve
import zipfile as zf
import math

def one_dimension(sequence):
    rw = {"ag-ct": [], "ac-gt": [], "at-cg": []}
    #rw["ag-ct"] = []
    tracker = 0
    for base in sequence:
        if base in ["A", "G"]:
            tracker += 1
        else:
            tracker -= 1
        rw["ag-ct"].append(tracker)
    #rw["ac-gt"] = []
    tracker = 0
    for base in sequence:
        if base in ["A", "C"]:
            tracker += 1
        else:
            tracker -= 1
        rw["ac-gt"].append(tracker)
    #rw["at-cg"] = []
    tracker = 0
    for base in sequence:
        if base in ["A", "T"]:
            tracker += 1
        else:
            tracker -= 1
        rw["at-cg"].append(tracker)
    return rw

def two_dimension(sequence):
    rw = {}
    rw["at-cg"] = [[],[]]
    rw["ag-ct"] = [[], []]
    rw["ac-tg"] = [[], []]
    atcg_x_tracker = 0
    atcg_y_tracker = 0
    actg_x_tracker = 0
    actg_y_tracker = 0
    agct_x_tracker = 0
    agct_y_tracker = 0
    for base in sequence:
        match base:
            case "A":
                atcg_x_tracker += 1
                actg_x_tracker += 1
                agct_x_tracker += 1
            case "G":
                atcg_y_tracker += 1
                actg_y_tracker -= 1
                agct_x_tracker -= 1
            case "C":
                atcg_y_tracker -= 1
                actg_x_tracker -= 1
                agct_y_tracker -= 1
            case "T":
                atcg_x_tracker -= 1
                actg_y_tracker += 1
                agct_y_tracker += 1
        rw["at-cg"][0].append(atcg_x_tracker)
        rw["at-cg"][1].append(atcg_y_tracker)
        rw["ag-ct"][0].append(agct_x_tracker)
        rw["ag-ct"][1].append(agct_y_tracker)
        rw["ac-tg"][0].append(actg_x_tracker)
        rw["ac-tg"][1].append(actg_y_tracker)
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
    k = math.log(n)/(math.log(n) + math.log(max_d/l))
    return k

def file_import(gen_id):
    url = (
        f"https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{gen_id}/download?include_annotation_type=FASTA_RNA"
    )
    filename = f"static/{gen_id}.zip"
    urlretrieve(url, filename)

    with zf.ZipFile(f"static/{gen_id}.zip", "r") as seq:
        seq.extract("ncbi_dataset/data/rna.fna", path=f"sequences/{gen_id}")

def seq_extract(gen_id):
    with open(f"sequences/{gen_id}/ncbi_dataset/data/rna.fna", "r") as rna:
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
