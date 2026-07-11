from urllib.request import urlretrieve
import zipfile as zf

def one_dimension(sequences):
    rw = {}
    for seq, code in sequences.items():
        rw[seq] = []
        tracker = 0
        for base in code:
            if base in ["A", "G"]:
                tracker += 1
            else:
                tracker -= 1
            rw[seq].append(tracker)
    return rw

def two_dimension(sequences):
    rw = {}
    for seq, code in sequences.items():
        rw[seq] = [[],[]]
        x_tracker = 0
        y_tracker = 0
        for base in code:
            match base:
                case "A":
                    x_tracker += 1
                case "G":
                    x_tracker -= 1
                case "C":
                    y_tracker -= 1
                case "T":
                    y_tracker += 1
            rw[seq][0].append(x_tracker)
            rw[seq][1].append(y_tracker)
    return rw

def three_dimension(sequences):
    rw = {}
    for seq, code in sequences.items():
        rw[seq] = [[], [], []]
        x_tracker = 0
        y_tracker = 0
        z_tracker = 0
        for base in code:
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
            rw[seq][0].append(x_tracker)
            rw[seq][1].append(y_tracker)
            rw[seq][2].append(z_tracker)
    return rw

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
        sequences = {}
        seq_num = 0
        for line in raw_seq:
            if line[0] == ">":
                seq_num += 1
                sequences[f"sequence {seq_num}"] = ""
            else:
                sequences[f"sequence {seq_num}"] = sequences[f"sequence {seq_num}"] + line
        return sequences
