from urllib.request import urlretrieve
import zipfile as zf
import matplotlib.pyplot as plt

id = input("GenBank ID: ")
url = (
    f"https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{id}/download?include_annotation_type=FASTA_RNA"
)
filename = "static/test_seq.zip"
urlretrieve(url, filename)

with zf.ZipFile("static/test_seq.zip", "r") as seq:
    seq.extract("ncbi_dataset/data/rna.fna", path = "sequences")

with open("sequences/ncbi_dataset/data/rna.fna", "r") as rna:
    raw_seq = rna.read().splitlines()
    sequences = {}
    seq_num = 0
    for line in raw_seq:
        if line[0] == ">":
            seq_num += 1
            sequences[f"sequence {seq_num}"] = ""
        else:
            sequences[f"sequence {seq_num}"] = sequences[f"sequence {seq_num}"] + line

rw = {}
for seq, code in sequences.items():
    rw[seq]=[]
    tracker = 0
    for base in code:
        if base in ["A", "G"]:
            tracker +=1
        else:
            tracker -=1
        rw[seq].append(tracker)

plt.plot(rw["sequence 1"])
plt.show()


