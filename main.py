from urllib.request import urlretrieve
import zipfile as zf
import matplotlib.pyplot as plt
from one_d import one_dimension
import os

while True:
    x = input("Please select which service you would like to use:\n1. View the time series graphs and analysis for one gene\n2. View the time series analysis for a group of genes "
              "(please input a list og RefSeq Gene IDs in ID.txt in the working directory)")
    if x==1 or x==2:
        break
    else:
        print("Invalid selection, please input either 1 or 2")

match x:
    case 1:
        id = input("GenBank ID: ")
        url = (
            f"https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{id}/download?include_annotation_type=FASTA_RNA"
        )
        filename = f"static/{id}.zip"
        urlretrieve(url, filename)

        with zf.ZipFile(f"static/{id}.zip", "r") as seq:
            seq.extract("ncbi_dataset/data/rna.fna", path = f"sequences/{id}")

        with open(f"sequences/{id}/ncbi_dataset/data/rna.fna", "r") as rna:
            raw_seq = rna.read().splitlines()
            sequences = {}
            seq_num = 0
            for line in raw_seq:
                if line[0] == ">":
                    seq_num += 1
                    sequences[f"sequence {seq_num}"] = ""
                else:
                    sequences[f"sequence {seq_num}"] = sequences[f"sequence {seq_num}"] + line

        one_d_series = one_dimension(sequences)

        plt.plot(one_d_series)      # get all graphs on one image to show it
        plt.show()