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

def file_import(gen_id):
    url = (
        f"https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/{gen_id}/download?include_annotation_type=FASTA_RNA"
    )
    filename = f"static/{gen_id}.zip"
    urlretrieve(url, filename)

    with zf.ZipFile(f"static/{gen_id}.zip", "r") as seq:
        seq.extract("ncbi_dataset/data/rna.fna", path=f"sequences/{gen_id}")