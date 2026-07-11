from urllib.request import urlretrieve
import zipfile as zf
import matplotlib.pyplot as plt
from funcs import one_dimension, file_import, seq_extract

while True:
    try:
        x = int(input("Please select which service you would like to use:\n"
                      "1. View the time series graphs and analysis for one gene\n"
                      "2. View the time series analysis for a group of genes (please upload a list of RefSeq Gene IDs in ID.txt in the working directory)\n"
                      "3. Exit\n"))
    except ValueError:
        print("Invalid selection, please input either 1 or 2")
        continue
    if x in [1, 2, 3]:
        match x:
            case 1:
                gen_id = input("GenBank ID: ")
                file_import(gen_id)

                sequences = seq_extract(gen_id)

                one_d_series = one_dimension(sequences)
                for seq_num, series in one_d_series.items():
                    plt.title = seq_num
                    plt.xlabel("base number")
                    plt.ylabel("static time series")
                    plt.plot(series)      # get all graphs on one image to show it
                    plt.show()

            case 2:
                print("triggered case 2")
                try:
                    with open("ID.txt", "r") as f:
                        gen_ids = f.read()
                        print("no error :) smiles so sneetly")
                except FileNotFoundError:
                    print("No file detected, please upload a comma separated list of gene IDs (e.g. <3208, 4288, 2764>) in a text file named ID.txt to the working directory")
                    continue
                gen_ids_list = gen_ids.split(", ")
                one_d_master = []
                two_d_master = []
                three_d_master = []
                for gen_id in gen_ids_list:
                    file_import(gen_id)
                    sequences = seq_extract(gen_id)
                    one_d_master.append(one_dimension(sequences))



                    # Don't want to have this many graphs show up, especially if its like 200 sequences

                    # for seq_num, series in one_d_series.items():
                    #     plt.title = seq_num
                    #     plt.xlabel("base number")
                    #     plt.ylabel("static time series")
                    #     plt.plot(series)  # get all graphs on one image to show it
                    #     plt.show()
            case 3:
                break
    else:
        print("Invalid selection, please input either 1 or 2")