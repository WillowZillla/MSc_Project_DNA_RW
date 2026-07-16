import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from funcs import *
import hurst


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
                two_d_series = two_dimension(sequences)
                three_d_series = three_dimension(sequences)
                for seq_num, series in one_d_series.items():
                    H, c, data = hurst.compute_Hc(series = series, kind = "random_walk", simplified = True)
                    print(f"{seq_num} Hurst Exponent: {H}")
                    if H>0.55:
                        print("Persistent trend")
                    elif H<0.45:
                        print("Mean reversal")
                    else:
                        print("Random Walk")
                    plt.title = seq_num
                    plt.xlabel("base number")
                    plt.ylabel("static time series")
                    plt.plot(series)      # get all graphs on one image to show it
                    plt.show()
                for seq_num, series in two_d_series.items():
                    plt.title = seq_num
                    plt.xlabel("A's and G's")
                    plt.ylabel("T's and C's")
                    plt.plot(series[0], series[1])      # get all graphs on one image to show it
                    plt.show()
                for seq_num, series in three_d_series.items():
                    fig = plt.figure()
                    ax = plt.axes(projection = "3d")
                    ax.plot3D(series[0], series[1], series[2])
                    ax.set_title(seq_num)
                    plt.show()


            case 2:
                try:
                    with open("ID.txt", "r") as f:
                        gen_ids = f.read()
                except FileNotFoundError:
                    print("No file detected, please upload a comma separated list of gene IDs (e.g. 3208, 4288, 2764) in a text file named ID.txt to the working directory")
                    continue
                gen_ids_list = gen_ids.split(", ")
                one_d_master = {}
                two_d_master = {}
                three_d_master = {}
                hursts = {}
                for gen_id in gen_ids_list:
                    file_import(gen_id)
                    sequences = seq_extract(gen_id)
                    one_d_master[gen_id]= one_dimension(sequences)
                    two_d_master[gen_id] = two_dimension(sequences)
                    three_d_master[gen_id] = three_dimension(sequences)
                for gen_id, sequences in one_d_master.items():
                    for seq_num, series in sequences.items():
                        H, c, data = hurst.compute_Hc(series=series, kind="random_walk", simplified=True)
                        hursts[gen_id][seq_num] = H

                for gen_id, sequences in hursts.items():
                    print(gen_id + ":")
                    for seq_num, h in sequences:
                        print(f"{seq_num}, {h}")

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