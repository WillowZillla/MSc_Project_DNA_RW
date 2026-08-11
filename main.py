import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from funcs import *
import hurst
import scipy.signal as sp
import csv
from fathon


while True:
    try:
        x = int(input("Please select which service you would like to use:\n"
                      "1. View the time series graphs and analysis for one gene\n"
                      "2. View the time series analysis for a group of genes "
                      "(please upload a list of RefSeq Gene IDs in ID.txt in the working directory)\n"
                      "3. Exit\n"))
    except ValueError:
        print("Invalid selection, please input either 1 or 2")
        continue
    if x in [1, 2, 3]:
        match x:
            case 1:
                gen_id = input("GenBank ID: ")
                file_import(gen_id)

                sequence = seq_extract(gen_id)

                one_d_series = one_dimension(sequence)
                two_d_series = two_dimension(sequence)
                three_d_series = three_dimension(sequence)
                for mapping, series in one_d_series.items():
                    H, c, data = hurst.compute_Hc(series = series, kind = "random_walk", simplified = True)
                    print(f"{mapping} Hurst Exponent: {H}")
                    if H>0.55:
                        print("Persistent trend")
                    elif H<0.45:
                        print("Mean reversal")
                    else:
                        print("Random Walk")
                    k = katz(series, 1)
                    print(f"{mapping} Katz dimension: {k}")


                    fig, (ax0, ax1) = plt.subplots(2, 1, layout="constrained")
                    plt.title = gen_id+" "+mapping
                    ax0.set_xlabel("base number")
                    ax0.set_ylabel("static time series")
                    ax0.plot(series)
                    ax0.set_xlabel("Frequency")
                    ax0.set_ylabel("V2/Hz")
                    ax1.psd(series, NFFT=16384)   # get all graphs on one image to show it
                    plt.show()


                for mapping, series in two_d_series.items():
                    plt.title = mapping
                    plt.xlabel(mapping[0]+mapping[1])
                    plt.ylabel(mapping[3]+mapping[4])
                    plt.plot(series[0], series[1])      # get all graphs on one image to show it
                    plt.show()
                fig = plt.figure()
                ax = plt.axes(projection = "3d")
                ax.plot3D(three_d_series[0], three_d_series[1], three_d_series[2])
                ax.set_title(gen_id)
                plt.show()


            case 2:
                try:
                    with open("ID.txt", "r") as f:
                        gen_ids = f.read()
                except FileNotFoundError:
                    print("No file detected, "
                          "please upload a comma separated list of gene IDs (e.g. 3208, 4288, 2764) "
                          "in a text file named ID.txt to the working directory")
                    continue

                gen_ids_list = gen_ids.split(", ")
                one_d = {}
                two_d = {}
                three_d = {}
                results = {}

                for gen_id in gen_ids_list:
                    file_import(gen_id)
                    sequence = seq_extract(gen_id)

                    # A, C, G, T = 0, 0, 0, 0
                    # for base in sequence:
                    #     match base:
                    #         case "A":
                    #             A+=1
                    #         case "G":
                    #             G+=1
                    #         case "C":
                    #             C+=1
                    #         case "T":
                    #             T+=1
                    # A = A / len(sequence)*100
                    # C = C / len(sequence) * 100
                    # T = T / len(sequence) * 100
                    # G = G / len(sequence) * 100
                    # print(f"Sequence Composition:\n"
                    #       f"\t* A: {round(A, 2)}%\n"
                    #       f"\t* C: {round(C, 2)}%\n"
                    #       f"\t* T: {round(T, 2)}%\n"
                    #       f"\t* G: {round(G, 2)}%\n"
                    #       f"Total: {A+G+C+T}%")

                    one_d[gen_id]= one_dimension(sequence)
                    two_d[gen_id] = two_dimension(sequence)
                    three_d[gen_id] = three_dimension(sequence)

                print("\n###############################\nStarting 1 dimensional analysis...\n###############################\n")

                for gen_id, mappings in one_d.items():
                    results[gen_id] = {"1D": {}}
                    for mapping, series in mappings.items():
                        H, c, data = hurst.compute_Hc(series=series, kind="random_walk", simplified=True)
                        K = katz(series, 1)
                        M = mean_pos(series, 1)
                        results[gen_id]["1D"][mapping]={f"Hurst Exponent": H,
                                                        f"Katz Dimension": K,
                                                        f"Mean Position": M}
                    print(f"{gen_id}\t\t{round(((gen_ids_list.index(gen_id)+1)/len(gen_ids_list))*100, 2)}%")

                print("1 dimensional analysis complete!")
                print("\n###############################\nStarting 2 dimensional analysis...\n###############################\n")

                for gen_id, mappings in two_d.items():
                    results[gen_id]["2D"] = {}
                    for mapping, series in mappings.items():
                        x_series = []
                        y_series = []
                        for i in range(len(series[0])):
                            x_series.append(series[0][i])
                            y_series.append(series[1][i])
                        xH, xc, xdata = hurst.compute_Hc(series=x_series, kind="random_walk", simplified=True)
                        yH, yc, ydata = hurst.compute_Hc(series=y_series, kind="random_walk", simplified=True)
                        K = katz(series, 2)
                        Mx, My = mean_pos(series, 2)
                        results[gen_id]["2D"][mapping] = {f"{mapping[0] + mapping[1]} Hurst Exponent": xH,
                                                          f"{mapping[3] + mapping[4]} Hurst Exponent": yH,
                                                          f"Katz Dimension": K,
                                                          f"Mean {mapping[0] + mapping[1]} Position": Mx,
                                                          f"Mean {mapping[3] + mapping[4]} Position": My}
                    print(f"{gen_id} \t\t{round(((gen_ids_list.index(gen_id)+1)/len(gen_ids_list))*100, 2)}%")

                print("2 dimensional analysis complete!")
                print("\n###############################\nStarting 3 dimensional analysis...\n###############################\n")

                for gen_id, series in three_d.items():
                    results[gen_id]["3D"] = {}
                    x_series = []
                    y_series = []
                    z_series = []
                    for i in range(len(series[0])):
                        x_series.append(series[0][i])
                        y_series.append(series[1][i])
                        z_series.append(series[2][i])
                    xH, xc, xdata = hurst.compute_Hc(series=x_series, kind="random_walk", simplified=True)
                    yH, yc, ydata = hurst.compute_Hc(series=y_series, kind="random_walk", simplified=True)
                    zH, zc, zdata = hurst.compute_Hc(series=x_series, kind="random_walk", simplified=True)
                    K = katz(series, 2)
                    Mx, My, Mz = mean_pos(series, 3)
                    results[gen_id]["3D"]["ACG"] = {f"Adenine Hurst Exponent": xH,
                                                    f"Guanine Hurst Exponent": yH,
                                                    f"Cytosine Hurst Exponent": zH,
                                                    f"Katz Dimension": K,
                                                    f"Mean A Position": Mx,
                                                    f"Mean G Position": My,
                                                    f"Mean C Position": Mz}
                    print(f"{gen_id}\t\t{round(((gen_ids_list.index(gen_id)+1)/len(gen_ids_list))*100, 2)}%")

                print("3 dimensional analysis complete!")
                print("Would you like to name your output directory? y/n")
                name_choice = input()
                if name_choice == "y" or name_choice == "Y":
                    name = input("Directory Name: ")
                else:
                    name = None

                print("Saving results...")

                file_data = {"1D": {"AC-GT": {}, "AG-CT": {}, "AT-CG": {}}, "2D": {"AC-GT": {}, "AG-CT": {}, "AT-CG": {}}, "3D": {"ACG": {}}}
                for gen_id, dimensions in results.items():
                    for dimension, mappings in dimensions.items():
                        for mapping, data in mappings.items():
                            file_data[dimension][mapping][gen_id] = {}
                            for label, datum in data.items():
                                file_data[dimension][mapping][gen_id][label] = datum

                iteration = len(os.listdir("results"))
                if name:
                    file_path = f"results/{name}"
                else:
                    file_path = f"results/output_{iteration}"
                i = 1
                while True:
                    try:
                        os.mkdir(file_path)
                        break
                    except FileExistsError:
                        if i > 1:
                            file_path = file_path[:-3]
                            file_path = file_path + f"({i})"
                        else:
                            file_path = file_path + f"({i})"
                        i += 1
                for dimension, mappings in file_data.items():
                    for mapping, gen_ids in mappings.items():
                        with open(f"{file_path}/{dimension}_{mapping}.csv", "w") as f:
                            first = True
                            w = csv.writer(f)
                            table = []
                            for gen_id, data in gen_ids.items():
                                if first:
                                    headers = ["Gene ID", "Gene"]
                                    for label in data.keys():
                                        headers.append(label)
                                    table.append(headers)
                                    first = False
                                row = [gen_id, get_name(gen_id)]                      # implement get_name() here
                                for datum in data.values():
                                    row.append(datum)
                                table.append(row)
                            w.writerows(table)

                print(f"\nYour results can be found in the working directory under <{file_path}> :)\n")


                    # for gen_id, dimensions in results.items():
                    #     print(gen_id + ":")
                    #     for dimension, mappings in dimensions.items():
                    #         print(f"* {dimension}:")
                    #         for mapping, exps in mappings.items():
                    #             print(f"\t* {mapping}:")
                    #             for exp, value in exps.items():
                    #                 print(f"\t\t-{exp} = {value}")
            case 3:
                break
    else:
        print("Invalid selection, please input either 1 or 2")