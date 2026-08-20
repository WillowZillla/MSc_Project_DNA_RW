import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
from funcs import *
import hurst
import scipy.signal as sp
from scipy import stats
import csv
import kneed as kn
from fooof import FOOOF

while True:
    try:
        x = int(input("Please select which service you would like to use:\n"
                      "1. View the time series graphs and analysis for one gene\n"
                      "2. View the time series analysis for a group of genes "
                      "(please upload a list of RefSeq Gene IDs in ID.txt in the working directory)\n"
                      "3. View the time series analysis for a group of sequences "
                      "(please upload a list of FASTA formatted sequences in SEQ.txt in the working directory)\n"
                      "4. Exit\n"))
    except ValueError:
        print("Invalid selection, please input either 1 or 2")
        continue
    if x in [1, 2, 3, 4]:
        match x:
            case 1:
                gen_id = input("GenBank ID: ")
                file_import(gen_id)

                #sequence = seq_extract(gen_id)

                sequence = 'ATGGTGCATCTGACTCCTGAGGAGAAGTCTGCCGTTACTGCCCTGTGGGGCAAGGTGAACGTGGATGAAGTTGGTGGTGAGGCCCTGGGCAGGCTGCTGGTGGTCTACCCTTGGACCCAGAGGTTCTTTGAGTCCTTTGGGGATCTGTCCACTCCTGATGCTGTTATGGGCAACCCTAAGGTGAAGGCTCATGGCAAGAAAGTGCTCGGTGCCTTTAGTGATGGCCTGGCTCACCTGGACAACCTCAAGGGCACCTTTGCCACACTGAGTGAGCTGCACTGTGACAAGCTGCACGTGGATCCTGAGAACTTCAGGCTCCTGGGCAACGTGCTGGTCTGTGTGCTGGCCCATCACTTTGGCAAAGAATTCACCCCACCAGTGCAGGCTGCCTATCAGAAAGTGGTGGCTGGTGTGGCTAATGCCCTGGCCCACAAGTACCACTAA'

                one_d_series = one_dimension(sequence)
                two_d_series = two_dimension(sequence)
                three_d_series = three_dimension(sequence)
                for mapping, series in one_d_series.items():
                    H, c, data = hurst.compute_Hc(series = series, kind = "random_walk", simplified = True)
                    print(f"Sequence length: {len(series)}")
                    print(f"{mapping} Hurst Exponent: {H}")
                    if H>0.55:
                        print("Persistent trend")
                    elif H<0.45:
                        print("Mean reversal")
                    else:
                        print("Random Walk")
                    k = katz(series, 1)
                    print(f"{mapping} Katz dimension: {k}")
                    dfa_H = dfa_hurst(series)
                    print(f"{mapping} DFA exponent: {dfa_H}")
                    pxx, freq = plt.psd(series, NFFT=16384)
                    plt.close()
                    fig, (ax0, ax1) = plt.subplots(2, layout="constrained")
                    ax0.title = get_name(gen_id)+" "+mapping
                    ax0.set_xlabel("base position")
                    ax0.set_ylabel(f"{mapping[0]+mapping[1]} against {mapping[3]+mapping[4]}")
                    ax0.plot(series)
                    #ax1.title(f"{get_name(gen_id)} Power Spectrum")
                    #ax1.psd(series, NFFT=16384)   # get all graphs on one image to show it

                    #ax1.set_xlabel("Frequency (log scale)")
                    #ax1.set_ylabel("Amplitude (log scale)")
                    # ax2.psd(series, NFFT=8192)
                    #ax1.set_xscale("log")
                    ax1.title(f"{get_name(gen_id)} Power Spectrum")
                    ax1.set_xlabel("log(Frequency)")
                    ax1.set_ylabel("log(Amplitude)")
                    type1_freq = np.arange(0, 1+1/(len(pxx)), 1/(len(pxx)))
                    type1_freq = np.delete(type1_freq, 0)
                    log_type1_freq = np.log(type1_freq)
                    log_type1_pxx = np.log(pxx)
                    ax1.plot(log_type1_freq, log_type1_pxx)

                    # ax2.set_xlabel("log(Frequency)")
                    # ax2.set_ylabel("V2/Hz")
                    # type2_freq = np.delete(freq, 0)
                    # type2_pxx = np.delete(pxx, 0)
                    # log_type2_freq = np.log(type2_freq)
                    # log_type2_pxx = np.log(type2_pxx)
                    # ax2.plot(log_type2_freq, log_type2_pxx)

                    # low_knee = len(freq)//100
                    # freq_knee = freq[low_knee:]
                    # pxx_knee = np.log(pxx)[low_knee:]
                    #
                    # print(pxx_knee)
                    #
                    # raw_kn = kn.KneeLocator(freq_knee, pxx_knee, curve = "convex", direction = "decreasing")
                    # log_kn = kn.KneeLocator(np.log(freq_knee), pxx_knee, curve = "convex", direction = "decreasing")
                    #
                    # print(f"Unlogged knee:\n"
                    #       f" x = {raw_kn.knee}, y = {raw_kn.knee_y}")
                    # print(f"Logged knee:\n"
                    #       f" x = {log_kn.knee}, y = {log_kn.knee_y}\n"
                    #       f" log(x) = {np.log(raw_kn.knee)}, log(y) = {np.log(raw_kn.knee_y)}")

                    x = log_type1_freq[0:len(log_type1_freq)//10]
                    y = log_type1_pxx[0:len(log_type1_freq)//10]


                    slope, intercept, r, p, std_err = stats.linregress(x, y)
                    print(intercept, slope)
                    beta = []
                    for x in log_type1_freq:
                        beta.append(slope*x+intercept)
                    ax1.plot(log_type1_freq, beta, color = "red")
                    ax1.grid()

                    ax1.legend(["PSD"], [f"Linear Regression, slope = {slope}"], loc = "lower left")

                    # fm = FOOOF()
                    # report = fm.report(freq, pxx, [freq[0], freq[-1]])
                    # print(type(report))
                    # print(report)

                    plt.show()

                for mapping, series in two_d_series.items():
                    plt.title = get_name(gen_id)+" "+mapping
                    plt.xlabel(mapping[0]+"/"+mapping[1])
                    plt.ylabel(mapping[3]+"/"+mapping[4])
                    plt.plot(series[0], series[1])      # get all graphs on one image to show it
                    plt.show()
                fig = plt.figure()
                ax = plt.axes(projection = "3d")
                ax.plot3D(three_d_series[0], three_d_series[1], three_d_series[2])
                ax.set_title(gen_id)
                ax.set_xlabel("A")
                ax.set_ylabel("G")
                ax.set_zlabel("C")
                plt.title(get_name(gen_id)+" AGC")
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
                gen_id_remove = []
                for gen_id in gen_ids_list:
                    failed = 0
                    print(gen_id)
                    if gen_id in os.listdir("sequences"):
                        print(f"{gen_id} already downloaded")
                    else:
                        failed = file_import(gen_id)
                        print(failed)
                        if failed:
                            gen_id_remove.append(gen_id)
                            print(f"{gen_id} skipped")
                            continue
                        else:
                            print(f"{gen_id} downloaded successfully")

                    sequence = seq_extract(gen_id)
                    for base in sequence:
                        if base not in ["A", "C", "G", "T"]:
                            print(f"Sequence error in {get_name(gen_id)} (gene ID {gen_id}), sequence skipped")
                            gen_id_remove.append(gen_id)
                            continue
                    one_d[gen_id] = one_dimension(sequence)
                    two_d[gen_id] = two_dimension(sequence)
                    three_d[gen_id] = three_dimension(sequence)
                    print(f"{gen_id} sequence extracted")
                for x in gen_id_remove:
                    gen_ids_list.remove(x)
                print(gen_ids_list)

                seq_analysis(gen_ids_list, one_d, two_d, three_d)


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

                # opt_n, opt_clusters, silhouette_score = k_means(results)
                # print(f"A total of {opt_n} clusters were found with a silhouette score of {silhouette_score}:")
                # print(f"opt_clusters = {opt_clusters}")

                    # for gen_id, dimensions in results.items():
                    #     print(gen_id + ":")
                    #     for dimension, mappings in dimensions.items():
                    #         print(f"* {dimension}:")
                    #         for mapping, exps in mappings.items():
                    #             print(f"\t* {mapping}:")
                    #             for exp, value in exps.items():
                    #                 print(f"\t\t-{exp} = {value}")
            case 3:
                try:
                    with open("SEQ.txt", "r") as f:
                        sequences = f.readlines()
                except FileNotFoundError:
                    print("No file detected, "
                          "please upload a list of comma separated gene names and sequences (<gene name>: <gene sequence>, etc.) "
                          "in a text file named SEQ.txt to the working directory")
                    continue
                sequences_list = {}
                seq_names = []
                seq_remove = []
                for sequence in sequences:
                    seqsplit = sequence.strip().split(": ")
                    for base in seqsplit[1]:
                        if base not in ["A", "C", "G", "T"]:
                            print(f"Sequence error in {seqsplit[0]}, sequence skipped")
                            seq_remove.append(seqsplit[0])
                    sequences_list[seqsplit[0]] = seqsplit[1]
                    seq_names.append(seqsplit[0])
                for name in seq_remove:
                    del sequences_list[name]
                one_d = {}
                two_d = {}
                three_d = {}
                results = {}
                for name, sequence in sequences_list.items():
                    one_d[name] = one_dimension(sequence)
                    two_d[name] = two_dimension(sequence)
                    three_d[name] = three_dimension(sequence)

                print(
                    "\n###############################\nStarting 1 dimensional analysis...\n###############################\n")

                for gen_id, mappings in one_d.items():
                    results[gen_id] = {"1D": {}}
                    for mapping, series in mappings.items():
                        H, c, data = hurst.compute_Hc(series=series, kind="random_walk", simplified=True)
                        K = katz(series, 1)
                        M = mean_pos(series, 1)
                        dfa_H = dfa_hurst(series)
                        results[gen_id]["1D"][mapping] = {f"R/S Hurst Exponent": H,
                                                          f"DFA Hurst Exponent": dfa_H,
                                                          f"Katz Dimension": K,
                                                          f"Mean Position": M}
                    print(f"{gen_id}\t\t{round(((seq_names.index(gen_id) + 1) / len(seq_names)) * 100, 2)}%")

                print("1 dimensional analysis complete!")
                print(
                    "\n###############################\nStarting 2 dimensional analysis...\n###############################\n")

                for gen_id, mappings in two_d.items():
                    results[gen_id]["2D"] = {}
                    for mapping, series in mappings.items():
                        x_series = []
                        y_series = []
                        for i in range(len(series[0])):
                            x_series.append(series[0][i])
                            y_series.append(series[1][i])
                        Hx, xc, xdata = hurst.compute_Hc(series=x_series, kind="random_walk", simplified=True)
                        Hy, yc, ydata = hurst.compute_Hc(series=y_series, kind="random_walk", simplified=True)
                        dfa_Hx = dfa_hurst(x_series)
                        dfa_Hy = dfa_hurst(y_series)
                        K = katz(series, 2)
                        Mx, My = mean_pos(series, 2)
                        results[gen_id]["2D"][mapping] = {f"{mapping[0] + mapping[1]} R/S Hurst Exponent": Hx,
                                                          f"{mapping[3] + mapping[4]} R/S Hurst Exponent": Hy,
                                                          f"{mapping[0] + mapping[1]} DFA Hurst Exponent": dfa_Hx,
                                                          f"{mapping[3] + mapping[4]} DFA Hurst Exponent": dfa_Hy,
                                                          f"Katz Dimension": K,
                                                          f"Mean {mapping[0] + mapping[1]} Position": Mx,
                                                          f"Mean {mapping[3] + mapping[4]} Position": My}
                    print(f"{gen_id} \t\t{round(((seq_names.index(gen_id) + 1) / len(seq_names)) * 100, 2)}%")

                print("2 dimensional analysis complete!")
                print(
                    "\n###############################\nStarting 3 dimensional analysis...\n###############################\n")

                for gen_id, series in three_d.items():
                    results[gen_id]["3D"] = {}
                    x_series = []
                    y_series = []
                    z_series = []
                    for i in range(len(series[0])):
                        x_series.append(series[0][i])
                        y_series.append(series[1][i])
                        z_series.append(series[2][i])
                    Hx, xc, xdata = hurst.compute_Hc(series=x_series, kind="random_walk", simplified=True)
                    Hy, yc, ydata = hurst.compute_Hc(series=y_series, kind="random_walk", simplified=True)
                    Hz, zc, zdata = hurst.compute_Hc(series=x_series, kind="random_walk", simplified=True)
                    dfa_Hx = dfa_hurst(x_series)
                    dfa_Hy = dfa_hurst(y_series)
                    dfa_Hz = dfa_hurst(z_series)
                    K = katz(series, 2)
                    Mx, My, Mz = mean_pos(series, 3)
                    results[gen_id]["3D"]["AGC"] = {f"Adenine R/S Hurst Exponent": Hx,
                                                    f"Guanine R/S Hurst Exponent": Hy,
                                                    f"Cytosine R/S Hurst Exponent": Hz,
                                                    f"Adenine DFA Hurst Exponent": dfa_Hx,
                                                    f"Guanine DFA Hurst Exponent": dfa_Hy,
                                                    f"Cytosine DFA Hurst Exponent": dfa_Hz,
                                                    f"Katz Dimension": K,
                                                    f"Mean A Position": Mx,
                                                    f"Mean G Position": My,
                                                    f"Mean C Position": Mz}
                    print(f"{gen_id}\t\t{round(((seq_names.index(gen_id) + 1) / len(seq_names)) * 100, 2)}%")

                print("3 dimensional analysis complete!")
                print("Would you like to name your output directory? y/n")
                name_choice = input()
                if name_choice == "y" or name_choice == "Y":
                    name = input("Directory Name: ")
                else:
                    name = None

                print("Saving results...")

                file_data = {"1D": {"AC-GT": {}, "AG-CT": {}, "AT-CG": {}},
                             "2D": {"AC-GT": {}, "AG-CT": {}, "AT-CG": {}},
                             "3D": {"AGC": {}}}
                for gen_id, dimensions in results.items():
                    for dimension, mappings in dimensions.items():
                        for mapping, data in mappings.items():
                            file_data[dimension][mapping][gen_id] = {}
                            for label, datum in data.items():
                                file_data[dimension][mapping][gen_id][label] = datum

                iteration = len(os.listdir("results"))
                results_json = json.dumps(results, indent=4)
                with open(f"results_{iteration}.json", "w") as r:
                    r.write(results_json)
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
                                    headers = ["Gene Name"]
                                    for label in data.keys():
                                        headers.append(label)
                                    table.append(headers)
                                    first = False
                                row = [gen_id]
                                for datum in data.values():
                                    row.append(datum)
                                table.append(row)
                            w.writerows(table)
                k_means_params = input("Which parameters would you like to use for the k-means analysis?\n"
                                       "1. All parameters\n"
                                       "2. 1D parameters\n"
                                       "3. 2D parameters\n"
                                       "4. 3D parameters\n"
                                       "5. 3D mean position\n")
                opt_n, opt_clusters, max_silhouette, max_silhouette_samples = k_means(results, k_means_params)
                k_means_results = {}
                for i in range(len(results.keys())):
                    k_means_results[list(results.keys())[i]] = [opt_clusters[i], max_silhouette_samples[i]]
                with open(f"{file_path}/k_means.csv", "w") as f:
                    w = csv.writer(f)
                    headers = ["Gene Name", "Cluster", "Silhouette Score"]
                    table = [headers]
                    for name, clusters in k_means_results.items():
                        table.append([name, clusters[0], clusters[1]])
                    w.writerows(table)

                print(f"\nYour results can be found in the working directory under <{file_path}> :)\n"
                      f"The optimum number of clusters was {opt_n}, with an average silhouette score of {max_silhouette}")



            case 4:
                break
    else:
        print("Invalid selection, please input either 1, 2, 3 or 4")
