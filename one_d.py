
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
