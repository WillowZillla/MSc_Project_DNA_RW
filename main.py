def main():
    # print("Welcome to the DNA Random Walk converter, please input the file name you would like analysed")
    # file_name = input()
    with open("static/genbank.txt", 'r') as f:
        seq_file = f.read()
        print(seq_file)
