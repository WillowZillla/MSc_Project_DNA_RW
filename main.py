def main():
    print("Welcome to the DNA Random Walk converter, please input the file name you would like analysed")
    file_name = input()
    with open(file_name, 'r') as f:
        seq_file = f.read()
        print
