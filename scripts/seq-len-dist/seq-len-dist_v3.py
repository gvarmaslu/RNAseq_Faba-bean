"""
Fa sequence of file processing
 # Get id and sequence information sequence
 Id corresponding to the length of each sequence statistics #
 # Sequence length statistics
 # Draw histogram
"""
import os
import pandas as pd
from collections import Counter
import matplotlib.pyplot as plt
import sys

def read_seq(file_path):
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                seq_id.append(line.rstrip('\n').replace('>', ''))
            else:
                seq.append(line.rstrip('\n'))
        seq_dic = dict (zip (seq_id, seq)) # sequence number and the sequence information list merge into two dictionaries
    for v in seq_dic.values():
        seq_len.append(len(v))
    return seq_id , seq, seq_len


def write_count(seq_id, seq, seq_len):
    # Pandas write sequence number, sequence information, sequence length
    data1 = pd.DataFrame({"Seq_ID": seq_id})
    data2 = pd.DataFrame({"Seq_Info": seq})
    data3 = pd.DataFrame({"Seq_Len": seq_len})

    #Using the writer = pd.ExcelWriter (abs_path + '\\' + "test1.xlsx") # windows
    data1.to_excel(writer,sheet_name="data",startcol=0,index=False)
    data2.to_excel(writer,sheet_name="data",startcol=1,index=False)
    data3.to_excel(writer,sheet_name="data",startcol=2,index=False)
    #Save # writer.save () # data to excel file

def count_bar(seq_len):
    """
         # The sequence length information obtained by step, its sort / uniq, matplotlib processed and histogrammed
         First, the data # ordered statistics counts the same length
         Videos for data cleansing bar # view of
    """
        #len_count = column length data to obtain a third step for cleaning the Counter (seq_len) # extract and count the number of each length
        # Matplotlib drawing
    x = []
    y = []
    for k ,v in len_count.items():
        x.append(k)
        y.append(v)
    print(x ,y)
    plt.bar(x,y)
    plt.xlabel("Sequence Length")
    plt.ylabel("Sequence Numbers")
    plt.show()

#def main():


if __name__ == "__main__":
        abs_path = os.getcwd () # get the current directory path
        #print(abs_path)
        file_name = sys.argv[0]
        file_path = abs_path + '\\' + '* .fa' # obtain files in the current directory 
        seq_id = [] # Create fasta file list storing sequence number information
        seq = [] # Create list storage file corresponding to the sequence information fasta
        seq_len = [] # Create list storing information corresponding to the length of the sequence
        read_seq(file_path)
        write_count(seq_id, seq, seq_len)
        count_bar(seq_len)


