# https://github.com/spro/char-rnn.pytorch
import re
import unidecode
import string
import random
import time
import math
import torch

# Reading and un-unicode-encoding data

#all_characters = string.printable

all_characters=[str(x) for x in list(range(4000))  ]     

n_characters = len(all_characters)

def read_file(filename):
    file = unidecode.unidecode(open(filename).read())

    pattern = re.compile(r"\S+")
    file = pattern.findall(file)

    #print("read_file",file[:100])
    file_info_dict = {}
    with open(filename,"r") as info:
        for line in info.readlines():
            
            file_info_dict["-".join(line.split(" ")[:-1])+"-0"]=0.0
            #print("-".join(line.split(" ")[:-1])+" 0")
            #print("line",line,"-".join(line.split(" ")[:-2]),float(line.split(" ")[-1]))
            #input()

    return file, len(file),file_info_dict

# Turning a string into a tensor

def char_tensor(string):
    #print(string)

    word_list = string
    #print("word_list",word_list)

    #print("word_list ",word_list,len(word_list))
    
    tensor = torch.zeros(len(word_list)).long()
    
    for c in range(len(word_list)):
        try:
            tensor[c] = all_characters.index(word_list[c])
        except:
            continue
    return tensor

# Readable time elapsed

def time_since(since):
    s = time.time() - since
    m = math.floor(s / 60)
    s -= m * 60
    return '%dm %ds' % (m, s)

if __name__=="__main__":
    #print(all_characters)
    string='''93 717 40 15 13 16 109 245 53 10 14 0'''.split(" ")
    
    print()
    
    print(string)
    
    print(char_tensor(string))

    filename = "ZINC250K_train.txt"
    
    print(read_file(filename)[1])

   # print(char_tensor(" 13 "))
    
