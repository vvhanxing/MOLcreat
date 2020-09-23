#!/usr/bin/env python
# https://github.com/spro/char-rnn.pytorch
import random
import torch
import os
import argparse
import pickle
from helpers import *
from model import *
import read_mlo_71_try_class3_sort3

def generate(decoder, prime_str='2', predict_len=25,E=0.0 ,temperature=0.8, cuda=False):
    #E  =  random.choice([-300.0,-350.0,-400.0,-450.0,-500.0])
    hidden = decoder.init_hidden(1)+(  E   )/(25)
    #print("hidden,E",hidden.size(),E)
    prime_input = Variable(char_tensor(prime_str).unsqueeze(0))

    if cuda:
        hidden = hidden.cuda()
        prime_input = prime_input.cuda()
    predicted = prime_str

    # Use priming string to "build up" hidden state
    for p in range(len(prime_str) - 1):
        _, hidden = decoder(prime_input[:,p], hidden)
        
    inp = prime_input[:,-1]
    #print("prime_input",prime_input.size())
    #print("inp",inp.size(),inp)
    #print("predict_len",predict_len)
    
    for p in range(predict_len):
        #print("inp",inp)
        output, hidden = decoder(inp, hidden)
        
        #print("output",output.size())
        #print("hidden",hidden.size())
        #input()
        
        # Sample from the network as a multinomial distribution
        output_dist = output.data.view(-1).div(temperature).exp()
        top_i = torch.multinomial(output_dist, 1)[0]
        #print("top_i",top_i)

        # Add predicted character to string and use as next input
        predicted_char = all_characters[top_i]
        predicted += ' '+predicted_char
        #print("predicted_char----------",predicted_char)
        inp = Variable(char_tensor([predicted_char]).unsqueeze(0))
        #print("inp2",inp)
        #input()
        if cuda:
            inp = inp.cuda()

    return predicted



# Run as standalone script
if __name__ == '__main__':

# Parse command line arguments
    argparser = argparse.ArgumentParser()
    argparser.add_argument('-filename', type=str,default="ZINC250K_train.pt")
    argparser.add_argument('-p', '--prime_str', type=str, default='2')
    argparser.add_argument('-l', '--predict_len', type=int, default=100)
    argparser.add_argument('-t', '--temperature', type=float, default=0.8)
    argparser.add_argument('--cuda', action='store_true')
    args = argparser.parse_args()



    decoder = torch.load(args.filename,map_location=torch.device('cpu'))
    del args.filename 
    #print(generate(decoder, **vars(args)))

    #r = file_info_dict.get("3-1253-2347-162-0-0-0-0-0-0","Not fond")
    #input(r)

    print("------------------")

    
    count = 0
    folder = "creat/"
    #txt_created_list=[]
    txt_created_list=[]
    while count<=10000:
        #E  =  random.choice([0.2,0.25,0.3,0.35,0.4])
        txt = generate(decoder, "2",25,0.0,0.8,cuda=args.cuda).split(" 0 ")[0]
        count+=1
        print(count," ",txt)
        
        
        order_dir_name = "order_info_dir_NPname_3000_.plk"
        with open(order_dir_name,"rb") as f:
            order_dir = pickle.load(f)
        moltools = read_mlo_71_try_class3_sort3.MolTools()
        moltools.write_fake_mol("creat/"+str(count)+".mol",txt,order_dir)
            
        






    
