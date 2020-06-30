#!/usr/bin/env python
# https://github.com/spro/char-rnn.pytorch

import torch
import torch.nn as nn
from torch.autograd import Variable
import argparse
import os

from tqdm import tqdm

from helpers import *
from model import *
from generate import *

# Parse command line arguments
argparser = argparse.ArgumentParser()
argparser.add_argument('--filename', type=str,default="mol_30000.txt")
argparser.add_argument('--model', type=str, default="gru")
argparser.add_argument('--n_epochs', type=int, default=6000)#2000
argparser.add_argument('--print_every', type=int, default=100)
argparser.add_argument('--hidden_size', type=int, default=200)
argparser.add_argument('--n_layers', type=int, default=2)
argparser.add_argument('--learning_rate', type=float, default=0.01)
argparser.add_argument('--chunk_len', type=int, default=10)
argparser.add_argument('--batch_size', type=int, default=200)
argparser.add_argument('--shuffle', action='store_true')
argparser.add_argument('--cuda', action='store_true')
args = argparser.parse_args()

if args.cuda:
    print("Using CUDA")

file, file_len,file_info_dict = read_file(args.filename)
#print(file_len /chunk_len)
#input()
def random_training_set(chunk_len, batch_size):
    inp = torch.LongTensor(batch_size, chunk_len)
    target = torch.LongTensor(batch_size, chunk_len)
    info = torch.FloatTensor(batch_size, 1)
    for bi in range(batch_size):
        start_index = random.randint(0, int(file_len /chunk_len)-2)*chunk_len
        end_index = start_index + chunk_len + 1
        #print("start_index,end_index",start_index,end_index)
        chunk = file[start_index:end_index]

        #i =
        i = chunk[:-1]
        try :
            inp[bi] = char_tensor(i)
        except RuntimeError:
            print(chunk[:-1])
        #print(inp[bi])

        t = chunk[1:-1]
        t.append('0')
        target[bi] = char_tensor(t)
        #print(t)

        info[bi]=(file_info_dict["-".join(i)]-0.325)/(0.65-0.325)

        
        #input()

    #print("inp,target",inp[0],target[0])
    #input()
    
    inp = Variable(inp)
    target = Variable(target)
    if args.cuda:
        inp = inp.cuda()
        target = target.cuda()
    #print(inp)
    #print(target)
    #print(info)
    #input()
    return inp, target,info

def train(inp, target,info):
    hidden = decoder.init_hidden(args.batch_size)
    hidden = hidden + info
    #print(hidden)
    #print(hidden.size())
    #print(info.size())
    #input()
    if args.cuda:
        hidden = hidden.cuda()
    decoder.zero_grad()
    loss = 0

    for c in range(args.chunk_len):
        #print("->",inp[:,c],inp[:,c].size())
        
        output, hidden = decoder(inp[:,c], hidden)

        loss += criterion(output.view(args.batch_size, -1), target[:,c])

        #print("inp[:,c]",inp[:,c],inp[:,c].size())
        #print("output.view(args.batch_size, -1)",output.view(args.batch_size, -1).size())
        #print("target[:,c]",target[:,c],target[:,c].size())
        #input()

    loss.backward()
    decoder_optimizer.step()

    return loss.item() / args.chunk_len

def save():
    save_filename = os.path.splitext(os.path.basename(args.filename))[0] + '.pt'
    torch.save(decoder, save_filename)
    print('Saved as %s' % save_filename)

# Initialize models and start training

decoder = CharRNN(
    n_characters,
    args.hidden_size,
    n_characters,
    model=args.model,
    n_layers=args.n_layers,
)
decoder_optimizer = torch.optim.Adam(decoder.parameters(), lr=args.learning_rate)
criterion = nn.CrossEntropyLoss()

if args.cuda:
    decoder.cuda()

start = time.time()
all_losses = []
loss_avg = 0

try:
    print("Training for %d epochs..." % args.n_epochs)
    for epoch in tqdm(range(1, args.n_epochs + 1)):
        loss = train(*random_training_set(args.chunk_len, args.batch_size))
        loss_avg += loss

        if epoch % args.print_every == 0:
            print('[%s (%d %d%%) %.4f]' % (time_since(start), epoch, epoch / args.n_epochs * 100, loss))
            E  =  random.choice([0.2,0.25,0.3,0.35,0.4])
            
            print(generate(decoder, '1', 10,E, 0.9,cuda=args.cuda), '\n')

    print("Saving...")
    save()

except KeyboardInterrupt:
    print("Saving before quit...")
    save()

