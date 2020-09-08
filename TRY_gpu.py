#import os
#os.environ['CUDA_VISIBLE_DEVICES'] = "0"

#def do_some_thing()
    #for i in range(1000):
        #print(i)


#conda install --channel https://conda.anaconda.org/conda-forge

import pandas as pd

from rdkit import Chem

from moses.metrics import get_all_metrics


gen = ['CNC', 'CCC',

                    'CC',

                    'CCC',

                    'CCCC']

print(get_all_metrics(gen))
input()
