
#!usr/bin/python3

# python sdftosmiles.py molecules.sdf

#conda activate my-rdkit-env
import pandas as pd
import matplotlib.pyplot as plt
from collections import OrderedDict
import numpy  as np
import os

import sys

from rdkit import Chem

from rdkit.Chem import Draw

from rdkit.Chem import AllChem

from rdkit.Chem import QED

from rdkit.Chem import Descriptors

import umap




def U_map(data):
    


    
    #fig, ax_array = plt.subplots(20, 20)
    #axes = ax_array.flatten()
    #for i, ax in enumerate(axes):
         #x.imshow(digits.images[i], cmap='gray_r')
    #plt.setp(axes, xticks=[], yticks=[], frame_on=False)
    #plt.tight_layout(h_pad=0.5, w_pad=0.01)
    #plt.show()


    
    reducer = umap.UMAP(random_state=42)
    embedding = reducer.fit_transform(data)
    #print(digits.data.shape)
    print(embedding.shape)

    plt.scatter(embedding[:, 0], embedding[:, 1], c=np.ones(embedding.shape[0]), cmap='Spectral', s=5)
    plt.gca().set_aspect('equal', 'datalim')
    plt.colorbar(boundaries=np.arange(11)-0.5).set_ticks(np.arange(10))
    plt.title('UMAP projection of the Digits dataset')
    plt.show()














def converter(file_name,save_name):

    mols = [ mol for mol in Chem.SDMolSupplier( file_name ) ]

    outname = save_name + ".smi"

    out_file = open( outname, "w" )

    for  mol in mols:

        smi = Chem.MolToSmiles(mol)
        #print(smi)

        name = mol.GetProp("_Name")

        out_file.write( "{}\t{}\n".format(smi, name ))

        m = Chem.MolFromSmiles(smi)

        m_qed = Chem.QED.qed(m)
        
        m_LogP = round(Descriptors.MolLogP(mol), 4)

        print(file_name,end = "  ")
        
        print("->",m_qed,m_LogP) #Chem.QED.properties(m)

        Draw.MolToImageFile(m,save_name+".png",size=(300, 300))

        m = Chem.AddHs(m)

        AllChem.EmbedMolecule( m,randomSeed=3 )

        try :
            AllChem.MMFFOptimizeMolecule(m)

            #Chem.MolToMolFile(m,file_name+".mol")

            #out_file.close()

            return smi,m_qed,m_LogP
            
        except ValueError:
            
            print("Rdkit not opt mol")
            return 0
            

    

 

if __name__=="__main__":


    folder = os.listdir("creat/")
    count = 0
    #folder.sort(key = lambda x : int(x[:-9]))
    dir_order = {}
    
    smi_list = []
    m_qed_list = []
    m_LogP_list = []
    m_file_name_list = []
    
    
    

            
    for file_name in folder[:20]:
        print(file_name)

                          
        
        count+=1

        info = converter("creat/"+file_name,"picture/"+file_name[:-4])
        if info !=0:
            smi,m_qed,m_LogP = info

            smi_list.append(smi)
            m_qed_list.append(m_qed)
            m_LogP_list.append(m_LogP )
            m_file_name_list.append(file_name)
        





    dir_order ["Smi"] = smi_list

    mols=list(map(lambda x: Chem.MolFromSmiles(x), smi_list))
    fingerprint=np.array(list(map(lambda x: AllChem.GetMorganFingerprintAsBitVect(x, 2, 2048), mols)))
    print(fingerprint.shape)
    U_map(fingerprint)
    input()
    
    dir_order ["LogP"]=m_LogP_list
    
    dir_order ["qed"]= m_qed_list 

    dir_order ["Mol name"] = m_file_name_list
    
        
    dataframe = pd.DataFrame(dir_order)
    dataframe.to_csv("creatMol.csv")


    
