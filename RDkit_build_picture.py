
#!usr/bin/python3

# python sdftosmiles.py molecules.sdf

#conda activate my-rdkit-env
import os

import sys

from rdkit import Chem

from rdkit.Chem import Draw

from rdkit.Chem import AllChem

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

        Draw.MolToImageFile(m,save_name+".png",size=(300, 300))

        m = Chem.AddHs(m)

        AllChem.EmbedMolecule( m,randomSeed=3 )

        AllChem.MMFFOptimizeMolecule(m)

        Chem.MolToMolFile(m,file_name[:-9]+".mol")
        print(file_name[:-9])

    out_file.close()

 

if __name__=="__main__":


    folder = os.listdir("creat_mol/")
    count = 0
    folder.sort(key = lambda x : int(x[:-9]))
    for file_name in folder[:]:
        print(file_name)
        count+=1
        #try :
        converter("creat_mol/"+file_name,"picture/"+file_name[:-4])
        #except KeyboardInterrupt: Boost.Python.ArgumentError :
            #print("--------------")
        #finally :
            #print("------------")
            #continue
            #break

    
