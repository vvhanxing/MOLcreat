import random
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole 
from rdkit.Chem.Draw.MolDrawing import MolDrawing, DrawingOptions
from rdkit.Chem import Descriptors

def randow_replace():

    replace_atom = [ "H","C","N","O","S" ,"F","Cl","Br","I"]   
    smi_unknow = "[*]c1c([*])c2c([*](c(C(NC)=O)c2[*])[*])c([*])c1[*]"

    smi = ""

    for fragment in smi_unknow.split("*"):

        selcet_atom =  random.choice(replace_atom) 
        smi +=fragment+selcet_atom


    print(smi )
    return smi



randow_replace()

#mol = Chem.MolFromSmiles(randow_replace())
#'CC1=C(C(C)=C(C2=C1[H])N(C(C(NC)=O)=C2Br)C)N'
#SMI  = "[H]C(C(CC)=C([H])C1=C2N(C(C(NC)=O)=C1Br)CN)=C2[H]"
#SMI = "[H]C(C(CN)=C([H])C1=C2N(C(C(NC)=O)=C1Br)CC)=C2[H]"
#smi = "[H]C(C(CN)=C([H])C1=C2N(C(C(NC)=O)=C1[H])CCBr)=C2[H]"
#smi = "[H]C(C(CNBr)=C([H])C1=C2N(C(C(NC)=O)=C1[H])CC)=C2[H]"
#smi = "[H]C(C(CCBr)=C([H])C1=C2N(C(C(NC)=O)=C1[H])CN)=C2[H]"
#smi = "[H]C(C([H])=C1[H])=C(CCBr)C2=C1C([H])=C(C(NC)=O)N2CN"
#smi = "[H]C(C([H])=C1[H])=CC2=C1C(CCBr)=C(C(NC)=O)N2CN"
#smi = "[H]C(C([H])=C1[H])=CC2=C1C(C(Br)C)=C(C(NC)=O)N2CN"
#smi = "[H]C(C([H])=C1[H])=C(Br)C2=C1C(CC)=C(C(NC)=O)N2CN"
#smi = "[H]C(C([H])=C1[H])=C(Br)C2=C1C(CN)=C(C(NC)=O)N2CC"
#smi = "[H]C(C([H])=C1[H])=CC2=C1C(CBr)=C(C(NC)=O)N2CCN"
#smi = "[H]C(C([H])=C1[H])=CC2=C1C(CBr)=C(C(NC)=O)N2CCN"
#smi = "[H]C(C([H])=C1)=CC2=C1C(C(C)CBr)=C(C(NC)=O)N2N"
#smi = "[H]C(C([H])=C1)=CC2=C1C(CCBr)=C(C(NC)=O)N2CN"
#smi = "[H]C(C(Br)=C1)=CC2=C1C(CC)=C(C(NC)=O)N2CN"
#smi = "[H]C(C([H])=C1)=C(C(Br)C)C2=C1C=C(C(NC)=O)N2CN"
#smi = "[H]C(C([H])=C1)=C(CN)C2=C1C=C(C(NC)=O)N2CCBr"
with open("smi.txt","r") as txt:
    for smi in    list(txt.readlines()):
        
        #print(smi)
        
        mol = Chem.MolFromSmiles(smi )
            


        a = round(Descriptors.MolWt(mol), 2)
        b = round(Descriptors.MolLogP(mol), 2)
        c = round(Descriptors.TPSA(mol), 2)
        d = Descriptors.NumRotatableBonds(mol)
        e = AllChem.CalcNumLipinskiHBD(mol)
        f = AllChem.CalcNumLipinskiHBA(mol)
        g = AllChem.CalcNumRings(mol)
        h = AllChem.CalcNumAromaticRings(mol)
        i = round(AllChem.CalcFractionCSP3(mol), 2)
        j = mol.GetNumHeavyAtoms()

        print(a,b,c,d,e,f,g,h,i,j)
print("310.2 2.11 60.05 4 0~3 4~7 2 2 0.31 17~20")
