from rdkit import Chem
import csv

def smi2mol(smi,isAddHs=True):
    """
    smi to sdf txt by rdkit 
    """
    mol = Chem.MolFromSmiles(smi)

    if isAddHs:
        
        mol = Chem.AddHs(mol)

    return Chem.MolToMolBlock(mol)


def writeSdf(sdfTxt,name):
    
    with open(name,"w") as txt:
        
        txt.writelines(sdfTxt)
        
    return True

if __name__ =="__main__":
    
    folder  = "mols/"
    
    csv_file_name = "250k_rndm_zinc_drugs_clean_3.csv"
    
    with open(csv_file_name) as f:
        
        reader = csv.reader(f)
        
        rows = [row for row in  reader]

        for row_count,row in enumerate(rows):
            
            if row_count>0:
                
                smi = row[0]

                if "+" not in smi and "-" not in  smi:
                
                    name =  "0"*(6-len(str(row_count)))+str(row_count)+".mol"
                
                    writeSdf(smi2mol(smi),folder+name)
                
                    print(name)
            
        
