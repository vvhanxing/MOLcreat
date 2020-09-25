
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









def U_map(smis,color):

    



    #fig, ax_array = plt.subplots(20, 20)

    #axes = ax_array.flatten()

    #for i, ax in enumerate(axes):

         #x.imshow(digits.images[i], cmap='gray_r')

    #plt.setp(axes, xticks=[], yticks=[], frame_on=False)

    #plt.tight_layout(h_pad=0.5, w_pad=0.01)

    #plt.show()



    mols=list(map(lambda x: Chem.MolFromSmiles(x), smis))    

    fingerprint=np.array(list(map(lambda x: AllChem.GetMorganFingerprintAsBitVect(x, 2, 2048), mols)))



    print(fingerprint.shape)

    reducer = umap.UMAP(random_state=42)

    embedding = reducer.fit_transform(fingerprint)

    #print(digits.data.shape)

    print(embedding.shape)

    print()



    plt.scatter(embedding[:, 0], embedding[:, 1], c=np.array(color), cmap='Spectral', s=5)

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

        

        m_logP = round(Descriptors.MolLogP(mol), 4) #油水系数

        m_Wt = round(Descriptors.MolWt(mol), 4) #分子质量

        m_TPSA  =  round(Descriptors.TPSA(mol), 4) #极性表面积

        m_RotatableBonds= Descriptors.NumRotatableBonds(mol)  #可旋转键数目

        m_HBD = AllChem.CalcNumLipinskiHBD(mol)  #氢键供体

        m_HBA = AllChem.CalcNumLipinskiHBA(mol)  #氢键受体

        m_NumRings =  AllChem.CalcNumRings(mol)    # 环的个数

        m_AromaticRings = AllChem.CalcNumAromaticRings(mol) #芳香环个数

        m_CSP3 = round(AllChem.CalcFractionCSP3(mol), 4)   #SP3杂化C原子书目
        
        m_HeavyAtoms = mol.GetNumHeavyAtoms() #重原子的数目

        print(file_name,end = "  ")
        print("->",m_qed,m_logP) #Chem.QED.properties(m)

        Draw.MolToImageFile(m,save_name+".png",size=(300, 300))

        m = Chem.AddHs(m)


        AllChem.EmbedMolecule( m,randomSeed=3 )

        try :

            #AllChem.MMFFOptimizeMolecule(m)



            #Chem.MolToMolFile(m,file_name+".mol")



            #out_file.close()



            return smi,m_qed,m_logP,m_Wt ,m_TPSA,m_RotatableBonds,m_HBD,m_HBA,m_NumRings,m_AromaticRings,m_CSP3,m_HeavyAtoms,file_name

            

        except ValueError:

            

            print("Rdkit not opt mol")

            return 0

            







def build_creat_csv():

    folder = os.listdir("creat/")

    count = 0

    #folder.sort(key = lambda x : int(x[:-9]))

    dir_order = {}


    smi_list,m_qed_list,m_logP_list,m_Wt_list ,m_TPSA_list,m_RotatableBonds_list,m_HBD_list,m_HBA_list,m_NumRings_list,m_AromaticRings_list,m_CSP3_list,m_HeavyAtoms_list,file_name_list = [],[],[],[],[],[],[],[],[],[],[],[],[]


    for file_name in folder[:]:

        print(file_name)

        count+=1



        info = converter("creat/"+file_name,"picture/"+file_name[:-4])

        if info !=0:

            smi,m_qed,m_logP,m_Wt ,m_TPSA,m_RotatableBonds,m_HBD,m_HBA,m_NumRings,m_AromaticRings,m_CSP3,m_HeavyAtoms,file_name  = info

            smi_list.append(smi)
            m_qed_list.append(m_qed)
            m_logP_list.append(m_logP)
            m_Wt_list.append(m_Wt)
            m_TPSA_list.append(m_TPSA)
            m_RotatableBonds_list.append(m_RotatableBonds)
            m_HBD_list.append(m_HBD)
            m_HBA_list.append(m_HBA)
            m_NumRings_list.append(m_NumRings)
            m_AromaticRings_list.append(m_AromaticRings)
            m_CSP3_list.append(m_CSP3)
            m_HeavyAtoms_list.append(m_HeavyAtoms)
            file_name_list .append(file_name)


    P = [smi_list,m_qed_list,m_logP_list,m_Wt_list ,m_TPSA_list,m_RotatableBonds_list,m_HBD_list,m_HBA_list,m_NumRings_list,m_AromaticRings_list,m_CSP3_list,m_HeavyAtoms_list,file_name_list ]
    for index, p in enumerate("smi,m_qed,m_logP,m_Wt,m_TPSA,m_RotatableBonds,m_HBD,m_HBA,m_NumRings,m_AromaticRings,m_CSP3,m_HeavyAtoms,file_name".split(",")):
        print(index, p)
        dir_order [p] = P[index]


    dataframe = pd.DataFrame(dir_order)

    

    dataframe.to_csv("creatMol.csv")

    

    return dir_order 



 



if __name__=="__main__":











    

    





    

    csv_data = pd.read_csv("250k_rndm_zinc_drugs_clean_3.csv")

    smi_train_list = list(csv_data["smiles"])

    smi_train_list  = [ smi[:-1] for smi in smi_train_list if "+" not in smi and "-" not in smi][:50000]
    

    color_train = [0 for x in  smi_train_list][:100]



    print(len(smi_train_list),len(color_train))





    
    
    d#ir_order  = build_creat_csv()

    #smi_list = dir_order ["smi"]

    #color = [1 for x in smi_list]

    #print(len(smi_list ),len(color))



    csv_data = pd.read_csv("creatMol.csv")

    smi_list = list(csv_data["smi"])[:10000]  

    color = [1 for x in smi_list][:10000]



    



    print()





    

    smi_list_cat = smi_list#生成，标记为1

    smi_list_cat.extend(smi_train_list)





    

    color_cat  = color

    color_cat.extend(color_train)







    print(smi_list_cat[:10], color_cat[:10])

    

    print(len(smi_list_cat),len(color_cat))

    U_map(smi_list_cat, color_cat )    

    



    
