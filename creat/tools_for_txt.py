import csv
import pickle

def get_nets_dir_info(nets_dir_name):
    

    with open(nets_dir_name,"rb") as f:

        nets_dir =pickle.load(f) #{ 矩阵名字：[矩阵，排序，出现次数] } #编码

    nets_dir_list = list(nets_dir.values() )

    for i in range(5000):
        print(nets_dir_list[i][1:]) 

#nets_dir_name = "nets_info_dir_NPname_3000_.plk"
#get_nets_dir_info(nets_dir_name)
#input()
    



        

def get_txt(txt_name,max_num ):

    with open(txt_name ,"r") as txt:
        all_count = 0
        count = 0
        lines = []

        l= []
        for line in txt.readlines():
            txt_list = line.split(" ")
            end_index = txt_list.index("0")
            mot_str = [x for x in txt_list[:end_index ]]
        
            if len( [x for x in mot_str if int(x)>max_num ])==0:
                
                count+=1
                
                l.append(len(mot_str))

                lines.append(" ".join(mot_str))

                    
            all_count+=1
            
        print(all_count,count,count/all_count,max(l),sum(l)/count)
        return lines
    

txt_name = "ZINC250K.txt"
get_txt(txt_name,4000)



def mol_str_pad(line):

        txt_list = line.split(" ")
        if "0" in txt_list:
            end_index = txt_list.index("0")
            mol_encode_list  = [str(x) for x in txt_list[:end_index ]]
            mol_str = " ".join(mol_encode_list )
            mol_str= mol_str + " 0"*len(mol_encode_list) +"\n"
        else:
            mol_encode_list  = [str(x) for x in txt_list]
            mol_str = " ".join(mol_encode_list )
            mol_str= mol_str + " 0"*(25-len(mol_encode_list)) +"\n"            

        return mol_str



def get_txt_for_train(txt_name):
    with open("ZINC250K_train.txt","a") as train_txt:

        for line in get_txt(txt_name,4000):
                
            train_txt.writelines(mol_str_pad(line))
            
txt_name = "ZINC250K.txt"
get_txt_for_train(txt_name)



print("finish")


input()
input()

def sort_txt():
    info_list = []
    with open("testQM9.txt","r") as txt_r:
        lines = txt_r.readlines()
        for line in lines:
            #print(  [  " ".join(line.split(" ")[:-2])   ,  line.split(" ")[-1]    ]  ,line.split(" ")[-1] [:-5] )
            info_list.append(   [  " ".join(line.split(" ")[:-2])   ,  line.split(" ")[-1]    ]   )
            #input()

    info_list.sort( key = lambda x:int(x[1][:-5])  )

    with open("testQM9_sort.txt","a") as txt_w:
        for line in info_list:
            
            txt_w.writelines(line[0]+" "+line[1])

#sort_txt()
#input("finish")
#input()




def cat_txtCode_info():
    with open("250k_smiles.csv","r") as csvfile:
        reader=csv.reader(csvfile)
        data = list(reader)[1:150003]
        LogP_data = [ x[1] for x in data]
        #lumo_data = [ float( x[4]) for x in data]
        print(LogP_data[-1])
        #input()



    with open("testQM9.txt","r") as txt_r:
        lines = txt_r.readlines()
        #lines.sort(key = lambda x :int(x.split(" ")[-1][:-5]))

        print(len(lines))
        print(len(LogP_data))

        with open("testQM9_LogP.txt","a") as txt_w:
            for line in list(zip(lines,LogP_data)):
                txt_w.writelines(  line[0][:-1] +" "+line[1]+"\n")
                
                
        


#cat_txtCode_info()














def bigger(line):
    #print(len(line.split(" ")))
    if line[:2]==" 0":
        #print(line,"-----------")
        return True    
    for i in  line.split(" ")[:-2]:

        try :
            if int(i)>3000:
            #print(i)
            #input()
                return True
        
        except ValueError:
            print(line.split(" ")[:-2])
            #input()

    return False

def creat_info_txt():
    
    with open("testQM9_LogP.txt","r") as txt_r:
        
    
        lines = txt_r.readlines()
        #str_count_list = []
        with open("mol_30000_info.txt","w") as txt_w:
            for line in lines:
                encode_str = line.split(" ")[:-2]
                encode_str = " ".join(encode_str)+" 0"*(26-len(encode_str))
                l =  encode_str   + " ".join(line.split(" ")[-2:] )
                #print(  l )
                
                if not bigger(l):
                    #print(l)
                    #input()
                    txt_w.writelines(l)
            #str_count_list.append(len(line.split(" ")[:-2]))
        
        #print(max(str_count_list ))
        #input("finish")
        #input()
        






            
creat_info_txt()
with open("mol_30000_info.txt","r") as txt_r:
    for i in txt_r.readlines():
        if i[:2] ==" 0":
            print(i)
        if i.split(" ")[-3]!="0":
            print(i)

                #print(line)


            #input()
print("finish")
#input()
def creat_txt():

    with open("mol_30000_info.txt","r") as txt_r:
        with open("mol_30000.txt","w") as txt_w:
            for i in txt_r.readlines():
                txt_w.writelines( " ".join(i.split(" ")[:-2])+"\n" )
    

creat_txt()
