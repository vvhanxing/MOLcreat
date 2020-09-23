
###################
#                 #
# 6.6生成错误处理 #
#                 #
###################

import os

import time

import re

from collections import Counter

from copy import deepcopy

import pickle

import random

import numpy as np

import matplotlib.pyplot as plt

from matplotlib import colors

from collections import OrderedDict

#from mpl_toolkits.axes_grid1 import make_axes_locatable



BANDencode = OrderedDict( [

    ("0",0.0),

    ("1",0.3),

    ("2",0.5),

    ("3",0.6),

    ("4",0.35),

    ("11",0.1),

    ("12",0.15),

    ("13",0.2),

    ("14",0.25)

    ])



ATOM_type_list =[

    (0.11,"h1"),    #数值大小也决定它在同层中的位置

    (0.12,"h2"),    

    (0.13,"h3"),    

    (0.14,"h4"),   

    (0.1,"H"),

    (0.2,"C"),

    (0.3,"N"),

    (0.4,"O"),

    (0.5,"F"),

    (0.6,"P"),

    (0.7,"S"),

    (0.8,"Cl"),

    (0.9,"Br"),

    (1.0,"I"),

    (1.1,"Other")]



ATOM_valence =OrderedDict([

    ("H",(1,)),

    ("C",(4,)),

    ("N",(3,)),

    ("O",(2,)),

    ("F",(1,)),

    ("P",(3,5)),

    ("S",(2,6)),

    ("Cl",(1,)),

    ("Br",(1,)),

    ("I",(1,))



    ])

colors_list=['white',

             'white',

             'white',

             'white',

            'ghostwhite',

            'grey',

            'blue',

            'red',

            'cyan',

            'darkorange',

            'yellow',

            'lime',

            'darkred',

            'purple',

            'black'

                     ]

ATOMdecode = OrderedDict(ATOM_type_list)



BANDdecode = OrderedDict(


    [(0.0,"0"),

     (0.3,"1"),

     (0.5,"2"),

     (0.6,"3"),

     (0.35,"4"),

     

    (0.1,"11"),

    (0.15,"12"),

    (0.2,"13"),

    (0.25,"14")]





    )

#print(ATOMdecode)

#input(ATOMdecode[0.2])



###############


def getLineElement(file_name):     #找到文件每行元素输出为二维列表



    mode=file_name[-3:]



    with open(file_name,"r") as txt:

        

        pattern = re.compile(r'\S+')

        

        lines_info_list=[]

        

        for line in txt.readlines():



            lines_info_list.append ( pattern.findall(line))

    

    return lines_info_list ,mode





def get_mol_info(info):     #将输入的二维文本列表转化为原子类型，原子坐标，键连关系

    Lines,mode = info[0],info[1]

    

    if mode == "sdf" or mode == "mol":

        atom_info = []

        atom_bond = []

        for cont, line in enumerate( Lines):

            #print(cont,line)

        

            if len(line)==16:

            

                atom_info .append(line)

            if len(line)==4 :  #4

                if "CHG" in line:

                    break

                        

                atom_bond .append(line[:3])

            

            if len(line)==3 : #3  #超过99个分子的mol文本粘连问题

                if "CHG" in line:

                    break

                

                #print(line)

                #input()

                if len(line[0])==4:           

                    

                    #print(bond_line_str[:3],"-----", bond_line_str[3:])

                    atom_bond .append( [   line[0][:-3], line[0][-3:]  ,line[1] ])

                if len(line[0])==5:           

                    

                    #print(bond_line_str[:3],"-----", bond_line_str[3:])

                    atom_bond .append( [   line[0][:-3], line[0][-3:] ,line[1]  ])

                if len(line[0])==6:           

                    

                    #print(bond_line_str[:3],"-----", bond_line_str[3:])

                    atom_bond .append( [   line[0][:-3], line[0][-3:] ,line[1]  ])



        atom_bond_out = []

        for c,i in enumerate(atom_bond):

            #print("atom_bond",c+1,i)

            atom_bond_out.append([str(c+1), i[0] ,i[1],i[2] ])





              

        atom_type = {}

        atom_pos = []

        for count,i in enumerate (atom_info):      

            atom_pos.append(  (float(i[0]),float(i[1]),float(i[2]))   )

            for ATOM_type in  ATOM_type_list:

                if   ATOM_type[1] == i[3]:

                    atom_type[str(count+1)] = (ATOM_type[1],ATOM_type[0])

                if i[3] not in list(zip(*ATOM_type_list))[1]:

                    print(i[3])

                    atom_type[str(count+1)] = ("Other",1.1)

                

                



        return atom_type , atom_pos, atom_bond_out






###############





class MolNets():

   



    def __init__(self,file_name=""):

        

        self.file_name = file_name

        self.atom_type,self.pos ,self.band  = get_mol_info(getLineElement(self.file_name))

        self.atom_num = len(self.pos)





    def get_C_list(self,file_name):      #为了防止初始结构种类太多，根节点从C出发

        

        C_list=[]

        atom_type,pos ,band  = self.atom_type,self.pos ,self.band 

        for i,j in atom_type.items():

            if j[0]=="C":

                C_list.append(int(i))

        C_list.sort()      

        return C_list

    

    

    def have_other_type(self,file_name):

        atom_type,pos ,band  = self.atom_type,self.pos ,self.band 

        for i,j in atom_type.items():

            #print(i,j)

            if j[0]=="Other":

                return True

    

        return False

        

    

    

    def branch(self,band,core):    #找到某个原子周围相键连的所有原子，core中心原子

        bud = []

        for line in band:

        

            #print(line[1],line[2])

            

            if core in (int(line[1]),int(line[2])):

                #print(line[1:3])

                bud.append(int(line[1]))

                bud.append(int(line[2]))

    

        while bud.count(core)>0:

            bud.remove(core)

        #print(bud)

        return bud

    

    

    


    #root=1

    

    ##########

    


    

    def mol2net(self,file_name,root):

        """

        Algorithm 2 placeholder

        """





        graph = {}

        atom_type,pos,band = self.atom_type,self.pos ,self.band





    

        for i in atom_type:

    

            graph[str(i)] =   [str(x) for x in self.branch(band ,int(i))]

            #print(i,[str(x) for x in self.branch(band ,int(i))])





        

        L = []

        placeholder_atom = {}

        placeholder_dir = {}

        placeholder_atom_num = self.atom_num

        width_list  = []

        

        l_1 = [str(root)]

        L.append(l_1)

        width_list.append(len(l_1))

        

        l_2 = graph[str(root)]

        #print(l_2,"l_2",graph)

        L.append(l_2)

        width_list.append(len(l_2))





        

        n = 2

    

        while True:

            

            s = set()

            for v in L[n-1]:

                if "-" not in v:

                    s.update(graph[v])

                

            l_n = [ v for v in s if v not in L[n-1] and v not in L[n-2]]

    

            ###########add placeholder

            l_inner_connect_vs = [ v for v in s if v in L[n-1] ]

            inner_connect_pairs = [  set([v,v_c])  for v in  l_inner_connect_vs for v_c in graph[v] if v_c in l_inner_connect_vs  ]

            connect_pairs_not_repeat = []

            for pair in  inner_connect_pairs:

                if pair not in connect_pairs_not_repeat:

                    pair = sorted(list(pair),key  =lambda x:int(x))

                    connect_pairs_not_repeat.append(pair)

                    

                    placeholder_atom_num+=1

                    

                    

                    #print(pair,pair)

                    placeholder_atom["-".join(pair)]=str(placeholder_atom_num)

                    

                    placeholder_dir["-".join(pair)]= [ [pair[0],str(placeholder_atom_num)],[pair[1],placeholder_atom_num]  ]

            

            placeholder_vs = [ "-".join( sorted(list(pair),key  =lambda x:int(x)) ) for pair in connect_pairs_not_repeat]

            l_n.extend(placeholder_vs)

            

            

            

            ###########add placeholder

            

            if l_n == [] :

                break



            #################

            ########### sort layer vertex by order,avoid edge intersection


            #################




    

            L.append(l_n)

            width_list.append(len(l_n))

            

            n +=1

            

        Widest_width = max(width_list)



        #print("placeholder_atom",placeholder_atom)

        #print("placeholder_dir",placeholder_dir)

        #print("L",L)

       


        nets = []

        for index in range(len(L)):

            net = []

            if index>0:

                #print(index,L[index-1],L[index])

                for index_v, v in enumerate(L[index-1]):

                    if "-" not in v:

                        for index_v_c,v_c in enumerate(L[index]):

                            if "-" not in v_c  :

                                if v_c in graph[v]:

                                    net.append([v,v_c])

                            elif  v in v_c.split("-"):

                                    net.append([v,placeholder_atom[v_c]])

                                    #print("----",[v,v_c],placeholder_atom[v_c])

    

    

                nets.append(net)

        #print("nets",nets)
        #input()
        return nets ,  Widest_width-1 , len(L)-1 ,placeholder_dir

        #return  Layers_net , max(x_size) ,Layers_long,placeholder_dir

    



    ##############################################################################################33

    

    #die = [root]
     #将树表示为矩阵
    #-------------------------------------------------

    




    #Layers_net

    

    def net2matrix(self,file_name,root='1',size=[20,15,15,5]):
        root = str(root)
        

    

        atom_type,pos ,band = self.atom_type,self.pos ,self.band   #获得分子mol文件属性

        nodes,_,_,placeholder_dir = self.mol2net(file_name,root)                                 #输出为连接关系的列表，二维列表：第一维是层，第二维是分子连接关系

    

        

    ############################为分子增加占位符信息

    

        band_type_dir = {}

        for i in band:                     #建立成键与成键类型字典  如"1-2":1.  1 2 3 4 分别为单 双 三 占位键

            if int(i[1])<int(i[2]):

                band_type_dir["-".join([i[1],i[2]])]=i[3]

            else:

                band_type_dir["-".join([i[2],i[1]])]=i[3]

    

        

        #print("band",band)

        #print(band_type_dir)

        

        all_bond_num = len(band)

        for key,item in placeholder_dir.items():

    

            

            if int(item[0][0])< int(item[1][0]):

                    band_name = str(item[0][0]) +"-"+ str(item[1][0])

    

            else:

                    band_name = str(item[1][0]) +"-"+ str(item[0][0])

    

                                    

            #print(band_name,item)

            if band_type_dir [ band_name  ]=="1":

                #print("item ",item )

            

                atom_type[str(item[0][1])] =  (ATOM_type_list[0][1],ATOM_type_list[0][0])   #('h1', 0.12)

                pos.append((0.0,0.0,0.0))#################应该是相加除二

            

                all_bond_num+=1

                

                band.append([str(all_bond_num),str(item[0][0]),str(item[0][1]),str(11)])

    

                all_bond_num+=1

                

                band.append([str(all_bond_num),str(item[1][0]),str(item[1][1]),str(11)])

    

    

            if band_type_dir [ band_name ]=="2":

                #print("item ",item )

            

                atom_type[str(item[0][1])] =(ATOM_type_list[1][1],ATOM_type_list[1][0])  #('h2', 0.13) #

                pos.append((0.0,0.0,0.0))

            

                all_bond_num+=1

                

                band.append([str(all_bond_num),str(item[0][0]),str(item[0][1]),str(12)])

    

                all_bond_num+=1

                

                band.append([str(all_bond_num),str(item[1][0]),str(item[1][1]),str(12)])

    

            if band_type_dir [ band_name  ]=="3":

                #print("item ",item )

            

                atom_type[str(item[0][1])] =(ATOM_type_list[2][1],ATOM_type_list[2][0])  #('h3', 0.14)

                pos.append((0.0,0.0,0.0))

            

                all_bond_num+=1

                

                band.append([str(all_bond_num),str(item[0][0]),str(item[0][1]),str(13)])

    

                all_bond_num+=1

                

                band.append([str(all_bond_num),str(item[1][0]),str(item[1][1]),str(13)])

    

            if band_type_dir [ band_name  ]=="4":

                #print("item ",item )

            

                atom_type[str(item[0][1])] =(ATOM_type_list[3][1],ATOM_type_list[3][0])  #('h4', 0.15)

                pos.append((0.0,0.0,0.0))

            

                all_bond_num+=1

                

                band.append([str(all_bond_num),str(item[0][0]),str(item[0][1]),str(14)])

    

                all_bond_num+=1

                

                band.append([str(all_bond_num),str(item[1][0]),str(item[1][1]),str(14)])

                

    

            

        band_type_dir = {}

        for i in band:                     #建立成键与成键类型字典  如"1-2":1.  1 2 3 4 分别为单 双 三 占位键

            if int(i[1])<int(i[2]):

                band_type_dir["-".join([i[1],i[2]])]=i[3]

            else:

                band_type_dir["-".join([i[2],i[1]])]=i[3]

                

            

            

            

            

            

    ###############################################

    

        

        #for key ,item in atom_type.items():

            

            #print("key ,item",key ,item)

    

                                            

    

        b=[[root]]                                      #创建每层含有的原子列表

        #b_order = [[[root,atom_type[str(root)][1],0]]]   

        

    

    

        #已知一个顶点，建立一个能获取连接该顶点的前一个顶点的序数的字典

        last_nodes_dicts = []
        next_nodes_dicts = []

        for i in nodes:

            

            last_nodes_dicts.append(dict([  [elem[1],[]]  for elem in i]))
            next_nodes_dicts.append(dict([  [elem[0],[]]  for elem in i]))


        


        for cont, i in enumerate (nodes):               #遍历连接关系列表

            #print(i)
            for elem in i :
                
                last_nodes_dicts[cont][elem[1]].append(elem[0]) #找到与上层连接的定点
               
                next_nodes_dicts[cont][elem[0]].append(elem[1]) #找到与下层连接的定点

            #print("next_nodes_dicts",next_nodes_dicts)
            

        #input()
            

      

    #########################################    对于不同元素的排序

    
        #print("nodes",nodes)
        for cont, i in enumerate (nodes):               #遍历连接关系列表



            y=[]

            

            for j in i:

    

                #print(j[1])

    

                y.append((j[1],atom_type[str(j[1])][1] ))

            

            y=list(set(y))                              #去除重复

            

            y.sort(key=lambda x:-x[1])   #以元素周期表排序

                       

            #y.sort(reverse=True)

            

            #if cont >-1: #从第二层开始对同种元素交换数序

            ys_order = []   #单层包含的节点

    

            

    ###############################################         对于同种元素的排序   

            for y_elem in y:

                #print("y_elem",y,y_elem)
                last_link_nodes = list(set(last_nodes_dicts[cont][y_elem[0]]))#找到与该顶点连接的上一层顶点
                #print(next_nodes_dicts)
                
                #print(cont,len(nodes),len(next_nodes_dicts))
                if cont+1<len(next_nodes_dicts)-1:
                    next_link_nodes = list(set(    next_nodes_dicts[cont+1].get(y_elem[0],"0")    ))#找到与该顶点连接的下一层顶点


                #print(b[cont],y_elem,[last_link_node])

                

                order_1 =  [ b[cont].index(last_link_node) for last_link_node in last_link_nodes ]  #按照与上层连接的顶点次序对该顶点设定等级，连接多个顶点就取平均

                order_2 =  [ 3**(atom_type[next_link_node][1]*10) for next_link_node in next_link_nodes if next_link_node!="0"]  #按照与上层连接的顶点次序对该顶点设定等级，连接多个顶点就取平均
                #print(order_2)

                    #print(y_elem[0],'-',order ,end=" ")


                ys_order.append( [y_elem[0],atom_type[str(y_elem[0])][1],sum(order_1)/len(order_1),sum(order_2)])

    

                    

                

                #y_order .sort(key=lambda x:x[1])

                

                

                #b_order.append(list(zip(*y_order))[0])

                #print("y_order",y_order)

                #print("y_order_",sort_atom(y_order))

                #print("b",b[cont])
            #print("ys_order",ys_order)
            ys_order =  sorted(ys_order,key = lambda x :-x[3]) #与下层连接的原子
            ys_order =  sorted(ys_order,key = lambda x :x[2])  #与上层连接的原子
            ys_order =  sorted(ys_order,key = lambda x :-x[1])  #此原子的种类
    
            #ys_order_0= list(zip(*   self.sort_atom(ys_order)   ))[0]
            #print("ys_order",ys_order)
            ys_order = [ y_order[0] for y_order  in ys_order ]

                #b_order.append(y_order_0)

    ##########################################################################

                #print("y_elem:")

            #b_order [0]=[root]

            #print("b_order",b_order)          #

            #y=list(zip(*y))[0]  #取每一对的第一个元素

            #print("y",y)

            b.append(ys_order)                               #每层有顺序的原子列表

            #print(b)

            

    

        #b为layers

            #print("--b--",b)

    

    

        #print()

        chain = np.zeros(size)                          #新建一个矩阵，函数会返回之。第一维是层数，第二维第三维是原子可以连接的最大数目X_size。第四维是前后分子信息及分子的连接向量 

        #b=b_order

        

        #cont统计层数

        for cont, i in enumerate (nodes):

    

            for j in i:

                for m in range(len(b[cont])):                            #定义了二维表，在这个有顺序的表里找元素打勾

                    for n in range(len(b[cont+1])):

    

                        if j==[b[cont][m],b[cont+1][n]]:

    

                            

    

                            if len(size)==3:

                                chain[cont,m,n] = 1.0

                            if len(size)==4 and size[3]==2:########彼此连接的元素类型

                                #print("j",j)

                                chain[cont,m,n] = [ atom_type[ str( j[0] ) ][1]  ,atom_type[ str(j[1]) ][1]     ]  #BANDencode

    

    

    

                            if len(size)==4 and size[3]==3:#######彼此连接的元素类型与成键类型

                                

                                if j[0]< j[1]:

                                    band_name = "-".join( [ str(j[0]), str(j[1]) ])

    

                                else:

                                    band_name = "-".join( [ str(j[1]), str(j[0]) ] )

                

    

                                 

                                ##print("band_type_encode ",band_type_encode )

    

    

                                band_name = "-".join(  sorted( [ x for x in  band_name.split("-")] ,key = lambda x :int(x) )  )

                                chain[cont,m,n] = [ atom_type[ str( j[0] ) ][1]  ,atom_type[ str(j[1]) ][1]  ,BANDencode[band_type_dir[band_name]]   ]

                                #print(atom_type[ str( j[0] ) ][1]  ,atom_type[ str(j[1]) ][1]  ,band_type_encode)

    

    

    

    

                                

                            if len(size)==4 and size[3]==5 :

                                p =   np.array(pos[j[1]-1])  -  np.array(pos[ j[0]-1] )

                                chain[cont,m,n] = [ atom_type[ str( j[0] ) ][1]  ,atom_type[ str(j[1]) ][1],p[0],p[1],p[2]     ]                  

    

    
        #print("b",b)
        return chain ,len(nodes),b  #矩阵 长度  层顶点

    
    def plot_nets(self,file_name,root,size,show_num=False,show_info=False):

    

    

        atom_type,pos ,band = self.atom_type,self.pos ,self.band   #获得分子mol文件属性

        

        nodes,_,_,placeholder_dir = self.mol2net(file_name,root)                                 #输出为连接关系的列表，二维列表：第一维是层，第二维是分子连接关系

    

    

        

    ############################为分子增加占位符信息

        band_type_dir = {}

        for i in band:                     #建立成键与成键类型字典  如"1-2":1.  1 2 3 4 分别为单 双 三 占位键

            if int(i[1])<int(i[2]):

                band_type_dir["-".join([i[1],i[2]])]=i[3]

            else:

                band_type_dir["-".join([i[2],i[1]])]=i[3]

    

        

        #print("band",band)

        #print(band_type_dir)

        

        all_bond_num = len(band)

        for key,item in placeholder_dir.items():

    

            

            if int(item[0][0])< int(item[1][0]):

                    band_name = str(item[0][0]) +"-"+ str(item[1][0])

    

            else:

                    band_name = str(item[1][0]) +"-"+ str(item[0][0])

    

                                    

    
            print(band_name)
            if band_type_dir [ band_name  ]=="1":

                #print("item ",item )

            

                atom_type[str(item[0][1])] =  (ATOM_type_list[0][1],ATOM_type_list[0][0])   #('h1', 0.12)

                pos.append((0.0,0.0,0.0))#################应该是相加除二

            

                all_bond_num+=1

                

                band.append([str(all_bond_num),str(item[0][0]),str(item[0][1]),str(11)])

    

                all_bond_num+=1

                

                band.append([str(all_bond_num),str(item[1][0]),str(item[1][1]),str(11)])

    

    

            if band_type_dir [ band_name ]=="2":

                #print("item ",item )

            

                atom_type[str(item[0][1])] =(ATOM_type_list[1][1],ATOM_type_list[1][0])  #('h2', 0.13) #

                pos.append((0.0,0.0,0.0))

            

                all_bond_num+=1

                

                band.append([str(all_bond_num),str(item[0][0]),str(item[0][1]),str(12)])

    

                all_bond_num+=1

                

                band.append([str(all_bond_num),str(item[1][0]),str(item[1][1]),str(12)])

    

            if band_type_dir [ band_name  ]=="3":

                #print("item ",item )

            

                atom_type[str(item[0][1])] =(ATOM_type_list[2][1],ATOM_type_list[2][0])  #('h3', 0.14)

                pos.append((0.0,0.0,0.0))

            

                all_bond_num+=1

                

                band.append([str(all_bond_num),str(item[0][0]),str(item[0][1]),str(13)])

    

                all_bond_num+=1

                

                band.append([str(all_bond_num),str(item[1][0]),str(item[1][1]),str(13)])

    

            if band_type_dir [ band_name  ]=="4":

                #print("item ",item )

            

                atom_type[str(item[0][1])] =(ATOM_type_list[3][1],ATOM_type_list[3][0])  #('h4', 0.15)

                pos.append((0.0,0.0,0.0))

            

                all_bond_num+=1

                

                band.append([str(all_bond_num),str(item[0][0]),str(item[0][1]),str(14)])

    

                all_bond_num+=1

                

                band.append([str(all_bond_num),str(item[1][0]),str(item[1][1]),str(14)])

    

    

                

    

            

        band_type_dir = {}

        for i in band:                     #建立成键与成键类型字典  如"1-2":1.  1 2 3 4 分别为单 双 三 占位键

            if int(i[1])<int(i[2]):

                band_type_dir["-".join([i[1],i[2]])]=i[3]

            else:

                band_type_dir["-".join([i[2],i[1]])]=i[3]

                

            

            

            

    ###############################################

            

            

            

            

    ###############################################

    

                

        

        b= self.net2matrix(file_name,root,size)[2]#########有顺序的每层顶点

    

    

            

        b_layers=[]                #带信息的每层排列好的顶点

        for count,i in enumerate(b):

            b_layer=[]

            #print(">>Layer ",end="")

            for j in i:

                #print(j,atom_type[str(j)][0],end="   ")#给每一层顶点追加原序数

                b_layer.append((j,atom_type[str(j)][1],atom_type[str(j)][0]))

            #b_layer.sort(key=lambda x:-x[1])

                

            b_layers.append(b_layer)

            #print()

        #print("b_layers",b_layers)

        #input()

        

    

    

    

            

        x_plot=[]

            

        y_plot=[]

            

            #type_order_0=0

            #type_order_1=0

    

    

     
        #print("nodes",nodes)
        for count, i in enumerate (nodes):

                #print("->",i)

    

                #print("x_0",x_0)

             

            y_plot_list_0 =list(  zip(b[count],list(range(len(b[count]))) )  )

            y_plot_list_1 =list(  zip(b[count+1],list(range(len(b[count+1]))) )  )

    

                #y_plot_point_0=0

                #y_plot_point_1=0

            for j in i:

                #print("j",j,end="")

                x_plot.append([count,count+1])

                    

            for y_0 in y_plot_list_0 :

                for y_1 in y_plot_list_1:

                    #print("1",[y_0[0],y_1[0]])

                    for node in nodes[ count]:

                        #print("2",node,[y_0[0],y_1[0]])

                        if [y_0[0],y_1[0]] == node:

                                #print("yes")
                            #print("---------------",[str(b[count][y_0[1]]),str(b[count+1][y_1[1]])]  )

                            

                            y_plot.append([ [y_0[1],y_1[1]] ,[ atom_type[str(y_0[0])][0],atom_type[str(y_1[0])][0] ]   ,[str(b[count][y_0[1]]),str(b[count+1][y_1[1]])] ] )

                            

    

        #print(   x_plot])

    

        #print("nodes",nodes)

        #print("b",b)

        #print("x_plot",x_plot)

        #print("y_plot",y_plot)

        #input()

        #for i in band:

            #print(i)

        #print(x_plot)
        #print("----")
        #print(y_plot)
        
        #print("----")

        fig, ax = plt.subplots()      

        for i in range(len(x_plot)):

            #print("y_plot[i][0]", y_plot[i][2])

    

            if int(y_plot[i][2][0])<int(y_plot[i][2][1]):

                band_name = "-".join(  [y_plot[i][2][0] , y_plot[i][2][1]]  )

    

            else:

                band_name = "-".join(  [y_plot[i][2][1] , y_plot[i][2][0]]  )

                

    

            if band_type_dir[band_name] == "1":

    

                plt.plot(   x_plot[i], y_plot[i][0], color=(i/(len(x_plot)),1-i/(len(x_plot)),1-i/(len(x_plot))),linestyle='-',lw=2,zorder=1)  #画线

    

            if band_type_dir[band_name] == "2":

    

                plt.plot(   x_plot[i], y_plot[i][0], color=(i/(len(x_plot)),1-i/(len(x_plot)),1-i/(len(x_plot))),linestyle='-',lw=4,zorder=1)  #画线

    

            if band_type_dir[band_name] == "3":

    

                plt.plot(   x_plot[i], y_plot[i][0], color=(i/(len(x_plot)),1-i/(len(x_plot)),1-i/(len(x_plot))),linestyle='-',lw=8,zorder=1)  #画线

    

            if band_type_dir[band_name] == "4":

    

                plt.plot(   x_plot[i], y_plot[i][0], color=(i/(len(x_plot)),1-i/(len(x_plot)),1-i/(len(x_plot))),linestyle='-.',lw=4,zorder=1)  #画线

                

            if band_type_dir[band_name] == "11":

    

                plt.plot(   x_plot[i], y_plot[i][0], color=(i/(len(x_plot)),1-i/(len(x_plot)),1-i/(len(x_plot))),linestyle=':',lw=0.5,zorder=1)  #画线

    

            if band_type_dir[band_name] == "12":

    

                plt.plot(   x_plot[i], y_plot[i][0], color=(i/(len(x_plot)),1-i/(len(x_plot)),1-i/(len(x_plot))),linestyle=':',lw=0.5,zorder=1)  #画线

    

            if band_type_dir[band_name] == "13":

    

                plt.plot(   x_plot[i], y_plot[i][0], color=(i/(len(x_plot)),1-i/(len(x_plot)),1-i/(len(x_plot))),linestyle=':',lw=0.5,zorder=1)  #画线

    

            if band_type_dir[band_name] == "14":

    

                plt.plot(   x_plot[i], y_plot[i][0], color=(i/(len(x_plot)),1-i/(len(x_plot)),1-i/(len(x_plot))),linestyle=':',lw=0.5,zorder=1)  #画线

                

        x_plot_=[]

        for i in x_plot:

            x_plot_.extend(i)

              

        y_plot_=[]

        y_plot_text=[]

        for i in y_plot:

            y_plot_.extend(i[0])

            if show_num==False:

                y_plot_text.extend(i[1])

            #print([i[1][0]+i[2][0],i[1][1]+i[2][1]])

            if show_num==True:

                y_plot_text.extend([i[1][0]+i[2][0],i[1][1]+i[2][1]])

                

        

        for i in range(len(x_plot_)):

            #print(y_plot_text[i])

            plt.text(   x_plot_[i], y_plot_[i], y_plot_text[i],color='black',fontsize=15,ha='center', va='center',fontdict={'family' : 'Times New Roman'}  )#va='bottom',

    

    

            for count ,  ATOM_type in  enumerate(ATOM_type_list):

                if  y_plot_text[i]  == ATOM_type[1]:

                    

                    plt.scatter(x_plot_[i], y_plot_[i], color=colors_list[count],marker='o', edgecolors='black', s=200,zorder=2) #画点

    

            

    

        ax.set_axis_off()

        plt.show()

                        

    #plot_plot_plot_plot_plot_plot_plot_plot_plot_plot_plot_plot_plot_

    

    def plot_matrix(self,file_name,root,size  ):

    

        

        

        matrix,layer_num,b= self.net2matrix(file_name,root,size)  #矩阵，长度，层顶点

    

        print("size",size)

        Nr = size[3]

        Nc = layer_num

        

    

            

        print("Nr,Nc",Nr,Nc)

        cmap ='viridis'# "cool"

    

        fig, axs = plt.subplots(Nr, Nc)

        fig.suptitle('Multiple images')

    

        images = []

        for i in range(Nr):

            for j in range(Nc):

                # Generate data with a range that varies from one plot to the next.

                #print(i,j)

                data = matrix[j,:,:,i].T

    

                if Nc ==1:

                

                    images.append(axs[i].imshow(data, cmap=cmap))

                    axs[i].label_outer()

                else:

                    images.append(axs[i, j].imshow(data, cmap=cmap))

                    axs[i, j].label_outer()

        # Find the min and max of all colors for use in setting the color scale.

        vmin = min(image.get_array().min() for image in images)

        vmax = max(image.get_array().max() for image in images)

        norm = colors.Normalize(vmin=vmin, vmax=vmax)

        for im in images:

            im.set_norm(norm)

    

        fig.colorbar(images[0], ax=axs, orientation='horizontal', fraction=.1)

    

    

    

    # Make images respond to changes in the norm of other images (e.g. via the

    # "edit axis, curves and images parameters" GUI on Qt), but be careful not to

    # recurse infinitely!

        def update(changed_image):

            for im in images:

                if (changed_image.get_cmap() != im.get_cmap()

                        or changed_image.get_clim() != im.get_clim()):

                    im.set_cmap(changed_image.get_cmap())

                    im.set_clim(changed_image.get_clim())

    

    

        for im in images:

            im.callbacksSM.connect('changed', update)

        

        plt.show()





    

    def renove_repeat(self,List):

        list_new=[]

        for i in List:

            if i not in list_new:

                list_new.append(i)

        return list_new

    

    def min_x_size_root(self,file_name):

        

        roots = self.get_C_list(file_name)

        root_size_long=[]

        for i in roots:

            nodes ,size, long,_ = self.mol2net(file_name,i)

            

    #        print(i, size, long)  #min size ,long

            root_size_long.append([i, size, long])

        

      

        nice_long_root =sorted(root_size_long ,key=lambda x :-x[2])     #取最长的序列对应根节点

        nice_size_root =sorted(nice_long_root ,key=lambda x :x[1])    #取最窄的“序列最大宽度”对应根节点

        #print(file_name,"nice_size_root",nice_size_root)

        #print(nice_size_root[0])

        nice_root = nice_size_root[0][0]

        return nice_root

        

    def min_x_size_roots(self,file_name):

        

        roots = self.get_C_list(file_name)

        root_size_long=[]

        for i in roots:

            nodes ,size, long,_ = self.mol2net(file_name,i)

            

    #        print(i, size, long)  #min size ,long

            root_size_long.append([i, size, long])

        nice_long_root =sorted(root_size_long ,key=lambda x :-x[2])     #取最长的序列对应根节点

        nice_size_root =sorted(nice_long_root ,key=lambda x :x[1])    #取最窄的“序列最大宽度”对应根节点

        #print(nice_size_root)

        #print(nice_size_root[0])

        if len(nice_size_root)>2:

            return nice_size_root[0][0],nice_size_root[1][0],nice_size_root[2][0]

        

        if len(nice_size_root)==2:

            return nice_size_root[0][0],nice_size_root[1][0] 

    

        if len(nice_size_root)<2:

            return [ nice_size_root[0][0] ]

    

    print("-----------------------------------------")##############################################################################

    

   



    



class MolTools():







    def branch(self,band,core):    #找到某个原子周围相键连的所有原子，core中心原子

        bud = []

        for line in band:

        

            #print(line[1],line[2])

            

            if core in (int(line[1]),int(line[2])):

                #print(line[1:3])

                bud.append(int(line[1]))

                bud.append(int(line[2]))

    

        while bud.count(core)>0:

            bud.remove(core)

        #print(bud)

        return bud

    ###############

    #file_name = "088602.mol"

    def get_valence(self,lines):

    

        pattern = re.compile(r'\S+')

            

        lines_info_list=[]

            

        for line in lines:

    

            lines_info_list.append ( pattern.findall(line))

    

        

        info = lines_info_list

        #print("info",info)

        atom_type,pos ,band  = get_mol_info([info,"mol"])

        #print("band",band)

    

        atoms_num = len(pos)

    

    

        band_type_dir = {}

        for i in band:                     #建立成键与成键类型字典  如"1-2":1.  1 2 3 4 分别为单 双 三 占位键

            if int(i[1])<int(i[2]):

                band_type_dir["-".join([i[1],i[2]])]=i[3]

            else:

                band_type_dir["-".join([i[2],i[1]])]=i[3]

    

        #for k,v in band_type_dir.items():

            #print(k,v)

    

    

    

        atom_valence_list = []

        

        

        for i in range(1,atoms_num+1):

            #print(i,atom_type[str(i)][0],branch(band,i),end = " ")

            #print(i,atom_type[str(i)][0],end = " ")

            atom_valence = 0

            for j in self.branch(band,i):

                if int(i)<int(j):

                    band_name = str(i)+"-"+str(j)

                else:

                    band_name = str(j)+"-"+str(i)

                

                atom_valence+=int(band_type_dir[band_name])

                #print('----------------',band_type_dir[band_name],end = " ")

            #print(atom_valence)

            atom_valence_list.append(  [str(i) , [atom_type[str(i)][0],atom_valence]]  )

                                            

        #for i in band :

            

            #print(i)

        #print("atom_valence_list",atom_valence_list)

        return atom_valence_list

    

    

                

    

    

    

    def have_right_atom_valence(self,lines):

        atom_valence_list = self.get_valence(lines)    

        for i in atom_valence_list:

        

            #print('-->',i,i[1][0],i[1][1])

            if i[1][1] not in ATOM_valence[i[1][0]]:

                return False

        return True

    #print(have_right_atom_valence(file_name))      

    #input()

 

    ###############

    ##########    

    

    #def build_train_test_data(save_name):

        #with open("gdb9.sdf_gdb_129135removed_.csv","r") as csvfile:

            #reader=csv.reader(csvfile)

            #floder="SDF/"

            #sdf_E = list( zip(  sorted(list(os.listdir(floder)) ,key=lambda x:int(x[:-4]))[0:133883]  ,   list(reader)[1:133884]) )

    

        #random.seed(0)

        #random.shuffle(sdf_E)

    

    

    

        #sdf_data_x=[]

        #sdf_data_y=[]

        #sdf_data_L=[]

        #count=0

    

        #dx = int(len(sdf_E)/5)

    

        #cat_num = -1

        

        #for i in range(5):

    

            #cat_num+=1

            #start = i*dx+cat_num

            #end = (i+1)*dx+cat_num

            #if end == 133880:

                #end = 133883

            

            #for sdf , E in sdf_E[ start :  end  ] :

            #print(floder+sdf,E[1:])

            #print(count)

            #input()

                #count+=1

             

                #if len(self.get_C_list(floder+sdf )) > 0:

                    #root=self.min_x_size_root(floder+sdf)

                    #L = self.mol2net(floder+sdf,root)[0]

                    #n = self.net2matrix(floder + sdf ,root,[10,12,12,3])[0]

    

    

                    #sdf_data_L.append( sdf )

    

                    #sdf_data_x.append( n.reshape(12*10,12,3) )

    

                    #sdf_data_y.append( np.array([float(e) for e in E[1:] ]) )

    

                    #print(np.array([float(e) for e in E[1:] ]))

                    

                    #print(count,sdf)

    

    

            #print(len(sdf_data_L))

    

            #sdf_data_train = [  sdf_data_L ,  sdf_data_x , sdf_data_y]

        

            #print(len(sdf_data_train))

    

            #with open(save_name+str(i),"wb") as f :

            

                #pickle.dump(sdf_data_train,f)

                

    

    

    #save_name = "sdf_data_train" 

    #start_cut = 0

    #end_cut = 80000

    #build_train_test_data(start_cut,end_cut ,save_name)

    

    #-----------------------------------------------remove   max(x_size)>10



                

    def mol_to_sequence(self,file_name,nets_dir,root=0):

        sequence=[]

        Net  = MolNets(file_name)

        if root==0:

            

            root=Net.min_x_size_root(file_name)

        L,w_num = Net.mol2net(file_name,root)[:2]

    

        if w_num<=15:

            mol =Net.net2matrix(file_name ,root,[len(L)+1,15,15,3])[0]

            #mol = np.load(to_files_name+file)



            for i in mol :


                s_r = str(i).replace(" [0.  0.  0. ]\n","").replace("\n","").replace("             ","")

                #print(nets_dir)
                if  s_r  in nets_dir.keys():

                    sequence.append(nets_dir[ s_r  ][1])

                    #if s[str(i)][1] >8282:

                        #print(s[str(i)][1],end="#")

                    #else:

                        #print(s[str(i)][1],end=" ")

    

                else:

                    #print("99999",end=" ")

                    sequence.append(99999)

        mol_str = " ".join([str(x) for x in sequence if x!=0])

        return mol_str 

    

    #files_name = "SDF/"

    #for i in range(1,5):

        

    #    print(i,mol_to_sequence(files_name+str(i)+".sdf","nets_info_dir_NPname_3000.plk"))

    #input("######################################")

    #######################################

    

    

    ######################################################################################

    def mol_to_nets_dir(self,files_list,dir_save_name,key_is_order=False):

        #files_list = os.listdir(folder)

        #files_list.sort(key=lambda x:int(x[:-4]))

        #file_dir = file_dir[:]

        

        

    

    

    #b= np.load( folder+"1..matrix.npy" )

        nets_str = []
        nets_list = []

        nets_dir = {}

        count = 0

        write_txt_file_list = []
        roots_list = []

        for file_name in files_list:

            #nets = np.load(folder + i )

            #print(self.get_C_list(folder+i))

            Net  = MolNets(file_name)

            if len(Net.get_C_list(file_name)) >0:


                root=Net.min_x_size_root(file_name)

                #root = root_info[cont][0]

                #count+=1

    

                
                
                print(file_name,root)

                L,w_num = Net.mol2net(file_name,root)[:2]

    

                if w_num<=15:

                    write_txt_file_list.append(file_name)
                    roots_list.append(root)

            

                    nets =Net.net2matrix(file_name ,root,[len(L)+1,15,15,3])[0]

                    nets_str .extend(list(map( lambda x:str(x).replace(" [0.  0.  0. ]\n","").replace("\n","").replace("             ",""),nets )))
                    
                    nets_list.extend( list(nets) )

                    

                    

                    



        nets_dir = dict( list(zip( nets_str,nets_list  )) )
        

        
        #print(len(nets_dir.keys()))
        #input(nets_dir[nets_str[1]])
        
        c = Counter(nets_str)

        count =0

        net_a=[]

        net_b=[]

        for i ,j in c.items():

            net_a.append(i)

            net_b.append(j)

    

        net_count=list(zip(net_a,net_b))

    

    

    

        #def take_second(elem):

            #return -elem[1]

    

        net_count .sort(key=lambda x :-x[1]) #元素 ，个数

    

        c_1=0

        c_2=0

        c_3=0

        c_4=0

        c_5=0

        nets_info_dir ={}

        for count,i in enumerate( net_count):

            print(i[1])

            if key_is_order==False:

                nets_info_dir[i[0] ]=[ nets_dir[i[0]],count ,i[1]]   #{ 矩阵名字：[矩阵，排序，出现次数] } #编码

            if key_is_order==True:

                nets_info_dir[ count]=[ nets_dir[i[0]],i[0] ,i[1]]   #{ 排序：[矩阵，矩阵名字，出现次数] }

                

            if i[1]==1:

                c_1+=1

            if i[1]==2:

                c_2+=1

            if i[1]==3:

                c_3+=1

            if i[1]==4:

                c_4+=1

            if i[1]==5:

                c_5+=1

        print(c_1,c_2,c_3,c_4,c_5)

        print(len(net_count))

        print("-----")

    

    

        with open(dir_save_name,"wb") as f :

            pickle.dump(nets_info_dir,f)
        return roots_list

            

    

    #print("-----------------------------------------数据库转化为文本")

    def write_txt(self,save_txt_name,nets_dir,files_list):

        

        #files_list = os .listdir(folder)

    

    

    

        

        

        #files_list.sort(key=lambda x:int(x[:-4]))

        #file_dir = file_dir[:]

        

        with open(save_txt_name,"a") as TXT :

            cont=0

            for file_name in files_list:

                #root = root_info[cont][0]

                

                Net  = MolNets(file_name)

                if len(Net.get_C_list(file_name)) >0:

                    root=Net.min_x_size_root(file_name)

                    #root = roots_list[cont]

                    print( file_name,root)

                    cont+=1

                    TXT.writelines(  str(x)+" " for x in    self.mol_to_sequence(file_name,nets_dir,root)  )

                    TXT.writelines("    "+file_name+"\n")

    

    #key为矩阵的字典转化成key为数序的字典

    

    def net_dir2order_dir(self,nets_dir_name,order_dir_name):

    

        

        with open(nets_dir_name,"rb") as f:

            

            nets_dir = pickle.load(f)

    

        #print(nets_dir.shape())

        order_dir={}

        

        for i,j in nets_dir.items():

    

            order_dir[j[1]]=j[0]

    

    

    

        with open(order_dir_name,"wb") as f :

            pickle.dump(order_dir,f)

    

    

    ##########################################################################################################

    #矩阵转化为二分图

    def matrix2layer(self,M):  #(13, 12, 12, 3)

        M=M.reshape(1,15,15,3)

        node_cont=0

        NETs = []

        LAYERs = []

        band_type_dir={}

        for layer_cont, m in enumerate(M):

            net = []

            #layer=[]

            #print("m",m==[0. , 0.  ,0. ])

            for index_x ,array in enumerate (m==[0. , 0.  ,0. ]):

                for index_y ,elem in enumerate(array):

                    #print("array ",array )

                    if elem[0]==False:

                        

                        type_array=M[layer_cont,index_x,index_y,:3] 

                        net.append([

                            (index_x ,ATOMdecode[type_array[0]])  ,

                            

                            (index_y ,ATOMdecode[type_array[1]])

    

                            ])

                        

                        band_type_dir[str(index_x)+"-"+str(index_y)]=BANDdecode[type_array[2]]

                        node_cont+=1

            #print(len(layer))            

            NETs.append(net)

     

            try :

                l1= list(set(list(zip(*net))[0]))  #-> IndexError: list index out of range

                l2= list(set(list(zip(*net))[1]))
                
            except IndexError:
                
                return  "IndexError","IndexError","IndexError"

            l1.sort(key = lambda x:x[0])

            l2.sort(key = lambda x:x[0])

            if layer_cont ==0:

                LAYERs .append(  l1     )

                

            LAYERs .append(  l2      )

                

    

    

        return NETs,LAYERs,band_type_dir

    

    #print()

    #input("finish!")

    #################################################################

    #文本转化为矩阵

    

    def txt2matrix(self,mol_str,order_dir):

    



    

            #print(order_dir[int(mol_str)].shape)

            #print(mol_str)

        return order_dir[int(mol_str)]

            

    #nets_dir_name = "nets_info_dir_NPname_3000.plk"

    #order_dir_name = "order_info_dir_NPname_3000.plk"

    

    #print("txt2matrix",txt2matrix("2",order_dir_name))

    #print("self.matrix2layer(M)",self.matrix2layer(txt2matrix("2",order_dir_name)))

    #input()

    

    def str_plot_nets(self,atom_type,b,nodes,band_type_dir,show_num=False):

    

    

    

        x_plot=[]

            

        y_plot=[]

            

            #type_order_0=0

            #type_order_1=0

    

    

     

        for count, i in enumerate (nodes):

                #print("->",i)

    

                #print("x_0",x_0)

             

            y_plot_list_0 =list(  zip(b[count],list(range(len(b[count]))) )  )

            y_plot_list_1 =list(  zip(b[count+1],list(range(len(b[count+1]))) )  )

    

                #y_plot_point_0=0

                #y_plot_point_1=0

            for j in i:

                #print("j",j,end="")

                x_plot.append([count,count+1])

                    

            for y_0 in y_plot_list_0 :

                for y_1 in y_plot_list_1:

                    #print("1",[y_0[0],y_1[0]])

                    for node in nodes[ count]:

                            #print("2",node)

                        if [y_0[0],y_1[0]] == node:

                                #print("yes")

                            

                            y_plot.append([ [y_0[1],y_1[1]] ,[ atom_type[str(y_0[0])][0],atom_type[str(y_1[0])][0] ]   ,[str(b[count][y_0[1]]),str(b[count+1][y_1[1]])] ] )

                            

    

        #print(   x_plot])

    

        #print("nodes",nodes)

        #print("b",b)

        #print("x_plot",x_plot)

        #print("y_plot",y_plot)

        #input()

        #for i in band:

            #print(i)

        

        fig, ax = plt.subplots()      

        for i in range(len(x_plot)):

            #print("y_plot[i][0]", y_plot[i][2])

    

            if int(y_plot[i][2][0])<int(y_plot[i][2][1]):

                band_name = "-".join(  [y_plot[i][2][0] , y_plot[i][2][1]]  )

    

            else:

                band_name = "-".join(  [y_plot[i][2][1] , y_plot[i][2][0]]  )

                

    

            if band_type_dir[band_name] == "1":

    

                plt.plot(   x_plot[i], y_plot[i][0], color=(i/(len(x_plot)),1-i/(len(x_plot)),1-i/(len(x_plot))),linestyle='-',lw=2,zorder=1)  #画线

    

            if band_type_dir[band_name] == "2":

    

                plt.plot(   x_plot[i], y_plot[i][0], color=(i/(len(x_plot)),1-i/(len(x_plot)),1-i/(len(x_plot))),linestyle='-',lw=4,zorder=1)  #画线

    

            if band_type_dir[band_name] == "3":

    

                plt.plot(   x_plot[i], y_plot[i][0], color=(i/(len(x_plot)),1-i/(len(x_plot)),1-i/(len(x_plot))),linestyle='-',lw=8,zorder=1)  #画线

    

            if band_type_dir[band_name] == "4":

    

                plt.plot(   x_plot[i], y_plot[i][0], color=(i/(len(x_plot)),1-i/(len(x_plot)),1-i/(len(x_plot))),linestyle='-.',lw=4,zorder=1)  #画线

                

            if band_type_dir[band_name] == "11":

    

                plt.plot(   x_plot[i], y_plot[i][0], color=(i/(len(x_plot)),1-i/(len(x_plot)),1-i/(len(x_plot))),linestyle=':',lw=0.5,zorder=1)  #画线

    

            if band_type_dir[band_name] == "12":

    

                plt.plot(   x_plot[i], y_plot[i][0], color=(i/(len(x_plot)),1-i/(len(x_plot)),1-i/(len(x_plot))),linestyle=':',lw=0.5,zorder=1)  #画线

    

            if band_type_dir[band_name] == "13":

    

                plt.plot(   x_plot[i], y_plot[i][0], color=(i/(len(x_plot)),1-i/(len(x_plot)),1-i/(len(x_plot))),linestyle=':',lw=0.5,zorder=1)  #画线

    

            if band_type_dir[band_name] == "14":

    

                plt.plot(   x_plot[i], y_plot[i][0], color=(i/(len(x_plot)),1-i/(len(x_plot)),1-i/(len(x_plot))),linestyle=':',lw=0.5,zorder=1)  #画线

    

        x_plot_=[]

        for i in x_plot:

            x_plot_.extend(i)

              

        y_plot_=[]

        y_plot_text=[]

        for i in y_plot:

            y_plot_.extend(i[0])

            if show_num==False:

                y_plot_text.extend(i[1])

            #print([i[1][0]+i[2][0],i[1][1]+i[2][1]])

            if show_num==True:

                y_plot_text.extend([i[1][0]+i[2][0],i[1][1]+i[2][1]])

                

            

        for i in range(len(x_plot_)):

            #print(y_plot_text[i])

            plt.text(   x_plot_[i], y_plot_[i], y_plot_text[i],color='black',fontsize=15,ha='center', va='center',fontdict={'family' : 'Times New Roman'}  )#va='bottom',

    

    

            for count ,  ATOM_type in  enumerate(ATOM_type_list):

                if  y_plot_text[i]  == ATOM_type[1]:

                    

                    plt.scatter(x_plot_[i], y_plot_[i], color=colors_list[count],marker='o', edgecolors='black', s=200,zorder=2) #画点

    

            

    

        ax.set_axis_off()

        plt.show()

    

    

    #原子类型字典,每一层的原子序数，每两层之间的连接关系nets,序数从每一层开始转到以整个分子开始

    

    

    def str2nets(self,mol_str,order_dir,show=False):

        

            

        atom_type = {}

        layers=[[1]]

        

        

        mol_str_list = mol_str_list = mol_str.split(" ")

        m=self.txt2matrix(mol_str_list[0],order_dir)
        #print(">>>",mol_str_list[0],m)

        

        #得到所有原子序数及对应的原子类型字典atom_type，得到每一层的原子layers

        atom_type["1"]= (self.matrix2layer(m.reshape((1, 15, 15, 3)))[1][0][0][1],)  #初始原子为C

        #print(self.matrix2layer(m.reshape((1, 12, 12, 2)))[1][0][0][1])

        atom_order = 2

        for i in mol_str_list:

            m=self.txt2matrix(i,order_dir)

            layer=[]
            LAYERs = self.matrix2layer(m.reshape((1, 15, 15, 3)))[1]
            if "IndexError" not in LAYERs:
                
                LAYERs_1 = LAYERs[1]
            else:
                return "IndexError","IndexError","IndexError","IndexError"
            
            for j in LAYERs_1:

    

                atom_type[str(atom_order)]=(j[1],)

                

    

                layer.append(atom_order)

    

                atom_order+=1

                

            layers.append(layer)

                

    

        #print("atom_type",atom_type)

        #print("layers",layers,len(layers))

        

    

    

    

    

        Layers_band_type_dir = {}

        #得到每两层之间的连接关系nets

        #print("---------------------------")   

        nets=[] 

        for m_, i in enumerate(mol_str.split(" ")):

            m=self.txt2matrix(i,order_dir)

            #print(self.matrix2layer(m.reshape((1, 12, 12, 2)))[0])

            #print(m_)

            NETs,_,layer_bannd_type_dir =self.matrix2layer(m.reshape((1, 15, 15, 3)))

     

            for n_, j in enumerate(NETs):

                #print(j)

                

                net=[]

                

                for l_,k in enumerate(j):

                    #print(k)

                    try :

                        Layers_band_type_dir [ str(layers[m_][k[0][0]])+"-"+str(layers[m_+1][ k[1][0] ])   ]= layer_bannd_type_dir[   str(k[0][0])+"-"+str( k[1][0])    ]

                    except IndexError:

                        #print("next mol")



                        return  "IndexError","IndexError","IndexError","IndexError"

                            

                    

                    net.append( [  layers[m_][k[0][0]],  layers[m_+1][ k[1][0] ]   ])  #前一层layers[m_][k[0][0]]原子序数  后一层layers[m_+1][k[1][0]]原子序数

                    atom_order+=1

                #print("--")

                nets.append(net)

                #print("---------------")

        #print("nets",nets)

        if show ==True:

            self.str_plot_nets(atom_type,layers,nets,Layers_band_type_dir)  #layers是b ，nets是nodes
        #print("---------------------")
        #print("atom_type",atom_type )
        #print("layers",layers )
        #print("nets",nets )
        #print("Layers_band_type_dir ",Layers_band_type_dir )
        #print("---------------------")

        return atom_type,layers,nets,Layers_band_type_dir 

    

    

    def write_fake_mol(self,mol_name,mol_encode_str,order_dir):

        """
        分子文件创建
        """
        atoms_type,_,_,bands_type =  self.str2nets(mol_encode_str,order_dir)


        if atoms_type=="IndexError":

            print("IndexError")
    
            return "IndexError"


        #print("去除冗余占位信息")#去除冗余占位信息

        atoms_type_new = atoms_type.copy()
        bands_type_new = bands_type.copy()
        remove_hodels = []
        hodels_link=[]
        replace_index_atoms = []
        all_atoms_num = len(atoms_type_new)
        for key ,value in atoms_type_new.items():
            hodel_link = []
            if value[0][0] =="h":
                #
                remove_hodels.append(key)
                for link in list(bands_type_new.keys()):
                    h_a  = link.split("-")
                    if key in h_a:
                        
                        h_a.remove(key)
                        hodel_link.append([key,h_a,bands_type_new[link]]) #占位点，与谁相连，代替什么键
                        #print("",bands_type_new[link],link )
                        del bands_type_new[link] #取代占位符
                try:
                    bands_type_new[hodel_link[0][1][0]+"-"+hodel_link[1][1][0]]=hodel_link[0][2][1] #取代占位符 -> IndexError: list index out of range
                except IndexError:
                    return "Error"
                #print("replace",hodel_link[0][1][0]+"-"+hodel_link[1][1][0]," ",hodel_link[0][2])
                hodels_link.append(hodel_link)
                replace_index_atoms.append([str(all_atoms_num),str(key)])  #键大于原子总数，值：占位符位置
                all_atoms_num-=1
            
        #去除冗余占位信息
        #print(atoms_type_new)
        #print("bands_type_new",bands_type_new)
        #print(replace_index_atoms)
        
        #从新排序
        
        atoms_key  = list(atoms_type_new.keys())
        for key in atoms_key :
            if int(key)>all_atoms_num:
                #print(all_atoms_num,str(key),dict(replace_index_atoms),dict(replace_index_atoms))
                atoms_type_new[ dict(replace_index_atoms)[str(key) ] ] = atoms_type_new[str(key) ]
                del atoms_type_new[str(key) ]
        #print("--")
        #print("atoms_type_new",atoms_type_new)
        #print("--")
        
        bands_key  = list(bands_type_new.keys())
        for key in bands_key:
            for be_replace in dict(replace_index_atoms).keys():
                if be_replace in key:
                    value = bands_type_new[key]
                    del bands_type_new[key]
                    key= key.replace(be_replace,dict(replace_index_atoms)[be_replace])
                    bands_type_new[key] = value
                    
                    
            
        #print(bands_type_new)
                
      
        #print("去除冗余占位信息")#去除冗余占位信息
        

        #print(bands_type_new)##############################################################

    

        mol_info_1=[[int(x[0]),x[1][0]] for x in atoms_type_new.items()]

        mol_info_1.sort(key  =lambda x:x[0])

    

        mol_info_2=[[x[0].split("-"),x[1]] for x in bands_type_new.items()]

        mol_info_2.sort(key  =lambda x:int(x[0][0]))

        

    

    

        
        #print(mol_encode_str)

        lines = ["-".join(mol_encode_str.split(" "))+"\n",

                 "fake_mol\n","\n",

                  " "+str(len(atoms_type_new))+" "+str(len(bands_type_new))+"  0  0  0  0  0  0  0  0999 V2000\n"]

        for i in  mol_info_1:

            lines.append("    "+str(random.random()*10)[:6]+"    "+str(random.random()*10)[:6]+"    0.0000 "+i[1]+"   0  0  0  0  0  0  0  0  0  0  0  0\n")

            

        for i in  mol_info_2: #a，b为连接两端的原子序号

            #print(i[0][0],i[0][1],i[1])

            

            if len(i[0][0])==1:

                a="  "+i[0][0]

            if len(i[0][0])==2:

                a=" "+i[0][0]

            if len(i[0][0])==3:

                a=""+i[0][0]

                

            if len(i[0][1])==1:

                b="  "+i[0][1]

            if len(i[0][1])==2:

                b=" "+i[0][1]

            if len(i[0][1])==3:

                b=""+i[0][1]

    

    
              

            lines.append(a+b+"  "+i[1]+"  0\n")

            

        lines.append("M  END\n")

    

        #with open(mol_name+"_fake.mol","w") as mol_txt:

            #mol_txt.writelines(lines)

    

        if self.have_right_atom_valence(lines):   

    

            with open(mol_name+"_fake.mol","w") as mol_txt:

                mol_txt.writelines(lines)

        #print("txt")

            return lines

        else:

            #with open(mol_name+"_fake.mol","w") as mol_txt:

                #mol_txt.writelines(lines)

            print("valenceError")

            return "valenceError"

            



###############





if __name__=="__main__":




    

    moltools = MolTools()

    nets_dir_name = "nets_info_dir_NPname_3000_.plk"

    order_dir_name = "order_info_dir_NPname_3000_.plk"

    #folder = "mol_zinc_150000/"

    #files_list = os.listdir(folder)

    #files_list.sort(key=lambda x:int(x[:-4]))

    #files_list = [folder+x for x in files_list if ".mol" in x][90817:]

    #moltools.mol_to_nets_dir(files_list,nets_dir_name,key_is_order=False)

    
    with open(nets_dir_name,"rb") as f:

        nets_dir =pickle.load(f)

    #moltools.write_txt("ZINC250K.txt",nets_dir,files_list)

    #moltools.net_dir2order_dir(nets_dir_name,order_dir_name)

    #input("finish")



    #编码到分子

    moltools = MolTools()

    with open(order_dir_name,"rb") as f:

        order_dir = pickle.load(f)
        
    file_name = "1.mol"

    mol_str = "2 122 10 4 47 81 183 816 349 4 3 13 15"
    #mol_str = "2 1 37 10 4 47 81 183 835 137 1256 1311 13 15"
    mol_str ="2 1 11 672 2586 171 12 5 40 6 30 32 5 138 10 18 7 19 16 123 39 27 219 616 180"
    mol_str ="2 4 82 172 16 8 9 1 20 12 5 148 43 114 308 33 1064 608 1257 12 1 18 7 19 15"

    moltools.write_fake_mol(file_name,mol_str,order_dir)

    moltools.str2nets(mol_str,order_dir,True)
    print("finish")
    
    


    input()
    input()





    #分子到分子
    
    folder = "mol_zinc_150000/"

    files_list = os.listdir(folder)

    files_list.sort(key=lambda x:int(x[:-4]))

    files_list = [x for x in files_list if ".mol" in x] [:50]  

    moltools = MolTools()

    with open(order_dir_name,"rb") as f:

        order_dir = pickle.load(f)

    for file_name in [x  for x in files_list ][:50]:

        Net  = MolNets(folder+file_name)

        root=Net.min_x_size_root(folder+file_name)

        mol_str= moltools.mol_to_sequence(folder+file_name ,nets_dir,root)

        
        if "99999" not in mol_str:

            print(file_name,mol_str)

            moltools.write_fake_mol(file_name,mol_str,order_dir) #encode





    




