import os
import time
import re
from collections import Counter
import pickle
import random
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors
from PIL import Image
from collections import OrderedDict

from copy import deepcopy
#from mpl_toolkits.axes_grid1 import make_axes_locatable
###############
#6.6生成错误处理#
###############
class Gragh():
    def __init__(self,nodes,sides):
        '''
        nodes 表示点
        sides 表示边

        '''
        # self.sequense是字典，key是点，value是与key相连接的点
        self.sequense = {}
        # self.side是临时变量，主要用于保存与指定点相连接的点
        self.side=[]
        for node in nodes:
            for side in sides:
                u,v=side
                # 指定点与另一个点在同一个边中，则说明这个点与指定点是相连接的点，则需要将这个点放到self.side中
                if node ==u:
                    self.side.append(v)
                elif node == v:
                    self.side.append(u)
            self.sequense[node] = self.side
            self.side=[]
        #print ("sequense",self.sequense)
        #input()




    '''
     breadth-First-Search
     BFS是从根节点开始，沿着树的宽度遍历树的节点。如果所有节点均被访问，则算法中止。
           广度优先搜索的实现一般采用open-closed表。
    '''
    def BFS(self,node0):
        #queue本质上是堆栈，用来存放需要进行遍历的数据
        #order里面存放的是具体的访问路径
        queue,order = [],[]
        #首先将初始遍历的节点放到queue中，表示将要从这个点开始遍历
        # 由于是广度优先，也就是先访问初始节点的所有的子节点，所以可以
        queue.append(node0)
        order.append(node0)


        #elements = []
        while queue:
            #queue.pop(0)意味着是队列的方式出元素，就是先进先出，而下面的for循环将节点v的所有子节点
            #放到queue中，所以queue.pop(0)就实现了每次访问都是先将元素的子节点访问完毕，而不是优先叶子节点
            #print(queue,"-",order)

            

            #print("queue",queue)
            
            v = queue.pop(0)
            #print("v:",v)
            
            #element=[]
            
            
            for w in self.sequense[v]:
                
                
                
                if w not in order:
                    #print("-----------",w)
                    # 这里可以直接order.append(w) 因为广度优先就是先访问节点的所有下级子节点，所以可以
                    # 将self.sequense[v]的值直接全部先给到order
                    order.append(w)
                    queue.append(w)
                    

            

                    
                    #element.append(w)
            #elements.append(element)

                    
            
                 
            
                
        return order


###############
#
###############

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
    (0.12,"h2"),    #数值大小也决定它在同层中的位置
    (0.13,"h3"),    #数值大小也决定它在同层中的位置
    (0.14,"h4"),    #数值大小也决定它在同层中的位置
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
            if len(line)==7 :
                if "CHG" in line:
                    break
                        
                atom_bond .append(line[:3])
            
            if len(line)==6 :   #超过99个分子的mol文本粘连问题
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






def get_C_list(file_name):      #为了防止初始结构种类太多，根节点从C出发
    
    C_list=[]
    atom_type,pos ,band  = get_mol_info(getLineElement(file_name))
    for i,j in atom_type.items():
        if j[0]=="C":
            C_list.append(int(i))
    C_list.sort()      
    return C_list


def have_other_type(file_name):
    atom_type,pos ,band  = get_mol_info(getLineElement(file_name))
    for i,j in atom_type.items():
        #print(i,j)
        if j[0]=="Other":
            return True

    return False
    


def branch(band,core):    #找到某个原子周围相键连的所有原子，core中心原子
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



def get_nodes_sides(atoms_num,band):
    nodes = [i+1 for i in range(atoms_num)]
    
    sides = []
    
    for i in range(1,atoms_num+1):
        #print(i,branch(band,i))
        for j in branch(band,i):
            sides.append([i,j])
    #sides = list(set(sides ))
    sides_ =[]  #剔除重复

    for i in sides:
        if (i not in sides_) and  (i.reverse() not in sides):
            sides_.append(i)
    
    ##print(sides_)
    
    return nodes,sides_


#root=1

##########
###############
#file_name = "088602.mol"
def get_valence(file_name):
    atom_type,pos ,band  = get_mol_info(getLineElement(file_name))

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
        for j in branch(band,i):
            if int(i)<int(j):
                band_name = str(i)+"-"+str(j)
            else:
                band_name = str(j)+"-"+str(i)
            
            atom_valence+=int(band_type_dir[band_name])
            #print(band_type_dir[band_name],end = " ")
        #print(atom_valence)
        atom_valence_list.append(  [str(i) , [atom_type[str(i)][0],atom_valence]]  )
                                        
    #for i in band :
        
        #print(i)
    return atom_valence_list


            



def have_right_atom_valence(file_name):
    atom_valence_list = get_valence(file_name)    
    for i in atom_valence_list:
        #print(i[1][0],i[1][1])
        if i[1][1] not in ATOM_valence[i[1][0]]:
            return False
    return True
#print(have_right_atom_valence(file_name))      
#input()
##########

def mol2net(file_name,root):
    atom_type,pos,band = get_mol_info(getLineElement(file_name))
    all_atom_num = len(pos)
    #print("all_atom_num",all_atom_num )
    #band = get_mol_info(getLineElement(file_name))[2]       

    nodes,sides = get_nodes_sides(len(pos),band)  #得到所有的点和边
    G = Gragh(nodes,sides)                        
    ##print()
    ##print()
    ##print ("DFS",G.DFS(root))
    ##print()
    ##print ("BFS",G.BFS(root))

    NODES = G.BFS(root)               #输出为遍历的点序列
    #print(">>",NODES)
    NODES.pop(0)
    NODES.append(0)
    #print("NODES",NODES)
    
    L=[]
    L.append([root])
    l=[]
    L_last=0

    ######
    placeholder_dir ={}
    ######
    #L_next=0
    placeholder_count =0
                    
    for count , node in enumerate(NODES):

        
        L_last_branch=[]  #获得上一层节点产生的所有分支节点    L_last负责计层数
        for i in L[L_last]:
            L_last_branch.extend(branch(band,i))
            #print("#",L_last_branch)
        if L_last>0:  #从第二层开始删除L_last_branch上上层的节点，这是由于branch（）不分方向
            
            for i in L[L_last-1]:
                if i in L_last_branch:
                    L_last_branch.remove(i)
            
        
        
        if node not in L_last_branch:  #node 如果不在上一层的下一分支中，则该node不属于该层
            L_last += 1
            L.append(l)
            #print(">l",l)
            l=[node]    #node作为新一层的开始
            #break
        else :
            l.append(node)  #node 在上一层产生的分支节点当中，说明这一层还没结束
            #print("--l",l)
           
            
            #print("l",l)
            #检查同层成键  在不断追加node 时，我们需要看看新追加到该层的原子有没有内部相连
            
            for m in l:
                for n in branch(band,m):
                    if node ==n:
                        #print("yes,get it!",L_last,"     ",node,m)
                        all_atom_num += 1
                        new_holder = all_atom_num
                        placeholder_count+=1
                        #print("placeholder_count---------------",placeholder_count)
                        placeholder_dir[str(L_last+1)+"-"+str(placeholder_count)] = [[node,new_holder],[m,new_holder]]  #创建字典保存新的placeholder键连关系

                    
                
                

        #print( ">L:",L[L_last],"last---------------------->",L_last)
    # L为每一层元素，原始顺序
    

    Layer = L  #删除i_重复


    #if show==True:
        #for i in Layer :
        
            #print(">>>Layer",i)

    Layer.append([])
    

    Layer_long = len(Layer)

    Layers_net = []

    for i in range(Layer_long):

        Layer_net = []
    
        for j in Layer[i]:
            for k in Layer[i+1]:
                if k in branch(band,j):
                    Layer_net.append([j,k])

        Layers_net.append(Layer_net)


    for i in range(Layers_net.count([])):         
        if [] in Layers_net:
            Layers_net.remove([])


#################################不记placeholder的边宽与图长
    x_size = []
    for i in Layers_net:
        #print(i)
        x_size.append(len(i))
    Layers_long=len(Layers_net)


    
        
    ########增加placeholder键连关系
            
    for key,item in placeholder_dir.items():
        for i in item:

            
            
            #print(int(key.split("-")[0]),Layers_long)
            
            if int(key.split("-")[0])>=Layers_long:  #如果placeholder在最后一层
                Layers_long+=1
                Layers_net.append([])
            

            
            Layers_net[int(key.split("-")[0])].append(i)
    #####################################
    


    #print("Layers_net",Layers_net)
    
    #print("placeholder_dir",placeholder_dir)
    #input("-------------##")

    return  Layers_net , max(x_size) ,Layers_long,placeholder_dir

##############################################################################################33

#die = [root]
#-------------------------------------------------



 #将树表示为矩阵


def lists_order(num,lists):
    for cont , i in enumerate(lists):
        if num == i[0]:
            #print("yes")
            return lists.index(i)
        if cont ==len(lists)-1:
            return 0


def sort_atom(t):#原子已经按照原子种类排好后同种原子种类的排序
    d=[]
    l=[t[0]]
    for cont,i in enumerate(t):
        if cont>0:
            if i[1] ==t[cont-1][1]:
                l.extend([i])
            else:
                d.append(l)
                l=[i]
    d.append(l)
    d=[sorted(i,key = lambda i:i[2]) for i in d]
    r=[]
    for i in d:
        for j in i:
            r.append(j)
    return r



#Layers_net

def net2matrix(file_name,root=1,size=[20,12,12,5]):
    

    atom_type,pos ,band = get_mol_info(getLineElement(file_name))  #获得分子mol文件属性
    nodes,_,_,placeholder_dir = mol2net(file_name,root)                                 #输出为连接关系的列表，二维列表：第一维是层，第二维是分子连接关系

    
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

        
        if item[0][0]< item[1][0]:
                band_name = str(item[0][0]) +"-"+ str(item[1][0])

        else:
                band_name = str(item[1][0]) +"-"+ str(item[0][0])

                                

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
    nodes_dicts = []
    for i in nodes:
        #print("--i",i)
        nodes_dicts.append(dict([  [elem[1],elem[0]]  for elem in i]))
        
  
#########################################    对于不同元素的排序

    for cont, i in enumerate (nodes):               #遍历连接关系列表

        y=[]
        
        for j in i:

            #print(j[1])

            y.append((j[1],atom_type[str(j[1])][1] ))
        
        y=list(set(y))                              #去除重复
        
        y.sort(key=lambda x:-x[1])   #以元素周期表排序
                   
        #y.sort(reverse=True)
        
        #if cont >-1: #从第二层开始对同种元素交换数序
        y_order = []   #单层包含的节点

        same_atom_type=[]
###############################################         对于同种元素的排序   
        for y_elem in y:



                
                
            last_link_node = nodes_dicts[cont][y_elem[0]]#找到与该顶点连接的上一个顶点，如果多个相连就取最后一个
                #print(b[cont],last_link_node)
                
            order = b[cont].index(last_link_node)  #按照上一个顶点在上一层中的位置来划分该顶点在此位置上的排序等级
                #print(y_elem[0],'-',order ,end=" ")
            y_order.append( [y_elem[0],atom_type[str(y_elem[0])][1],order])

                
            
            #y_order .sort(key=lambda x:x[1])
            
            
            #b_order.append(list(zip(*y_order))[0])
            #print("y_order",y_order)
            #print("y_order_",sort_atom(y_order))
            #print("b",b[cont])
        y_order_0= list(zip(*   sort_atom(y_order)   ))[0]
            #b_order.append(y_order_0)
##########################################################################
                


                
            #print("y_elem:")
        
        #b_order [0]=[root]
        #print("b_order",b_order)          #
        #y=list(zip(*y))[0]  #取每一对的第一个元素
        #print("y",y)
        b.append(y_order_0)                               #每层有顺序的原子列表
        

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


                            
                            chain[cont,m,n] = [ atom_type[ str( j[0] ) ][1]  ,atom_type[ str(j[1]) ][1]  ,BANDencode[band_type_dir[band_name]]   ]
                            #print(atom_type[ str( j[0] ) ][1]  ,atom_type[ str(j[1]) ][1]  ,band_type_encode)




                            
                        if len(size)==4 and size[3]==5 :
                            p =   np.array(pos[j[1]-1])  -  np.array(pos[ j[0]-1] )
                            chain[cont,m,n] = [ atom_type[ str( j[0] ) ][1]  ,atom_type[ str(j[1]) ][1],p[0],p[1],p[2]     ]                  


    return chain ,len(nodes),b  #矩阵 长度  层顶点







def plot_nets(file_name,root,size,show_num=False,show_info=False):


    atom_type,pos ,band = get_mol_info(getLineElement(file_name))  #获得分子mol文件属性
    
    nodes,_,_,placeholder_dir = mol2net(file_name,root)                                 #输出为连接关系的列表，二维列表：第一维是层，第二维是分子连接关系


    
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

        
        if item[0][0]< item[1][0]:
                band_name = str(item[0][0]) +"-"+ str(item[1][0])

        else:
                band_name = str(item[1][0]) +"-"+ str(item[0][0])

                                

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

            
    
    b= net2matrix(file_name,root,size)[2]#########有顺序的每层顶点


        
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
        
#plot_plot_plot_plot_plot_plot_plot_plot_plot_plot_plot_plot_plot_





def plot_matrix(file_name,root,size  ):

    
    
    matrix,layer_num,b= net2matrix(file_name,root,size)  #矩阵，长度，层顶点

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

def renove_repeat(List):
    list_new=[]
    for i in List:
        if i not in list_new:
            list_new.append(i)
    return list_new







def min_x_size_root(file_name):
    
    roots = get_C_list(file_name)
    root_size_long=[]
    for i in roots:
        nodes ,size, long,_ = mol2net(file_name,i)
        
#        print(i, size, long)  #min size ,long
        root_size_long.append([i, size, long])
    nice_long_root =sorted(root_size_long ,key=lambda x :-x[2])     #取最长的序列对应根节点
    nice_size_root =sorted(nice_long_root ,key=lambda x :x[1])    #取最窄的“序列最大宽度”对应根节点
    #print(file_name,"nice_size_root",nice_size_root)
    #print(nice_size_root[0])
    nice_root = nice_size_root[0][0]
    return nice_root
    

def min_x_size_roots(file_name):
    
    roots = get_C_list(file_name)
    root_size_long=[]
    for i in roots:
        nodes ,size, long,_ = mol2net(file_name,i)
        
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





def build_train_test_data(save_name):
    with open("gdb9.sdf_gdb_129135removed_.csv","r") as csvfile:
        reader=csv.reader(csvfile)
        floder="SDF/"
        sdf_E = list( zip(  sorted(list(os.listdir(floder)) ,key=lambda x:int(x[:-4]))[0:133883]  ,   list(reader)[1:133884]) )

    random.seed(0)
    random.shuffle(sdf_E)
    #for i in sdf_E:
        #print(i[0],i[1][0])


    sdf_data_x=[]
    sdf_data_y=[]
    sdf_data_L=[]
    count=0

    dx = int(len(sdf_E)/5)

    cat_num = -1
    
    for i in range(5):

        cat_num+=1
        start = i*dx+cat_num
        end = (i+1)*dx+cat_num
        if end == 133880:
            end = 133883
        
        for sdf , E in sdf_E[ start :  end  ] :
        #print(floder+sdf,E[1:])
        #print(count)
        #input()
            count+=1
         
            if len(get_C_list(floder+sdf )) > 0:
                root=min_x_size_root(floder+sdf)
                L = mol2net(floder+sdf,root)[0]
                n = net2matrix(floder + sdf ,root,[10,12,12,3])[0]


                sdf_data_L.append( sdf )

                sdf_data_x.append( n.reshape(12*10,12,3) )

                sdf_data_y.append( np.array([float(e) for e in E[1:] ]) )

                #print(np.array([float(e) for e in E[1:] ]))
                
                print(count,sdf)


        #print(len(sdf_data_L))

        sdf_data_train = [  sdf_data_L ,  sdf_data_x , sdf_data_y]
    
        #print(len(sdf_data_train))

        with open(save_name+str(i),"wb") as f :
        
            pickle.dump(sdf_data_train,f)
            


#save_name = "sdf_data_train" 
#start_cut = 0
#end_cut = 70000
#build_train_test_data(start_cut,end_cut ,save_name)




#-----------------------------------------------remove   max(x_size)>10
def remove(files_name,xy_size):
#files_name="train_3-0/"
    file_dir = os.listdir(files_name)
    count=0
    for i in  file_dir:
        if "mol" in i:
            #L = mol2net(files_name+i,root)[0]
            #print(i)
            root=min_x_size_root(files_name+i)
            if mol2net(files_name+i ,root)[1]>xy_size:#len(L)+1 多一层表示隔断
                os.remove(files_name+i)
                count+=1
                print(count,": ",i,mol2net(files_name+i ,root)[1])
      

#files_name="train_3/"
#to_files_name="net2matrix__3/"
#remove(files_name,16)
#input()

def mol_to_sequence(file_name,nets_dir_name,root=0):
    sequence=[]
    if root==0:
        root=min_x_size_root(file_name)
    L,w_num = mol2net(file_name,root)[:2]

    if w_num<=12:
        mol =net2matrix(file_name ,root,[len(L)+1,12,12,3])[0]
        #mol = np.load(to_files_name+file)
        with  open(nets_dir_name,"rb") as f:
            s =pickle.load(f)
        for i in mol :
            if str(i) in s.keys():
                sequence.append(s[str(i)][1])
                #if s[str(i)][1] >8282:
                    #print(s[str(i)][1],end="#")
                #else:
                    #print(s[str(i)][1],end=" ")

            else:
                #print("99999",end=" ")
                sequence.append(99999)
            
    return sequence

#files_name = "SDF/"
#for i in range(1,5):
    
#    print(i,mol_to_sequence(files_name+str(i)+".sdf","nets_info_dir_NPname_3000.plk"))
#input("######################################")
#######################################


######################################################################################
def mol_to_nets_dir(folder,dir_save_name,key_is_order=False):
    file_dir = os.listdir(folder)
    file_dir.sort(key=lambda x:int(x[:-4]))
    file_dir = file_dir[:]
    


#b= np.load( folder+"1..matrix.npy" )
    nets_str = []
    nets_dir = {}
    count = 0
    for   i in file_dir:
        #nets = np.load(folder + i )
        #print(get_C_list(folder+i))
        if len(get_C_list(folder+i)) >0:

            
            
            root=min_x_size_root(folder+i)
            #root = root_info[cont][0]
            #count+=1

            
            print(i,root)
            L,w_num = mol2net(folder+i,root)[:2]

            if w_num<=12:
        
                nets =net2matrix(folder+i ,root,[len(L)+1,12,12,3])[0]
                for net in nets:
                    nets_str .append(str(net))
                    #print(str(net))
                    nets_dir [str(net)] = net

            
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
            nets_info_dir[i[0] ]=[ nets_dir[i[0]],count ,i[1]]   #{ 矩阵名字：[矩阵，排序，出现次数] }
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
        






#print("-----------------------------------------数据库转化为文本")
def write_txt(TXT_name,folder,nets_dir_name):
    
    files_list = os .listdir(folder)



    
    
    #files_list.sort(key=lambda x:int(x[:-4]))
    #file_dir = file_dir[:]
    
    with open(TXT_name,"w") as TXT :
        cont=0
        for file_name in files_list:
            #root = root_info[cont][0]
            

            if len(get_C_list(folder+file_name)) >0:
                root=min_x_size_root(folder+file_name)
                print( file_name,root)
                cont+=1
                TXT.writelines(  str(x)+" " for x in    mol_to_sequence(folder+file_name,nets_dir_name,root)  )
                TXT.writelines("    "+file_name+"\n")















#key为矩阵的字典转化成key为数序的字典



def net_dir2order_dir(nets_dir_name,order_dir_name):

    
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
def matrix2layer(M):  #(13, 12, 12, 3)
    M=M.reshape(1,12,12,3)
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
 

        l1= list(set(list(zip(*net))[0]))
        l2= list(set(list(zip(*net))[1]))
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


def txt2matrix(mol_str,order_dir_name):

    with open(order_dir_name,"rb") as f:
        
        order_dir = pickle.load(f)

        #print(order_dir[int(mol_str)].shape)
        return order_dir[int(mol_str)]
        


#nets_dir_name = "nets_info_dir_NPname_3000.plk"
#order_dir_name = "order_info_dir_NPname_3000.plk"

#print("txt2matrix",txt2matrix("2",order_dir_name))
#print("matrix2layer(M)",matrix2layer(txt2matrix("2",order_dir_name)))
#input()

def str_plot_nets(atom_type,b,nodes,band_type_dir,show_num=False):



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


def str2nets(mol_str,order_dir_name,show=False):
    
        
    atom_type = {}
    layers=[[1]]
    
    
    mol_str_list = mol_str.split(" ")
    m=txt2matrix(mol_str_list[0],order_dir_name)
    
    #得到所有原子序数及对应的原子类型字典atom_type，得到每一层的原子layers
    atom_type["1"]= (matrix2layer(m.reshape((1, 12, 12, 3)))[1][0][0][1],)  #初始原子为C
    #print(matrix2layer(m.reshape((1, 12, 12, 2)))[1][0][0][1])
    atom_order = 2
    for i in mol_str_list:
        m=txt2matrix(i,order_dir_name)
        layer=[]
        for j in matrix2layer(m.reshape((1, 12, 12, 3)))[1][1]:

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
        m=txt2matrix(i,order_dir_name)
        #print(matrix2layer(m.reshape((1, 12, 12, 2)))[0])
        #print(m_)
        layer_bannd_type_dir =matrix2layer(m.reshape((1, 12, 12, 3)))[2]
 
        for n_, j in enumerate(matrix2layer(m.reshape((1, 12, 12, 3)))[0]):
            #print(j)
            
            net=[]
            
            for l_,k in enumerate(j):
                #print(k)
                try :
                    Layers_band_type_dir [ str(layers[m_][k[0][0]])+"-"+str(layers[m_+1][ k[1][0] ])   ]= layer_bannd_type_dir[   str(k[0][0])+"-"+str( k[1][0])    ]
                except IndexError:
                    print("next mol")
                    return  "IndexError","IndexError","IndexError","IndexError"
                        
                
                net.append( [  layers[m_][k[0][0]],  layers[m_+1][ k[1][0] ]   ])  #前一层layers[m_][k[0][0]]原子序数  后一层layers[m_+1][k[1][0]]原子序数
                atom_order+=1
            #print("--")
            nets.append(net)
            #print("---------------")
        
    #print("nets",nets)
    if show ==True:
        str_plot_nets(atom_type,layers,nets,Layers_band_type_dir)  #layers是b ，nets是nodes
    return atom_type,layers,nets,Layers_band_type_dir 


def write_fake_mol(mol_name,mol_encode_str,order_dir_name):

    
    atoms_type,_,_,bands_type =  str2nets(mol_encode_str,order_dir_name)
    if atoms_type=="IndexError":
        return "IndexError"

    atoms_type_new = atoms_type.copy()
    bands_type_new = bands_type.copy()

    place_holder_list  = []
    
    is_replaced_node = []
    
    for key,item in atoms_type.items():  #发现占位符
        
        if "h1" in item or "h2" in item or"h3" in item or"h4" in item :
            
            is_replaced_node.append(key)
            place_holder_list.append([key,item])
            atoms_type_new.pop(key)
            #print("key ",key)
            
        


    for holder in place_holder_list:  #寻找被取代的键
        new_band = []
        band_type = ""
        for key,item in bands_type.items():
        
            
            if holder[0] in key.split("-"):

                #print(key.split("-"),item)
                try:
                    bands_type_new.pop(key)
                except KeyError:
                    return "KeyError"
                new_band.append(key.split("-")[0])
                if item=="11":
                    band_type = "1"
                if item=="12":
                    band_type = "2"
                if item=="13":
                    band_type = "3"
                if item=="14":
                    band_type = "4"

            
                    
        
        #print(new_band,band_type)
        bands_type_new["-".join(new_band)]=band_type
        

    #print("is_replaced_node-",is_replaced_node)
    #print('all_atoms',len(atoms_type_new))

    is_replaced_node_s = []
    for i in is_replaced_node:
        if int(i)<=len(atoms_type_new):#没必要用大于原子数的数值替换
            is_replaced_node_s.append(i)
            

    #解决序数大于总原子数    
    replace_count = 0
    #print("atoms_type_new",atoms_type_new)    
    replaced_node = {} #值来替换键
    for key,item in atoms_type_new.items():  
        if int(key)>len(atoms_type_new):
            
            #print("bigger ",key)
            
                
            #print(key ,"replace ",is_replaced_node_s[replace_count] )
            atoms_type_new[ is_replaced_node_s[replace_count]  ]=item
                    
            replaced_node[key]= is_replaced_node_s[replace_count]#为下一步band的修改做准备
            atoms_type_new.pop(key)#删除大于原子数的键
            
            replace_count+=1
                    
            
    #print("atoms_type_new",atoms_type_new)
    #print("replaced_node",replaced_node)
    for r_key,r_item in replaced_node.items():
        
        for key,item in bands_type_new.items():
            nodes = key.split("-")
            if r_key in nodes:
                
                if r_key == nodes[0]:
                    #print("0node_r ",r_key)
                    bands_type_new.pop(key)
                    bands_type_new[ replaced_node[r_key]+"-"+nodes[1]]=item
                    
                
                elif r_key == nodes[1]:
                    #print("1node_r ",r_key)
                    bands_type_new.pop(key)
                    bands_type_new[ replaced_node[r_key]+"-"+nodes[0]]=item
                    


    #print(bands_type_new)

    mol_info_1=[[int(x[0]),x[1][0]] for x in atoms_type_new.items()]
    mol_info_1.sort(key  =lambda x:x[0])

    mol_info_2=[[x[0].split("-"),x[1]] for x in bands_type_new.items()]
    mol_info_2.sort(key  =lambda x:int(x[0][0]))
    


    
    lines = [mol_encode_str+"\n",
             "fake_mol\n","\n",
              " "+str(len(atoms_type_new))+" "+str(len(bands_type_new))+"  0  0  0  0  0  0  0  0999 V2000\n"]
    for i in  mol_info_1:
        lines.append("    "+str(random.random()*10)[:6]+"    "+str(random.random()*10)[:6]+"    0.0000 "+i[1]+"   0  0  0  0  0  0  0  0  0  0  0  0\n")
        
    for i in  mol_info_2:
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




            
        lines.append(a+b+"  "+i[1]+"  0  0  0  0\n")
        
    lines.append("M  END\n")
   

    with open(mol_name+"_fake.mol","w") as mol_txt:
        mol_txt.writelines(lines)
    #print("txt")
    return lines
        




if __name__=="__main__":

    nets_dir_name = "nets_info_dir_NPname_3000.plk"
    order_dir_name = "order_info_dir_NPname_3000.plk"


    #folder = "SDF/"
    #dir_save_name="nets_info_dir_NPname_3000.plk"
    #mol_to_nets_dir(folder,dir_save_name,key_is_order=False)
    #write_txt("testQM9.txt",folder,"nets_info_dir_NPname_3000.plk")
    #net_dir2order_dir(nets_dir_name,order_dir_name)
    #file_name= "creat_mol/56_fake.sdf"
    #print(min_x_size_roots(file_name))
    #root=min_x_size_root(file_name)
    #print(root)
    #root=16
    #print(mol2net(file_name,root))
    #L = mol2net(file_name,root)[0]
    #print(mol2net(file_name,root)[1:])
    #print(net2matrix(file_name,root,size=[len(L),12,12,3])[1:])
    #M=net2matrix(file_name,root,size=[len(L),12,12,3])[0]
    #print(M.shape)
    #plot_nets(file_name,root,[len(L)+1,12,12,3],False,False)
    #mol_encode_str =  " ".join( [ str(x) for x in mol_to_sequence(file_name,nets_dir_name)[:-1]] )
    #print("---------------------------------")
    #print(str2nets(mol_encode_str,order_dir_name))


    #str2nets(mol_encode_str,order_dir_name,True)

    #write_fake_mol("1",mol_encode_str,order_dir_name)




    

    #input("Finish")



    #nets_dir_name = "nets_info_dir_NPname_3000.plk"
    #order_dir_name = "order_info_dir_NPname_3000.plk"

    #mol_str = "22 9 4 117 71 94"
    
    mol_encode_str = "11 2 4 8 44 41 1 10 7 304 300"  #56
    
    #str2nets(mol_encode_str,order_dir_name,True)

    write_fake_mol("1",mol_encode_str,order_dir_name)

    #input()
    print("finish")










    input()

    input()







    file_name= "107760.sdf"


    
    #print(min_x_size_roots(file_name))
    root=min_x_size_root(file_name)
    #print(root)
    #root=1

    
    L = mol2net(file_name,root)[0]
    print(mol2net(file_name,root)[1:])
    M=net2matrix(file_name,root,size=[len(L),12,12,3])[0]
    
    print(M.shape)

    
    plot_nets(file_name,root,[len(L)+1,12,12,3],False,False)
    plot_matrix(file_name,root,[len(L),12,12,3]  )



    input()
    input()
    input()
#net_dir2order_dir(nets_dir_name,order_dir_name)
######################################################################################################
#    print("Finish")
#    with  open(dir_save_name,"rb") as f:
#        s =pickle.load(f)
    
#    count = 0 
#    for i ,j in s.items(): 
#        count +=1
        #print(i,j[1],j[2])
#        if j[1]<20:
#            print(j[1],j[2])
	    #break
#        if count >20:
#            break
#nets_dir(to_files_name,"nets_info_dir.plk",key_is_order=False)
#nets_dir(to_files_name,"nets_info_dir_NPname.plk",key_is_order=True)

#to_files_name="net2matrix+0/"
#file_dir = os.listdir(to_files_name)
#for count , file in enumerate(file_dir):
    
#    mol = np.load(to_files_name+file)
#    with  open("nets_info_dir_NPname_1.plk","rb") as f:
#        s =pickle.load(f)
#    for i in mol :
#        if s[str(i)][1] >8282:
#            print(s[str(i)][1],end="#")
#        else:
#            print(s[str(i)][1],end=" ")
#    print()
#    if count>100:
#        break

