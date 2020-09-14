
G1 = {
    "A":["B","C"],
    "B":["A","C","D"],
    "C":["A","B","D","E"],
    "D":["B","C","E","F"],
    "E":["C","D"],
    "F":["D","G"],
    "G":["F","H","I"],
    "H":["G"],
    "I":["G"]
    }


G2 = {
    "A":["B","C"],
    "B":["E","C","D","A"],
    "C":["T","B","F"],
    "D":["B","G","E"],
    "E":["B","D","H","I"],
    "F":["C","I","J"],
    "G":["D"],
    "H":["E"],
    "I":["E","F","J"],
    "J":["F","I"]
    }


G3 = {
    "A":["B","C"],
    "B":["A","F","D","E"],
    "C":["A","G","H","I"],
    "D":["B","J"],
    "E":["B","K","P"],
    "F":["B","L"],
    "G":["C","M"],
    "H":["C","N","P"],
    "I":["C","O"],
    "J":["D"],
    "K":["E"],
    "L":["F"],
    "M":["G"],
    "N":["H"],
    "O":["I"],
    "P":["E","H"]
    
    }



G4 = {
    
    "A":["T"],
    "T":["B","C","D","A"],
    "B":["T","E"],
    "C":["T","F"],
    "D":["E","G","T"],
    "E":["B","D"],
    "F":["C"],
    "G":["D"],


    }





G5 = {
    
    "A":["T"],
    "T":["B","C","D","A"],
    "B":["T","E","F"],
    "C":["T","F"],
    "D":["E","G","T"],
    "E":["B","D"],
    "F":["C","B"],
    "G":["D"],

    }


G6={
    '1':['2','3'],
    '2':['1','4','5','6','7','8','9','10','19'],
    '3':['1','17','11'],
    '4':['2','12'],
    '5':['2','13'],
    '6':['2','12','14'],
    '7':['2','13'],
    '8':['2','14'],
    '9':['2','16'],
    '10':['2','15'],
    '11':['3','16'],
    '12':['4','6'],
    '13':['5','7'],
    '14':['6','8'],
    '15':['10'],
    '16':['11','9'],
    '17':["3","18"],
    '18':['17'],
    '19':['2']

    }


G7 = {


    '1':['2','3','4'],
    '2':['1','5','6','7'],
    '3':['1','8'],
    '4':['1','9'],
    '5':['2','12'],
    '6':['2','11'],
    '7':['2','10'],
    '8':['3','11'],
    '9':['4','12'],
    '10':['7'],
    '11':['6','8'],
    '12':['5','9']



    }





vertex_info = {
    "A":0.1,
    "B":0.2,
    "C":0.3,
    "D":0.4,
    "E":0.5,
    "F":0.6,
    "G":0.7,
    "H":0.8,
    "I":0.9,
    "J":1.0
    }


def BFS(graph,s):
    queue =[]
    queue.append(s)
    seen = []
    seen.append(s)
    while (len(queue)>0):
        vertex = queue.pop(0)
        nodes = graph[vertex]
        for w in nodes:
            if w not in seen:
                queue.append(w)
                seen.append(w)
    return seen




def BFS_layered(graph,root):
    queue =[]
    queue.append(root)
    seen = []
    seen.append(root)
    #--------------------------
    count  = 0               
    layers = [[root]]
    layer = []
    layer_index = 0
    #-------------------------- 
    while len(queue)>0:
        vertex = queue.pop(0)      
        vertex_neighbor_nodes = graph[vertex]
        for node in vertex_neighbor_nodes:
            if node not in seen:
                queue.append(node)
                seen.append(node)
        #--------------------------        
        if count>0:#except root
            layer_neighbor_nodes = []
            for layer_node in  layers[layer_index]:
                layer_neighbor_nodes .extend(graph[layer_node])
            if vertex in layer_neighbor_nodes:
                layer.append(vertex)
            else :
                layers.append(layer)
                layer = [vertex]
                layer_index+=1   
        count+=1
        #--------------------------          
    layers.append(layer)#append last layer   
    return layers



def Layered_(graph,root):
    """
    Algorithm 1
    """
    L = [None]
    width_list  = []
    
    l_1 = [root]
    L.append(l_1)
    width_list.append(len(l_1))
    
    l_2 = graph[root] 
    L.append(l_2)
    width_list.append(len(l_2))
    
    n = 3

    while True:
        
        s = set()
        for v in L[n-1]:
            s.update(graph[v])
            
        l_n = [ v for v in s if v not in L[n-1] and v not in L[n-2]]
        
        if l_n == [] :
            break
        
        L.append(l_n)
        width_list.append(len(l_n))
        
        
        n +=1
        
    L.pop(0)
    Widest_width = max(width_list)

    return L,Widest_width







def Layered(graph,root):
    """
    Algorithm 2 placeholder
    """
    L = []
    #L_inner_connect_vs = [[]]
    width_list  = []
    
    l_1 = [root]
    L.append(l_1)
    width_list.append(len(l_1))
    
    l_2 = graph[root] 
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
                connect_pairs_not_repeat.append(pair)   
        placeholder_vs = [ "-".join( list(pair)) for pair in connect_pairs_not_repeat]
        l_n.extend(placeholder_vs)
        ###########add placeholder
        
        if l_n == [] :
            break

        ########### sort layer vertex by order,avoid edge intersection
        vs_order_dir = {}
        
        for v in l_n:
            if "-" not in v:
                upper_layer_connect_vs_index =[L[n-1].index(v_c) for v_c in graph[v] if v_c in L[n-1]] #get v upper layer connect vs index
                v_order = sum(upper_layer_connect_vs_index)/len(upper_layer_connect_vs_index ) +1#get v order
                vs_order_dir[v] =v_order
            else:
                v_h_order = sum([ L[n-1].index(v_c) for v_c in v.split("-") ])/2
                vs_order_dir[v] = v_h_order

        
        l_n.sort(key = lambda x :vs_order_dir[x])
        
        vs_order_dir_copy = vs_order_dir.copy()
        
        print(l_n)
        print(vs_order_dir_copy)
        

        #同层顶点，如果它们在下一层有公共连接顶点,则在不影响order下靠近
        


        v_link_n = []
        for index_v, v in enumerate(l_n):
            if "-" not in v   :
                connect_by_next_l_vs = [v_c_c for v_c in graph[v] if v_c not in L[n-1] and v_c not in l_n for v_c_c in graph[v_c] if v_c_c in l_n and v_c_c!=v]#'\./'        
                v_link_n.append(len(connect_by_next_l_vs))
        l_n = list(zip(l_n,v_link_n))
        print("sort",l_n)
        l_n.sort( key = lambda x:-x[1])#从大到小，从以间接连接数排序
        l_n = list(zip(*l_n))[0]
        l_n = list(l_n)
        #input()
        

                
        for index_v, v in enumerate(l_n):
            
            if "-" not in v :
                
                connect_by_next_l_vs = [v_c_c for v_c in graph[v] if v_c not in L[n-1] and v_c not in l_n for v_c_c in graph[v_c] if v_c_c in l_n and v_c_c!=v]#'\./'

                #点V通过下层连接到同层的点的集合为connect_by_next_l_vs，该点V与此集合的关系有三种，集合中所有点/没有;全都是;部分是/与V所在级别的点
                
                if connect_by_next_l_vs!=[]:
                    
                    
                    vs_same_order = [ v_same_order for v_same_order in l_n if vs_order_dir[v]==vs_order_dir[v_same_order] ]
                    
                    
                    if len(set(connect_by_next_l_vs) - set(vs_same_order)&set(connect_by_next_l_vs))>0:#存在

                        be_affected_by_connect_vs_order = (-l_n.index(v)+sum( [ l_n.index( connect_by_next_l_v) for  connect_by_next_l_v in connect_by_next_l_vs ] )/len(connect_by_next_l_vs))/1000
                        #print('---',be_affected_by_connect_vs_order)
                        vs_order_dir_copy[v] = vs_order_dir[v]+be_affected_by_connect_vs_order


                    #print("--------------------",v,set(vs_same_order)&set(connect_by_next_l_vs))
                    #print("set",set(vs_same_order),set(connect_by_next_l_vs))
                    #print(v not in exchanged,v , exchanged)

                    #从连接最多的开始



        print(vs_order_dir_copy)
        
        print()
        
        l_n.sort(key = lambda x : vs_order_dir_copy[x])
        
        print("-->",l_n)
     
        




                
            
            

#如果一个点连接了被改变的点，则与这个被改变的点靠近。同时它也被改变了
#如果一个点有连接但没有连接被改变的点。它和它所连接的点彼此靠近



        for index_v, v in enumerate(l_n):
            
            if "-" not in v  and vs_order_dir_copy[v]%1==0 :
                connect_by_next_l_vs = [v_c_c for v_c in graph[v] if v_c not in L[n-1] and v_c not in l_n for v_c_c in graph[v_c] if v_c_c in l_n and v_c_c!=v]#'\./'

                #点V通过下层连接到同层的点的集合为connect_by_next_l_vs，该点V与此集合的关系有三种，集合中所有点/没有;全都是;部分是/与V所在级别的点
                
                if connect_by_next_l_vs!=[]:
                    
                    vs_same_order = [ v_same_order for v_same_order in l_n if vs_order_dir[v]==vs_order_dir[v_same_order] ]
                    #print("--------------------",v,set(vs_same_order)&set(connect_by_next_l_vs))
                    
                    
                    if len(set(vs_same_order)&set(connect_by_next_l_vs))==len(set(connect_by_next_l_vs)) : #间接连接点全都是同级别点
                        #mode .append(2)

                        
                        vs_connect_exchanged = [ v_connect_exchanged  for v_connect_exchanged  in  connect_by_next_l_vs if vs_order_dir_copy[v_connect_exchanged ]%1!=0]

                        
                        if len(vs_connect_exchanged)!=0 :
                        #如果一个点连接了被改变的点，则与这个被改变的点靠近。同时它order也被改变了
                            pass
                            #vs_order_dir_copy[v]  = [ v_order_exchanged  for v_order_exchanged in  connect_by_next_l_vs if vs_order_dir_copy[v_order_int]%1!=0]

                        if len(vs_connect_exchanged)==0 :
                        #如果一个点有连接但没有连接被改变的点。它和它所连接的点彼此靠近

                        
                            vs_order_dir_copy[v] = vs_order_dir_copy[v]+(1+l_n.index(v))/100000
                            #print(v,vs_order_dir_copy[v],l_n.index(v))
                            #保持这个中心顶点原先位置不变
                        
                            for connect_by_next_l_v in connect_by_next_l_vs:

                                vs_order_dir_copy[connect_by_next_l_v]=  vs_order_dir_copy[v]+( 1+l_n.index(connect_by_next_l_v)  -  l_n.index(v)  )/10000000
                                #被连接的点放在中心顶点附近




                
                
                
        print(vs_order_dir_copy)
        

        
        l_n.sort(key = lambda x : vs_order_dir_copy[x])
        

        print("-->>",l_n)

        input()
        
        ########### sort layer vertex by order,avoid edge intersection
        

        L.append(l_n)
        width_list.append(len(l_n))
        
        n +=1
        
    Widest_width = max(width_list)

    return L,Widest_width




def Bipartite_nets(graph,L):
    nets = []
    for index in range(len(L)):
        net = []
        if index>0:
            print(L[index-1],L[index])
            for index_v, v in enumerate(L[index-1]):
                if "-" not in v:
                    for index_v_c,v_c in enumerate(L[index]):
                        if "-" not in v_c  :
                            if v_c in graph[v]:
                                net.append([v,v_c])
                                #print(index_v,index_v_c)
                else:
                    net.extend( [ [x,v] for x in  v.split("-")] )
            nets.append(net)
            
    return nets
        


    
def Longest_L_root(graph):
    """
    Get longest layer. Input graph, output longest_layer,longest_root
    """
    L_info_dict = {}

    for root in graph.keys():
        L,width = Layered(graph,root)
        L_info_dict[root]=(len(L),width)


    L_info_list  = list(L_info_dict.items())
    L_info_list.sort( key =lambda x: x[1][1]) #narrowest
    L_info_list.sort( key =lambda x:-x[1][0]) #longest

    L_longest_narrowest_root   = L_info_list[0][0]
    L_longest_narrowest_length = L_info_list[0][1][0]
    L_longest_narrowest_width  = L_info_list[0][1][1]
    
    print(L_info_list)
    
    return  L_longest_narrowest_root,L_longest_narrowest_length,L_longest_narrowest_width

            



if __name__ =="__main__":

    graph = G6
    #root = Longest_L_root(graph)[0]
    root = '1'
    L = Layered(graph,root)[0]
    print(L)
    print(Bipartite_nets(graph,L))
    
    

    #for root in graph.keys():
        #print(root," ",BFS_layered(graph,root))
        #print(root," ",Layered(graph,root)[0])
        #input("----")











