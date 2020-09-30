import numpy as np

#A = np.random.randint(0,2,(5,5))

a1 = np.array([[1,0,0],
              [ 1,0,0],
              [ 0,0,0]])

a2 = np.array([[1,1,0],
              [ 0,1,0],
              [ 0,0,0]])

a3 = np.array([[1,0,0],
              [ 0,1,0],
              [ 0,1,0]])


x_1 = np.zeros([3,3,8])
x_1[0][0][0]  = 1
x_1[1][0][1]  = 1 

x_2 = np.zeros([3,3,8])
x_2[0][0][2]  = 1
x_2[0][1][3]  = 1
x_2[1][1][4]  = 1

x_3 = np.zeros([3,3,8])
x_2[0][0][5]  = 1
x_2[1][1][6]  = 1
x_2[2][1][7]  = 1


def F(a,b):
    

    #print(a[0],a[1],a[2])
    #print(b[0],b[1],b[2])
    #b = b.T
    W0=np.ones([3,])
    W1=np.ones([3,])
    W2=np.ones([3,])
    W3=np.ones([3,])
    W4=np.ones([3,])
    W5=np.ones([3,])
    
    c = np.zeros([3,3])
    cat = (b[0]*W0+b[1]*W1+b[2]*W2)
    print(b[0].shape,a[:,0].shape)
    c[:,0] =a.T[0]*W3+ cat
    c[:,1] =a.T[1]*W4+ cat
    c[:,2] =a.T[2]*W5+ cat
    #print( (b[0]+b[1]+b[2])*0.1)


    print(a*c)
    return  a*c

    
    

F(a2,a3)
print()
F(a1,a2)

print()

an = F(a2,a3)
F(a1,an)




