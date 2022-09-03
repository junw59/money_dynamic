import numpy as np
import ctypes

sotest = ctypes.cdll.LoadLibrary("./pcc.so")
# sotest = ctypes.WinDLL("./pcc3.so")
Num1 = ctypes.c_int(10)
Num2 = ctypes.c_int(20)
sotest.prip(Num1,Num2)
sotest.add(Num1,Num2)

def vivj(a):
    n=len(a)
    s=0
    num=0
    for i in range(n):
        for j in range(i+1,n):
            s+=a[i]*a[j]
            num+=1
    return s/num


def calcij(a):
    l=len(a[0])
    l2=len(a)
    s1=[]
    for k in range(l2):
        suma=0
        n=0
        b=a[k]
        for i in range(l):
            for j in range(i+1,l):
                if i == j : continue
                suma += b[i]*b[j]
                n+=1
        s1.append(suma/n)
        print(k,end=' ')
    return np.mean(s1)

# sotest = ctypes.cdll.LoadLibrary("./pcc.so")
# sotest = ctypes.WinDLL("./pcc.so")
# x=[1,2,3,4.5,10, 3,4,5,6.5,7]
# n=len(x)
# cbarr=ctypes.c_float * n
# cba=cbarr()
# for i in range(n):
#     cba[i] = x[i]

# sotest.cal_vivj.restype=ctypes.c_float
# sotest.mean_cij.restype=ctypes.c_float
# cbn=ctypes.c_int(n)
# sotest.cal_vivj(cba,cbn)
# sotest.mean_cij(cba,5,2)


a=np.array([[1,2,3,4,5,6,7,8,9,10],[4,5,6,7,8,9,10,11,12,13]])
a=np.random.normal(0,1,size=(100,100))
print("pycal ",calcij(a))


sotest = ctypes.cdll.LoadLibrary("./pcc.so")
n1,n2=a.shape
x=a.flatten()
# n1,n2=v3[0].shape
# x=v3[0].flatten()

n=n1*n2
cbarr=ctypes.c_float * n
cba=cbarr()
for i in range(n):
    cba[i] = x[i]


sotest.cal_vivj.restype=ctypes.c_float
sotest.mean_cij.restype=ctypes.c_float
print("py to c++ ")
# sotest.cal_vivj(cba,n1)
cr=sotest.mean_cij(cba,n1,n2)
print(cr)
