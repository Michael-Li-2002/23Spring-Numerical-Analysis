from math import sqrt,sin,cos,exp
import numpy as np
# 数值分析 第三次上机作业 2100010793 李佳
# 1. 前面是使用的函数
# 2. 运行后输出各例子在数值报告中所使用的数据
# 例3在第171行可手动更改n的大小, 之后单击运行即可

### 使用的函数
def Powell(m,x):
    fi = [0 for _ in range(4*m)]
    for i in range(m):
        fi[4*i] = x[4*i] + 10 * x[4*i+1]
        fi[4*i+1] = sqrt(5) * (x[4*i+2]-x[4*i+3])
        fi[4*i+2] = (x[4*i+1]-2*x[4*i+2])**2
        fi[4*i+3] = sqrt(10) * (x[4*i]-x[4*i+3])**2
    return fi

def tri(n,x):
    fi = [n for _ in range(n)]
    s = 0
    for j in range(n):
        s += cos(x[j])
    for i in range(n):
        fi[i] -= s - n*sin(x[i]) + n*(i+1)*(1-cos(x[i]))
    return fi

def wood(x):
    fi = [0 for _ in range(4)]
    fi[0] = 400*x[0]*(x[0]**2-x[1]) + 2*(x[0]-1)
    fi[1] = 200*(x[1]-x[0]**2) + 20.2*(x[1]-1) + 19.8*(x[3]-1)
    fi[2] = 360*x[2]*(x[2]**2-x[3]) + 2*(x[2]-1)
    fi[3] = 180*(x[3]-x[2]**2) + 20.2*(x[3]-1) + 19.8*(x[1]-1)
    return fi

def myinner(x,y):
    n = len(x)
    sum = 0
    for i in range(n):
        sum += x[i]*y[i]
    return sum

def myprod(A,x):
    n = len(x)
    ans = [0 for _ in range(n)]
    for i in range(n):
        for j in range(n):
            ans[i] += A[i][j] * x[j]
    return ans

def mydot(x,y):
    n = len(x)
    ans = [[x[i]*y[j] for j in range(n)] for i in range(n)]
    return ans

def mytrans(A):
    n = len(A)
    ans = [[A[j][i] for j in range(n)]for i in range(n)]
    return ans

def myminus(x,y):
    ans = [x[i]-y[i] for i in range(len(x))]
    return ans

def myminus2(A,B):
    ans = [myminus(A[i],B[i]) for i in range(len(A))]
    return ans

def mymult1(c,x):
    ans = [c*x[i] for i in range(len(x))]
    return ans

def mymult2(c,A):
    ans = [mymult1(c,A[i]) for i in range(len(A))]
    return ans

def Broyden(f,x0,invA0):
    # x1 = x0 - np.dot(invA0,f(x0))
    # k = 0
    # y0 = x1 - x0
    # while np.dot(y0,y0.T)[0,0] > 10**(-8):
    #     k += 1
    #     y1 = f(x1)
    #     invA0 -= 1/(np.dot(y0,y0.T)[0,0] + np.dot(y0,np.dot(invA0,y1.T))[0,0]) * np.dot(np.dot(invA0.T,y0).T,np.dot(invA0,y1))
    #     x2 = x1 - inv
    x1 = myminus(x0,myprod(invA0, f(x0)))
    k = 0
    y0 = myminus(x1,x0)
    while myinner(y0,y0) > 10 ** (-8):
        k += 1
        y1 = f(x1)
        invA0 = myminus2( invA0 , mymult2( 1 / (myinner(y0,y0) + myinner(y0, myprod(invA0, y1))), mydot(myprod(mytrans(invA0), y0),myprod(invA0, y1)) ) )
        x2 = myminus( x1, myprod(invA0,f(x1)) )
        x0 = x1; x1 = x2
    return x1

### 2.Powell Function
  ## Newton Method  只需计算 m=1, n=4的情形
m = 1
def f2(x):
    return Powell(m, x)
def JacobPowell(x):
    A = [[0 for _ in range(4)]for _ in range(4)]
    A[0][0] = 1; A[0][1]=10; A[1][2] = sqrt(5); A[1][3] = -sqrt(5)
    A[2][1] = 2 * (x[1] - 2 * x[2])
    A[2][2] = 4 * (2 * x[2] - x[1])
    A[3][0] = 2 * sqrt(10) * (x[0] - x[3])
    A[3][3] = 2 * sqrt(10) * (x[3] - x[0])
    return A
x0 = [3,-1,0,1]
k = 0
result = []
print('2.Powell Function: Newton Method')
print('(k:第k步迭代, 误差):')
#print(k,x1)
while myinner(x0,x0) > 10 ** (-16):#myinner(y0,y0) > 10 ** (-16):
    k += 1
    A0 = JacobPowell(x0)
    inverse = np.linalg.inv(np.mat(A0))
    invA0 = [[inverse[i, j] for j in range(4)] for i in range(4)]
    x0 = myminus( x0, myprod(invA0,Powell(1,x0)) )
    #print('第',k,'步, 误差:',sqrt(myinner(x0,x0)))
    result.append((k,sqrt(myinner(x0, x0))))
    #print(k,x1)
for c in result:
    print(c,end=' ')
print()
print('迭代次数:',k)
print('收敛的点:',x0)
print('函数值2-范数:',sqrt(myinner(f2(x0),f2(x0))))
print()

  ## Broyden Method
x0 = [0 for _ in range(4*m)]
for i in range(m):
    x0[4*i] = float(3.00)
    x0[4*i+1] = float(-1.00)
    x0[4*i+3] = float(1.00)
invA0 = [[0 for _ in range(4*m)]for _ in range(4*m)]
for i in range(m):
    invA0[4*i][4*i:4*i+4] = [ 0.04761905,-0.42591771, 0.23809524, 0.07529233]
    invA0[4 * i+1][4 * i:4 * i + 4] =[ 0.0952381,  0.04259177,-0.02380952,-0.00752923]
    invA0[4 * i+2][4 * i:4 * i + 4] =[ 0.04761905, 0.02129589, 0.23809524,-0.00376462]
    invA0[4 * i+3][4 * i:4 * i + 4] =[ 0.04761905,-0.42591771, 0.23809524,-0.00376462]
x1 = myminus(x0,myprod(invA0, f2(x0)))
k = 1
y0 = myminus(x1,x0)
result = [(k,sqrt(myinner(x1, x1)))]
print('2.Powell Function: Broyden Method')
print('(k:第k步迭代, 误差):')
#print('第', k, '步, 误差:', sqrt(myinner(x1, x1)))
while myinner(x1,x1) > 10 ** (-16):#myinner(y0,y0) > 10 ** (-20):
    k += 1
    y1 = f2(x1)
    invA0 = myminus2( invA0 , mymult2( 1 / (myinner(y0,y0) + myinner(y0, myprod(invA0, y1))), mydot(myprod(invA0, y1),myprod(mytrans(invA0), y0)) ) )
    x2 = myminus( x1, myprod(invA0,y1) )
    x0 = x1; x1 = x2
    y0 = myminus(x1,x0)
    #print('第', k, '步, 误差:', sqrt(myinner(x1, x1)))
    result.append((k,sqrt(myinner(x1, x1))))
    #print(k,x1)
for c in result:
    print(c,end=' ')
print()
print('迭代次数:',k)
print('收敛的点:',x1)
print('函数值2-范数:',sqrt(myinner(f2(x1),f2(x1))))
print()

### 3. Trigonometric Function
  ## Newton Method:
n = 10  # 此例子可更改方程个数
def f3(x):
    return tri(n,x)
def JacobTri(x):
    n = len(x)
    A = [[sin(x[j]) for j in range(n)] for i in range(n)]
    for i in range(n):
        A[i][i] += n*cos(x[i]) - n*(i+1)*sin(x[i])
    return A
x0 = [1/n for _ in range(n)]
# A0 = [[sin(1/n) for _ in range(n)]for _ in range(n)]
# for i in range(n):
#     A0[i][i] += n*cos(1/n) -n*(i+1)*sin(1/n)
# inverse = np.linalg.inv(np.mat(A0))
# invA0 = [[inverse[i,j] for j in range(n)]for i in range(n)]
# x1 = myminus(x0,myprod(invA0, f3(x0)))
k = 0
y0 = x0
result = []
print('3.Trignometric Function: Newton Method(初值(1/n,...,1/n))')
print('(k:第k步迭代, 相邻两步步长):')
#print(k,x1)
while myinner(y0,y0) > 10 ** (-16):
    k += 1
    A0 = JacobTri(x0)
    inverse = np.linalg.inv(np.mat(A0))
    invA0 = [[inverse[i, j] for j in range(n)] for i in range(n)]
    x1 = myminus( x0, myprod(invA0,f3(x0)) )
    y0 = myminus(x1,x0)
    x0 = x1
    #print('第', k, '步, 误差:', sqrt(myinner(y0, y0)))
    result.append((k, sqrt(myinner(y0, y0))))
    #print(k,x1)
# print(x1)
# print(sqrt(myinner(f3(x1),f3(x1))))
for c in result:
    print(c,end=' ')
print()
print('迭代次数:',k)
print('收敛的点:',x1)
print('函数值2-范数:',sqrt(myinner(f3(x1),f3(x1))))
print()

## Broyden Method:
x0 = [1/n for _ in range(n)]
A0 = [[sin(1/n) for _ in range(n)]for _ in range(n)]
for i in range(n):
    A0[i][i] += n*cos(1/n) -n*(i+1)*sin(1/n)
inverse = np.linalg.inv(np.mat(A0))
invA0 = [[inverse[i,j] for j in range(n)]for i in range(n)]
x1 = myminus(x0,myprod(invA0, f3(x0)))
k = 1
y0 = myminus(x1,x0)
result = [(k,sqrt(myinner(y0, y0)))]
#print(k,x1)
print('3.Trignometric Function: Broyden Method(初值(1/n,...,1/n))')
print('(k:第k步迭代, 相邻两步步长):')
while myinner(y0,y0) > 10 ** (-16) and k<=39:
    k += 1
    y1 = f3(x1)
    invA0 = myminus2( invA0 , mymult2( 1 / (myinner(y0,y0) + myinner(y0, myprod(invA0, y1))), mydot(myprod(invA0, y1),myprod(mytrans(invA0), y0)) ) )
    x2 = myminus( x1, myprod(invA0,y1) )
    x0 = x1; x1 = x2
    y0 = myminus(x1,x0)
    #print('第', k, '步, 相邻两步步长:', sqrt(myinner(y0, y0)))
    result.append((k,sqrt(myinner(y0, y0))))
    #print(k,x1)
# print(x1)
# print(sqrt(myinner(f3(x1),f3(x1))))
for c in result:
    print(c,end=' ')
print()
# print('迭代次数:',k)
# print('收敛的点:',x1)
# print('函数值2-范数:',sqrt(myinner(f3(x1),f3(x1))))
print()

x0 = [1/n**3 for _ in range(n)]
A0 = [[sin(1/n) for _ in range(n)]for _ in range(n)]
for i in range(n):
    A0[i][i] += n*cos(1/n) -n*(i+1)*sin(1/n)
inverse = np.linalg.inv(np.mat(A0))
invA0 = [[inverse[i,j] for j in range(n)]for i in range(n)]
x1 = myminus(x0,myprod(invA0, f3(x0)))
k = 1
y0 = myminus(x1,x0)
result = [(k,sqrt(myinner(y0, y0)))]
#print(k,x1)
print('3.Trignometric Function: Broyden Method(初值(1/n^3,...,1/n^3))')
print('(k:第k步迭代, 相邻两步步长):')
while myinner(y0,y0) > 10 ** (-16):# and k<=39:
    k += 1
    y1 = f3(x1)
    invA0 = myminus2( invA0 , mymult2( 1 / (myinner(y0,y0) + myinner(y0, myprod(invA0, y1))), mydot(myprod(invA0, y1),myprod(mytrans(invA0), y0)) ) )
    x2 = myminus( x1, myprod(invA0,y1) )
    x0 = x1; x1 = x2
    y0 = myminus(x1,x0)
    #print('第', k, '步, 相邻两步步长:', sqrt(myinner(y0, y0)))
    result.append((k,sqrt(myinner(y0, y0))))
    #print(k,x1)
# print(x1)
# print(sqrt(myinner(f3(x1),f3(x1))))
for c in result:
    print(c,end=' ')
print()
print('迭代次数:',k)
print('收敛的点:',x1)
print('函数值2-范数:',sqrt(myinner(f3(x1),f3(x1))))
print()

### 5.Wood Function
  ## Newton Method:
def JacobWood(x):
    # fi[0] = 400 * x[0] * (x[0] ** 2 - x[1]) + 2 * (x[0] - 1)
    # fi[1] = 200 * (x[1] - x[0] ** 2) + 20.2 * (x[1] - 1) + 19.8 * (x[3] - 1)
    # fi[2] = 360 * x[2] * (x[2] ** 2 - x[3]) + 2 * (x[2] - 1)
    # fi[3] = 180 * (x[3] - x[2] ** 2) + 20.2 * (x[3] - 1) + 19.8 * (x[1] - 1)
    A = [[0 for j in range(4)] for i in range(4)]
    A[0][0] = 1200 * x[0]**2 - 400 * x[1] + 2
    A[0][1] = -400 * x[0]
    A[1][0] = -400*x[0]; A[1][1] = 220.2; A[1][3] = 19.8
    A[2][2] = 1080 * x[2]**2 - 360 * x[3] + 2
    A[2][3] = -360 * x[2]
    A[3][1] = 19.8; A[3][2] = -360 * x[2]; A[3][3] = 200.2
    return A
def f5(x):
    return wood(x)
x0 = [-3,-1,-3,-1]
A0 = JacobWood(x0)
inverse = np.linalg.inv(np.mat(A0))
invA0 = [[inverse[i,j] for j in range(4)]for i in range(4)]
x1 = myminus(x0,myprod(invA0, f5(x0)))
k = 1
y0 = myminus(x1,x0)
print('5.Wood Function: Newton Method')
print('(第k步, 优化函数值f(x_k)): ')
#print(k,100*(x1[0]**2-x1[1])**2+(1-x1[0])**2 + 90*(x1[2]**2-x1[3])**2+(1-x1[2])**2 + 10.1*(1-x1[1])**2+10.1*(1-x1[3])**2 + 19.8*(1-x1[1])*(1-x1[3]))
result = [(k,100*(x1[0]**2-x1[1])**2+(1-x1[0])**2 + 90*(x1[2]**2-x1[3])**2+(1-x1[2])**2 + 10.1*(1-x1[1])**2+10.1*(1-x1[3])**2 + 19.8*(1-x1[1])*(1-x1[3]))]
while myinner(y0,y0) > 10 ** (-16):
    k += 1
    A0 = np.mat(JacobWood(x1))
    #A0 += 0.1*np.eye(4)
    inverse = np.linalg.inv(A0)
    invA0 = [[inverse[i, j] for j in range(4)] for i in range(4)]
    x2 = myminus( x1, myprod(invA0,f5(x1)) )
    x0 = x1; x1 = x2
    y0 = myminus(x1,x0)
    result.append((k,100*(x1[0]**2-x1[1])**2+(1-x1[0])**2 + 90*(x1[2]**2-x1[3])**2+(1-x1[2])**2 + 10.1*(1-x1[1])**2+10.1*(1-x1[3])**2 + 19.8*(1-x1[1])*(1-x1[3])))
    #print(k,100*(x1[0]**2-x1[1])**2+(1-x1[0])**2 + 90*(x1[2]**2-x1[3])**2+(1-x1[2])**2 + 10.1*(1-x1[1])**2+10.1*(1-x1[3])**2 + 19.8*(1-x1[1])*(1-x1[3]))
# print(100*(x1[0]**2-x1[1])**2+(1-x1[0])**2 + 90*(x1[2]**2-x1[3])**2+(1-x1[2])**2 + 10.1*(1-x1[1])**2+10.1*(1-x1[3])**2 + 19.8*(1-x1[1])*(1-x1[3]))
# print(x1)
for c in result:
    print(c,end=' ')
print()
print('迭代次数:',k)
print('收敛的点:',x1)
print('优化函数值:',100*(x1[0]**2-x1[1])**2+(1-x1[0])**2 + 90*(x1[2]**2-x1[3])**2+(1-x1[2])**2 + 10.1*(1-x1[1])**2+10.1*(1-x1[3])**2 + 19.8*(1-x1[1])*(1-x1[3]))
print()

x0 = [-3,-1,-3,-1]
A0 = JacobWood(x0)
inverse = np.linalg.inv(np.mat(A0))
invA0 = [[inverse[i,j] for j in range(4)]for i in range(4)]
x1 = myminus(x0,myprod(invA0, f5(x0)))
k = 1
y0 = myminus(x1,x0)
print('5.Wood Function: 改进的Newton Method(\mu_k=0.1)')
print('(第k步, 优化函数值f(x_k)):  ')
# print(k,100*(x1[0]**2-x1[1])**2+(1-x1[0])**2 + 90*(x1[2]**2-x1[3])**2+(1-x1[2])**2 + 10.1*(1-x1[1])**2+10.1*(1-x1[3])**2 + 19.8*(1-x1[1])*(1-x1[3]))
result = [(k,100*(x1[0]**2-x1[1])**2+(1-x1[0])**2 + 90*(x1[2]**2-x1[3])**2+(1-x1[2])**2 + 10.1*(1-x1[1])**2+10.1*(1-x1[3])**2 + 19.8*(1-x1[1])*(1-x1[3]))]
while myinner(y0,y0) > 10 ** (-16):
    k += 1
    A0 = np.mat(JacobWood(x1))
    A0 += 0.1*np.eye(4)
    inverse = np.linalg.inv(A0)
    invA0 = [[inverse[i, j] for j in range(4)] for i in range(4)]
    x2 = myminus( x1, myprod(invA0,f5(x1)) )
    x0 = x1; x1 = x2
    y0 = myminus(x1,x0)
    result.append((k,100*(x1[0]**2-x1[1])**2+(1-x1[0])**2 + 90*(x1[2]**2-x1[3])**2+(1-x1[2])**2 + 10.1*(1-x1[1])**2+10.1*(1-x1[3])**2 + 19.8*(1-x1[1])*(1-x1[3])))
    #print(k,100*(x1[0]**2-x1[1])**2+(1-x1[0])**2 + 90*(x1[2]**2-x1[3])**2+(1-x1[2])**2 + 10.1*(1-x1[1])**2+10.1*(1-x1[3])**2 + 19.8*(1-x1[1])*(1-x1[3]))
# print(100*(x1[0]**2-x1[1])**2+(1-x1[0])**2 + 90*(x1[2]**2-x1[3])**2+(1-x1[2])**2 + 10.1*(1-x1[1])**2+10.1*(1-x1[3])**2 + 19.8*(1-x1[1])*(1-x1[3]))
# print(x1)
# print(invA0)
for c in result:
    print(c,end=' ')
print()
print('迭代次数:',k)
print('收敛的点:',x1)
print('优化函数值:',100*(x1[0]**2-x1[1])**2+(1-x1[0])**2 + 90*(x1[2]**2-x1[3])**2+(1-x1[2])**2 + 10.1*(1-x1[1])**2+10.1*(1-x1[3])**2 + 19.8*(1-x1[1])*(1-x1[3]))
print()

print('  改进Newton法取值各\mu_k时达到极小值点所需迭代次数:')
print('  (\mu_k,迭代次数):')
for eps in [0.01,0.02,0.05,0.06,0.07,0.08,0.09,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,2.0,3.0,4.0]:
    x0 = [-3, -1, -3, -1]
    A0 = JacobWood(x0)
    inverse = np.linalg.inv(np.mat(A0))
    invA0 = [[inverse[i, j] for j in range(4)] for i in range(4)]
    x1 = myminus(x0, myprod(invA0, f5(x0)))
    k = 1
    y0 = myminus(x1,x0)
    while myinner(y0, y0) > 10 ** (-16):
        k += 1
        A0 = np.mat(JacobWood(x1))
        A0 += eps * np.eye(4)
        inverse = np.linalg.inv(A0)
        invA0 = [[inverse[i, j] for j in range(4)] for i in range(4)]
        x2 = myminus(x1, myprod(invA0, f5(x1)))
        x0 = x1; x1 = x2
        y0 = myminus(x1, x0)
    print((eps,k),end=' ')
print()
print()

## Broyden Method 1
x0 = [-3,-1,-3,-1]
A0 = JacobWood(x0)
inverse = np.linalg.inv(np.mat(A0))
invA0 = [[inverse[i,j] for j in range(4)]for i in range(4)]
x1 = myminus(x0,myprod(invA0, f5(x0)))
k = 1
y0 = myminus(x1,x0)
print('5.Wood Function: Broyden Method 1')
print('(第k步, 优化函数值f(x_k)):  (数量太多, 仅取部分展示)')
#print(k,x1)
while myinner(y0,y0) > 10 ** (-16):
    k += 1
    y1 = f5(x1)
    invA0 = myminus2( invA0 , mymult2( 1 / (myinner(y0,y0) + myinner(y0, myprod(invA0, y1))), mydot(myprod(invA0, y1),myprod(mytrans(invA0), y0)) ) )
    # for i in range(4):
    #     invA0[i][i] += 0.02
    x2 = myminus( x1, myprod(invA0,y1) )
    x0 = x1; x1 = x2
    y0 = myminus(x1,x0)
    if k % 50 == 2:
        print('(',k,',',100 * (x1[0] ** 2 - x1[1]) ** 2 + (1 - x1[0]) ** 2 + 90 * (x1[2] ** 2 - x1[3]) ** 2 + (
                    1 - x1[2]) ** 2 + 10.1 * (1 - x1[1]) ** 2 + 10.1 * (1 - x1[3]) ** 2 + 19.8 * (1 - x1[1]) * (
                          1 - x1[3]),')',end='')
    #print(k,x1)
print()
print('迭代次数:',k)
print('收敛的点:',x1)
print('优化函数值:',100*(x1[0]**2-x1[1])**2+(1-x1[0])**2 + 90*(x1[2]**2-x1[3])**2+(1-x1[2])**2 + 10.1*(1-x1[1])**2+10.1*(1-x1[3])**2 + 19.8*(1-x1[1])*(1-x1[3]))
print()
## Broyden Method 2
x0 = [-3,-1,-3,-1]
A0 = JacobWood(x0)
inverse = np.linalg.inv(np.mat(A0))
invA0 = [[inverse[i,j] for j in range(4)]for i in range(4)]
x1 = myminus(x0,myprod(invA0, f5(x0)))
k = 1
y0 = myminus(x1,x0)
print('5.Wood Function: Broyden Method 2')
print('(第k步, 优化函数值f(x_k)):')
print('(', k, ',', 100 * (x1[0] ** 2 - x1[1]) ** 2 + (1 - x1[0]) ** 2 + 90 * (x1[2] ** 2 - x1[3]) ** 2 + (
                             1 - x1[2]) ** 2 + 10.1 * (1 - x1[1]) ** 2 + 10.1 * (1 - x1[3]) ** 2 + 19.8 * (1 - x1[1]) * (
                                   1 - x1[3]),')',end='')
#print(k,x1)
while myinner(y0,y0) > 10 ** (-16):
    k += 1
    y1 = f5(x1)
    invA0 = myminus2( invA0 , mymult2( 1 / (myinner(y1,y0) + myinner(y1, myprod(invA0, y1))), mydot(myprod(invA0, y1), myprod(invA0, y1)) ) )
    x2 = myminus( x1, myprod(invA0,y1) )
    x0 = x1; x1 = x2
    y0 = myminus(x1,x0)
    print('(', k, ',', 100 * (x1[0] ** 2 - x1[1]) ** 2 + (1 - x1[0]) ** 2 + 90 * (x1[2] ** 2 - x1[3]) ** 2 + (
                             1 - x1[2]) ** 2 + 10.1 * (1 - x1[1]) ** 2 + 10.1 * (1 - x1[3]) ** 2 + 19.8 * (1 - x1[1]) * (
                                   1 - x1[3]),')',end='')
    #print(k,x1)
print()
print('迭代次数:',k)
print('收敛的点:',x1)
print('优化函数值:',100*(x1[0]**2-x1[1])**2+(1-x1[0])**2 + 90*(x1[2]**2-x1[3])**2+(1-x1[2])**2 + 10.1*(1-x1[1])**2+10.1*(1-x1[3])**2 + 19.8*(1-x1[1])*(1-x1[3]))
