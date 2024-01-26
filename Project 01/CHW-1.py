from math import *
import matplotlib.pyplot as plt
# 2100010793 李佳
# 前面部分为使用的函数
# 后面的1.给出一个三次样条插值的数值算例（图1）
# 2.给出收敛阶的检验（n^4*e_h ~ O(1)）（图2）
# 3.检验收敛阶中的数据输出 用于LaTeX绘图


# 插值函数
def f(x):
    return e**(sin(x)) + cos(4*x)


# 三次样条插值 f:插值函数 points:插值节点 x:要计算的插值多项式上的点 输出x上对应插值多项式函数值y
def cubicspline(f,points,x):
    yi = [f(c) for c in points]
    mi = thissolver(f,points)  # 求解导数值
    index = 0
    y = [0 for _ in range(len(x))]
    for i in range(len(x)):
        while x[i] > points[index+1]:
            index += 1
            if index >= len(points) - 1:
                index = len(points) - 2
                break
        # 按照插值基函数求插值多项式上的值
        y[i] = yi[index]*(1+2*(x[i]-points[index])/(points[index+1]-points[index]))*((x[i]-points[index+1])\
               /(points[index]-points[index+1]))**2 + yi[index+1]*(1+2*(x[i]-points[index+1])/(points[index]-points[index+1]))*((x[i]-points[index])\
               /(points[index+1]-points[index]))**2\
               + mi[index]*(x[i]-points[index])*((x[i]-points[index+1])/(points[index]-points[index+1]))**2\
               + mi[index+1]*(x[i]-points[index+1])*((x[i]-points[index])/(points[index+1]-points[index]))**2
    return y


# 追赶法求解并输出插值多项式各点导数值 func: 插值函数 points: 插值节点
def thissolver(func,points):
    n = len(points) - 1
    yi = [func(points[i]) for i in range(n+1)]
    mat = [[0 for _ in range(n)] for _ in range(n)]
    lam = [0 for _ in range(n)]; mu = [0 for _ in range(n)]
    lam[1:] = [(points[i]-points[i-1])/(points[i+1]-points[i-1]) for i in range(1,n)]
    lam[0] = (points[n]-points[n-1])/(points[n]-points[n-1]+points[1]-points[0])
    mu[1:] = [3*((1-lam[i])/(points[i]-points[i-1])*(yi[i]-yi[i-1]) + (lam[i])/(points[i+1]-points[i])*(yi[i+1]-yi[i])) for i in range(1,n)]
    mu[0] = 3*((1-lam[0])/(points[n]-points[n-1])*(yi[n]-yi[n-1]) + (lam[0])/(points[1]-points[0])*(yi[1]-yi[0]))
    mat[0][0]=4; mat[0][1]=1; mat[0][n-1] = 1
    mat[n-1][n-1] = 2; mat[n-1][n-2] = 1-lam[n-1]; mat[n-1][0] = lam[n-1]
    for i in range(0,n-1):
        mat[i][i] = 2; mat[i][i+1] = lam[i]; mat[i][i-1] = 1-lam[i]
    for i in range(n-2):
        mat[i+1][i] /= mat[i][i]
        mat[i+1][i+1] -= mat[i+1][i] * mat[i][i+1]
        mat[i+1][-1] -= mat[i+1][i] * mat[i][-1]
        mat[-1][i] /= mat[i][i]
        mat[-1][i + 1] -= mat[-1][i] * mat[i][i + 1]
        mat[-1][-1] -= mat[-1][i] * mat[i][-1]
    mat[n-1][n-2] /= mat[n-2][n-2]
    mat[-1][-1] -= mat[n-1][n-2] * mat[n-2][n-1]

    for i in range(n-2):
        mu[i+1] -= mat[i+1][i] * mu[i]
        mu[-1] -= mat[-1][i] * mu[i]
    mu[n-1] -= mat[n-1][n-2] * mu[n-2]

    mu[n-1] /= mat[-1][-1]
    for j in range(n-1):
        mu[j] -= mat[j][-1] * mu[n-1]
    for i in range(n-2,0,-1):
        mu[i] /= mat[i][i]
        mu[i-1] -= mat[i-1][i] * mu[i]
    mu[0] /= mat[0][0]

    mu.append(mu[0])
    return mu


# 1. 三次样条插值的数值样例 对应画出的第一张图
n = 15
grid = 100
points = [2*pi*k/n for k in range(n+1)]
xaxis = [2*pi*i/grid for i in range(grid+1)]
yi = [f(x) for x in points]
realy = [f(x) for x in xaxis]
cubicy = cubicspline(f,points,xaxis)

plt.rcParams['font.sans-serif'] = ['STSong']
plt.plot(xaxis,realy,color='orange',label='函数')
plt.plot(xaxis,cubicy,color= 'red', label='周期三次样条插值')
plt.legend(loc="upper right")
plt.show()

# 2. 检验收敛阶 对应生成的第二张图
grid = 100
err = []
xaxis = [2*pi*i/grid for i in range(grid+1)]
realy = [f(x) for x in xaxis]
for n in range(4,91,2):
    points = [2*pi*k/n for k in range(n+1)]
    yi = [f(x) for x in points]
    cubicy = cubicspline(f,points,xaxis)
    error = [abs(cubicy[i]-realy[i]) for i in range(len(xaxis))]
    err.append(max(error)*n**4)
plt.plot(list(range(4,91,2)),err,label='n^4*e_h')
plt.legend(loc="upper right")
plt.show()

# 3.“检验收敛阶”中得到的数据 用于在LaTeX中绘制折线图
for n in range(1,45):
    print('(',2*n+2,',',err[n-1],')',end=' ')
