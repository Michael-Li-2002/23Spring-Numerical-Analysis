from math import pi
# 数值分析 第二次上机作业 2100010793 李佳
# 1. 使用的函数
# 2. 输出各方法误差除以误差阶后的数据(即为上机报告中图的数据来源 可直接用于LaTex画图)
# 单击运行即可

## 1. 使用的函数：
def f(x):
    return 4/(1+x**2)

def int_mid(f,n,a,b): # 复合中点
    I = 0
    points = [a+i*(b-a)/n for i in range(n+1)]
    for i in range(n):
        I += f((points[i]+points[i+1])/2) * (b-a)/n
    return I

def int_ti(f,n,a,b): # 复合梯形
    I = 0
    points = [a+i*(b-a)/n for i in range(n+1)]
    for i in range(n):
        I += (f(points[i]) + f(points[i+1]))/2 * (b-a)/n
    return I

def int_simpson(f,n,a,b): # 复合Simpson
    I = 0
    points = [a + i * (b - a) / n for i in range(n + 1)]
    for i in range(n):
        I += (f(points[i]) + f(points[i+1]) + 4*f((points[i] + points[i + 1]) / 2) ) /6 * (b - a) / n
    return I

def int_romberg(f,n,k,a,b): # Romberg求积方法
    I = [int_ti(f,2**i*n,a,b) for i in range(k)]
    for i in range(k-1):
        for j in range(k-i-1):
            I[j] = (4**(i+1)*I[j+1]-I[j])/(4**(i+1)-1)
    return I[0]

def adaptive(f,a,b,eps): # 自适应方法
    Sleft = a; Sright = a; Aleft = a; Aright = b #;  Nleft = a; Nright = b
    integral = 0
    while Sright != b:
        mid = (Aleft + Aright)/2
        S1,S2,S = int_simpson(f,1,Aleft,mid),int_simpson(f,1,mid,Aright),int_simpson(f,1,Aleft,Aright)
        err = 1/10 * (S1 + S2 - S)
        if abs(err) < eps * (Aright - Aleft)/(b-a):
            integral += S1 + S2
            Sright = Aright
            Aleft = Aright; Aright = b
        else:
            Aright = mid
    return integral

## 2. 输出误差除以误差阶后的数据(即为上机报告中图的数据来源 可直接用于LaTex画图)
a,b = 0, 1
print('复合中点数据：')
for j in range(5,1001,5):
    print('(',j,',',abs(int_mid(f,j,a,b)-pi) * j**2,')',end=' ')
print()

print('复合梯形数据：')
for j in range(5,1001,5):
    print('(',j,',',abs(int_ti(f,j,a,b)-pi) * j**2,')',end=' ')
print()

print('复合Simpson数据：')
for j in range(5,301,5):
    print('(',j,',',abs(int_simpson(f,j,a,b)-pi) * j**4,')',end=' ')
print()

print('Romberg(k=2)数据：')
for j in range(5,301,5):
    print('(',j,',',abs(int_romberg(f,j,2,a,b)-pi) * j**4,')',end=' ')
print()

print('Romberg(k=3)数据：')
for j in range(5,301,5):
    print('(',j,',',abs(int_romberg(f,j,3,a,b)-pi) * j**6,')',end=' ')
print()

print('Romberg(k=4)数据：')
for j in range(2,91,2):
    print('(',j,',',abs(int_romberg(f,j,4,a,b)-pi) * j**8,')',end=' ')
print()

print('自适应求解结果')
for i in range(12):
    print('(',10**(-i-1),',',abs(adaptive(f,a=0,b=1,eps=10**(-i-1))-pi),')',end=' ')



