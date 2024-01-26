# 数值分析 第四次上机作业 2100010793 李佳
# 1. 使用的函数: 4级4阶Runge-Kutta古典显式公式
# 2. 主程序: 第30行开始设置参数: 输入时长、步长、步数(时长/步长)及相应初值,
#           输出数值解在各时刻的取值、x分量在各时刻的取值(均用于上机报告内的作图)

## 1. 函数
def RK4(f,gap,y):  # 自治方程的4级4阶Runge-Kutta古典显式
    n = len(y)
    K1 = f(y)
    y1 = [y[i] for i in range(n)]
    for i in range(n):
        y1[i] += gap/2*K1[i]
    K2 = f(y1)
    y2 = [y[i] for i in range(n)]
    for i in range(n):
        y2[i] += gap / 2 * K2[i]
    K3 = f(y2)
    y3 = [y[i] for i in range(n)]
    for i in range(n):
        y3[i] += gap * K3[i]
    K4 = f(y3)
    res = [y[i] for i in range(n)]
    for i in range(n):
        res[i] += gap/6*(K1[i]+2*K2[i]+2*K3[i]+K4[i])
    return res

## 2.主程序
# 参数设置
t = 20          # 时长
gap = 0.02      # 步长
num = 1000      # 步数：num = t//gap
p = [10,28,8/3] # 初值: [\sigma,\rho,\beta]

def f(y):  # 设置方程右端函数
    res = [p[0]*(y[1]-y[0]),p[1]*y[0]-y[1]-y[0]*y[2],y[0]*y[1]-p[2]*y[2]]
    return res

x = [0 for _ in range(num)]
y = [0 for _ in range(num)]
z = [0 for _ in range(num)]
(x[0],y[0],z[0]) = (0.01,0,0)
for i in range(1,num):
    (x[i],y[i],z[i]) = RK4(f,gap,[x[i-1],y[i-1],z[i-1]])

print('0-'+str(t)+'时间内, 以初值[sigma,rho,beta] =',p,'的数值解(轨迹上的点的坐标)为:')
print('( 输出格式: (x(t),y(t),z(t)) ... )')
for i in range(num):
    print('(',x[i],',',y[i],',',z[i],')',end=' ')
print()
print()

print('0-'+str(t)+'时间内, 以初值[sigma,rho,beta] =',p,'的数值解的x分量随时间的变化为')
print('( 输出格式: (t,x(t)) ... )')
for i in range(num):
    print('(',gap*i,',',x[i],')',end=' ')
print()

