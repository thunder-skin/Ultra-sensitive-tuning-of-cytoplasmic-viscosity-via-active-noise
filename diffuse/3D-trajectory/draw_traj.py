import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize


x = [0]
y = [0]
# 读取数据文件
for i in range(10):
    data_file = str(i)+'.txt'
    with open(data_file, 'r') as f:
        lines = f.readlines()
        
    #重归化    
    xp=x[-1]
    yp=y[-1]
    print(xp,yp)
    # 提取 x 和 y 坐标
    for line in lines:
        parts = line.strip().split()
        if len(parts) >= 2:
            x.append(float(parts[0])+xp)
            y.append(float(parts[1])+yp)
    
    
    
        
#数据居中
cx=min(x)/2+max(x)/2
cy=min(y)/2+max(y)/2
for i in range(len(x)):
    x[i]-=cx
    y[i]-=cy

# 创建颜色渐变
segments = []
colors = []

for i in range(len(x) - 1):
    segments.append([(x[i], y[i]), (x[i+1], y[i+1])])
    colors.append(i)

# 绘制轨迹连线
plt.figure(figsize=(8, 6))
plt.ylim(-100, 100)
plt.xlim(-100, 100)

# 调整颜色映射范围，从红色开始到紫色结束
lc = LineCollection(segments, cmap='plasma', norm=Normalize(vmin=0, vmax=len(segments)))

lc.set_array(colors)
plt.gca().add_collection(lc)


plt.xlabel('X Coordinate')
plt.ylabel('Y Coordinate')
plt.title('Particle Trajectory with Gradient Line')

plt.colorbar(lc, label='Time Step')
plt.show()
