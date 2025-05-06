import numpy as np
import matplotlib.pyplot as plt

# 设置随机变量A的采样数量
num_samples = 10000000

# 生成随机变量A的样本，A在0.1到1之间均匀分布
A_samples = np.random.uniform(0, 1, num_samples)

# 计算B = exp(A)
B_samples = np.exp(2.3*(A_samples-1))
C_samples = 0.1+0.9*A_samples
D_samples = np.log(1.105+(2.718-1.105)*A_samples)

hist, bin_edges = np.histogram(B_samples, bins=500, density=True)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
plt.plot(bin_centers, hist, linestyle='-', color='blue', linewidth=1.5,label='Type 1')

hist, bin_edges = np.histogram(C_samples, bins=500, density=True)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
plt.plot(bin_centers, hist, linestyle='-', color='green', linewidth=1.5,label='Type 2')

hist, bin_edges = np.histogram(D_samples, bins=500, density=True)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
plt.plot(bin_centers, hist, linestyle='-', color='red', linewidth=1.5,label='Type 3')


plt.title('Distribution of Radius',fontsize=16)
plt.xlabel('Radius',fontsize=16)
plt.ylabel('Density',fontsize=16)
plt.legend(fontsize=16)
plt.xticks(np.arange(0, 1.2, 0.1),fontsize=12)
plt.yticks(np.arange(0, 5, 0.5),fontsize=12)
plt.grid(True)
plt.show()
print(min(B_samples))