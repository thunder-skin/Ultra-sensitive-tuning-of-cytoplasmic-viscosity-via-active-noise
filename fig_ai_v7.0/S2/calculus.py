import numpy as np
from scipy.integrate import quad

# 定义积分的函数
def integrand(x, t):
    return np.exp(-x * t) * x**(-0.5)

# 定义计算定积分的函数
def calculate_integral(t):
    result, error = quad(integrand, 0.00000001, 1, args=(t,))
    return result


for i in np.arange(4,5,0.1):
    t=10**i
    integral_value = calculate_integral(t)
    print(t,integral_value)
    
    
