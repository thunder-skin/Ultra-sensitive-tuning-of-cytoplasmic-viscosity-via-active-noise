import statistics
import matplotlib.pyplot as plt

def read(id):
    data_collect=[0 for i in range(1000)]
    filename=str(1000+id)+"_0_0_luni.txt"
    file_path = filename  # 替换为你的文件路径
    count=-1
    with open(file_path, 'r') as file:
        for line in file:
            if count==-1:
                count+=1
            else:
                line = line.rstrip('\n')
                data=line.split(" ")
                data_collect[count]=float(data[1])
                count+=1
    return data_collect

data_all=[]
for j in range(4):
    for i in range(8):
        data_all.append(read(10*j+i))

mean=[]
dev=[]
for i in range(1000):
    test=[]
    for j in range(len(data_all)):
        if data_all[j][i]:
            test.append(data_all[j][i])
    average = sum(test) / len(test)
    std_dev = statistics.stdev(test)/5.7
    mean.append(average)
    dev.append(std_dev)
    
mean_=[]
for i in range(len(mean)-1):
    mean_.append(mean[i]/2+mean[i+1]/2)
mean=mean_[0:500]
dev=dev[0:500]


A=max(mean)
B=min(mean)
P=mean.index(A)
Q=mean.index(B)
R=dev[P]+dev[Q]

print(A-B,R)
f=(A+B)/2
g=(A-B)/2
for i in range(len(mean)):
    mean[i]=(mean[i]-f)/g
            
x=[i+0.5 for i in range(len(mean))]
plt.plot(x,mean)
plt.title('Output Curve with Period T=1e5')
#plt.errorbar(x, mean, yerr=dev, fmt='o', capsize=5, label='Data with Error Bars')
plt.show()

def write_file( a,b):
    with open("1e5.txt", 'a') as file:
        file.write(f"{a}, {b}\n")
            
for i in range(len(mean)):
    write_file(x[i],mean[i])
            
            
            

