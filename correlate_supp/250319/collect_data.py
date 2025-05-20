import statistics

file_read=8
name="3.5_1000"

time=[0.02,0.03,0.05,0.07,
      0.1,0.2,0.3,0.5,0.7,
      1,2,3,5,7,
      10,15,20,30,50,70,
      1e2,1.5e2,2e2,3e2,5e2,7e2,
      1e3,1.5e3,2e3,3e3,5e3,7e3,
      1e4,1.5e4,2e4,3e4,4e4,5e4,6e4,7e4,8e4,9e4,
      1e5,1.5e5,2e5,3e5,4e5,5e5,6e5,7e5,8e5,9e5,
      1e6,1.5e6,2e6,3e6,4e6,5e6,6e6,7e6,8e6,9e6,
      1e7,1.5e7,2e7,3e7,4e7,5e7,6e7,7e7,8e7,9e7]

def read(file_path):
    luni=[]
    visc=[]
    ener=[]
    take=0
    with open(file_path, 'r') as file:
        for line in file:
            if take==0:
                take=1
                
            elif take==1:
                take=2
                
            else:
                # 去掉首尾空格并按空格分割
                parts = line.strip().split()
                # 取得第一个数据并添加到列表
                if parts:  # 确保行不为空
                    luni.append(float(parts[1]))
                    visc.append(float(parts[2]))
    return luni,visc

lunii=[]
viscc=[]

for j in range(file_read):
    file_path="./"+str(name)+"/"+str(1000+j)+"_0_luni.txt"
    a,b=read(file_path)
    lunii.append(a)
    viscc.append(b)

luniii=[[] for i in range(100)]
visccc=[[] for i in range(100)]

luni_fin=[]
visc_fin=[]

for d in lunii:
    for j in range(len(d)):
        luniii[j].append(d[j])
        
for d in viscc:
    for j in range(len(d)):
        visccc[j].append(d[j])
        
for i in luniii:
    if i:
        luni_fin.append(sum(i)/len(i))
        
for i in visccc:
    if i:
        visc_fin.append(sum(i)/len(i))

for i in range(len(luni_fin)):
    print(((luni_fin[i]*1e6)//1)/1e6)
    
    

