import statistics

data_re=1
file_re=16
name="3.4"

def read(file_path):
    luni=[]
    visc=[]
    ener=[]
    take=0
    with open(file_path, 'r') as file:
        for line in file:
            if take==0:
                take=1
            else:
                # 去掉首尾空格并按空格分割
                parts = line.strip().split()
                # 取得第一个数据并添加到列表
                if parts:  # 确保行不为空
                    luni.append(float(parts[1]))
                    visc.append(float(parts[2]))
                    ener.append(float(parts[4]))
    return luni,visc,ener

lunii=[]
viscc=[]
enerr=[]
for i in range(data_re):
    for j in range(file_re):
        file_path="./"+str(name)+"/"+str(7000+j)+"_0_luni.txt"
        a,b,c=read(file_path)
        lunii.append(a)
        viscc.append(b)
        enerr.append(c)

length_l=[]
for i in lunii:
    length_l.append(len(i))
length_l.sort()
numm_l=length_l[10]

luniii=[]
visccc=[]
enerrr=[]
vis_vr=[]

for i in range(numm_l):
    a=[]
    b=[]
    c=[]
    for j in lunii:
        if len(j)>i:
            a.append(j[i])
    for j in viscc:
        if len(j)>i:
            b.append(j[i])
    for j in enerr:
        if len(j)>i:
            c.append(j[i])
            
    ave_l=sum(a)/len(a)
    ave_v=sum(b)/len(b)
    ave_e=sum(c)/len(c)
    var_v=sample_std_dev = statistics.stdev(b)
    luniii.append(ave_l)
    visccc.append(ave_v)
    enerrr.append(ave_e)
    vis_vr.append(var_v)

i_c=0
for i in range(len(luniii)):
    print((luniii[i]*1e6)//1/1e6)
    i_c+=1
    
    

