import statistics

def read(id):
    collect_1=[0 for i in range(50)]
    collect_2=[0 for i in range(50)]
    filename=str(id)+"_0_luni.txt"
    file_path = filename  # 替换为你的文件路径
    with open(file_path, 'r') as file:
        n=0
        count=0
        for line in file:
            if n==0:
                n=1
            else:
                line = line.rstrip('\n')
                data=line.split(" ")
                collect_1[count]=float(data[1])
                collect_2[count]=float(data[3])
                count+=1
    return collect_1,collect_2

data_all_1=[]
data_all_2=[]

for i in range(3):
    for j in range(8):
        a=1000+i*10+j
        q,r=read(a)
        data_all_1.append(q)
        data_all_2.append(r)

force=[]
force_dev=[]
visco=[]
visco_dev=[]

for i in range(50):
    means=[]
    for j in range(len(data_all_1)):
        if data_all_1[j][i] != 0:
            means.append(data_all_1[j][i])
    if means:
        value=statistics.mean(means)
        devia=statistics.stdev(means)/4.9
    else:
        value=0
        devia=0
    
    force.append(value)
    force_dev.append(devia)

for i in range(50):
    means=[]
    for j in range(len(data_all_2)):
        if data_all_2[j][i] != 0:
            means.append(data_all_2[j][i])
    if means:
        value=statistics.mean(means)
        devia=statistics.stdev(means)/4.9
    else:
        value=0
        devia=0
    
    visco.append(value)
    visco_dev.append(devia)

for i in range(50):
    print(force[i],force_dev[i],visco[i],visco_dev[i])
