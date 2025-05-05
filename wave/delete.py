def read():
    data_collect=[]
    filename="1000_0_0_luni.txt"
    file_path = filename  # 替换为你的文件路径
    with open(filename, 'r') as file:
            lines = file.readlines()
            count=-1
            for line in lines:
                if count==-1:
                    data_collect.append(line)
                    count+=1
                # 如果不包含要删除的字符串，则保留该行
                else:
                    line_data = line.rstrip('\n')
                    data=line.split(" ")
                    if int(data[0])==count:
                        count+=1
                        data_collect.append(line)
    with open(filename, 'w') as file:
        file.writelines(data_collect)
        
read()
                        
                    
