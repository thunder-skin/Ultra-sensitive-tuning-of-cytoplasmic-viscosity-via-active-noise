% 创建一个新的图形窗口
figure;

x_a=3;
x_b=7;

tor=1;
torr=0.5;


% 读取粒子数据
function particles = read_particles(file_path)
    fid = fopen(file_path, 'r');
    if fid == -1
        error('File not found');
    end
    particles = [];
    while ~feof(fid)
        line = fgetl(fid);
        data = str2num(line);
        if length(data) == 5
            particles = [particles; data];
        end
    end
    fclose(fid);
end

file_path = '1002.txt'; % 替换为你的txt文件路径
particles = read_particles(file_path);

% 设置空间范围
fig=figure
set(fig, 'Position', [50, 50, 800, 700]);
axis([x_a-tor x_b+tor x_a-tor x_b+tor x_a-tor x_b+tor]); % 设置坐标轴范围为10x10x10的空间
hold on; % 保持当前图形


c1=[1.0 0.3 0.3];
c2=[0.3 1.0 0.3];
c3=[0.3 0.3 1.0];


for i = 1:size(particles, 1)
    x = particles(i, 1);
    y = particles(i, 2);
    z = particles(i, 3);
    r = particles(i, 4);
    id = particles(i, 5);

    if x<x_b && y<x_b && z<x_b && x>x_a && y>x_a && z>x_a
        [X1, Y1, Z1] = sphere(200); % 生成球体的坐标
        P=(r-0.1)/0.9
        if P<0.5
            N=2*P;
            c4=N*c2+(1-N)*c1;
        end
        if P>0.5
            N=2*P-1;
            c4=N*c3+(1-N)*c2;
        end
        surf(r*X1+x, r*Y1+y, r*Z1+z, 'FaceColor',c4, 'EdgeColor', 'none'); % 蓝色
    else
        if x<x_b+torr && y<x_b+torr && z<x_b+torr && x>x_a-torr && y>x_a-torr && z>x_a-torr && r<0.5
            [X1, Y1, Z1] = sphere(200); % 生成球体的坐标
            P=(r-0.1)/0.9
            if P<0.5
                N=2*P;
                c4=N*c2+(1-N)*c1;
            end
            if P>0.5
                N=2*P-1;
                c4=N*c3+(1-N)*c2;
            end
            surf(r*X1+x, r*Y1+y, r*Z1+z, 'FaceColor',c4, 'EdgeColor', 'none'); % 蓝色
        end
    end
end

line([x_a x_b], [x_a x_a], [x_a,x_a], 'Color', [0 0 0], 'LineStyle', '-.', 'LineWidth', 2);
line([x_a x_b], [x_b x_b], [x_a,x_a], 'Color', [0 0 0], 'LineStyle', '-.', 'LineWidth', 2);
line([x_a x_b], [x_b x_b], [x_b,x_b], 'Color', [0 0 0], 'LineStyle', '-.', 'LineWidth', 2);
line([x_a x_b], [x_a x_a], [x_b,x_b], 'Color', [0 0 0], 'LineStyle', '-.', 'LineWidth', 2);
line([x_a x_a], [x_b x_a], [x_a,x_a], 'Color', [0 0 0], 'LineStyle', '-.', 'LineWidth', 2);
line([x_a x_a], [x_b x_a], [x_b,x_b], 'Color', [0 0 0], 'LineStyle', '-.', 'LineWidth', 2);
line([x_b x_b], [x_b x_a], [x_b,x_b], 'Color', [0 0 0], 'LineStyle', '-.', 'LineWidth', 2);
line([x_b x_b], [x_b x_a], [x_a,x_a], 'Color', [0 0 0], 'LineStyle', '-.', 'LineWidth', 2);
line([x_a x_a], [x_a x_a], [x_a,x_b], 'Color', [0 0 0], 'LineStyle', '-.', 'LineWidth', 2);
line([x_b x_b], [x_a x_a], [x_a,x_b], 'Color', [0 0 0], 'LineStyle', '-.', 'LineWidth', 2);
line([x_b x_b], [x_b x_b], [x_a,x_b], 'Color', [0 0 0], 'LineStyle', '-.', 'LineWidth', 2);
line([x_a x_a], [x_b x_b], [x_a,x_b], 'Color', [0 0 0], 'LineStyle', '-.', 'LineWidth', 2);


% 添加光照
camlight; % 添加摄像机光源
material shiny; % 设置材质为光亮

% 添加图例和坐标轴标签

grid on; % 显示网格

xticks(0:2:10);
yticks(0:2:10);
zticks(0:2:10);


%axis equal; % 设置坐标轴等比例
view(3); % 以三维视角查看
view(135,15);
hold off; % 释放图形
