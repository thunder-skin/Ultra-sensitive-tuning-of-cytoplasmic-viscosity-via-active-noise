% 创建一个新的图形窗口
figure;

x_a=10;
x_b=25;

tor=2;
torr=1;


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

file_path = '1000.txt'; % 替换为你的txt文件路径
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
        if r<1.2
            c4=c1;
        end
        if r>1.2
            c4=c2;
        end
        r=r*0.9;
        surf(r*X1+x, r*Y1+y, r*Z1+z, 'FaceColor',c4, 'EdgeColor', 'none'); % 蓝色
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

xticks(10:5:25);
yticks(10:5:25);
zticks(10:5:25);


%axis equal; % 设置坐标轴等比例
view(3); % 以三维视角查看
view(135,15);
hold off; % 释放图形
