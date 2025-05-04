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

% 筛选出在指定范围内的粒子
filtered_particles = particles( ...
    particles(:,1) > x_a-tor & particles(:,1) < x_b+tor & ...
    particles(:,2) > x_a-tor & particles(:,2) < x_b+tor & ...
    particles(:,3) > x_a-tor & particles(:,3) < x_b+tor, :);

% 设置空间范围
fig=figure;
set(fig, 'Position', [50, 50, 800, 700]);
axis([x_a-tor x_b+tor x_a-tor x_b+tor x_a-tor x_b+tor]); % 设置坐标轴范围
hold on; % 保持当前图形

c1=[1.0 0.3 0.3];
c2=[0.3 1.0 0.3];
c3=[0.3 0.3 1.0];

% 绘制接触的球之间的连线
for i = 1:size(filtered_particles, 1)
    for j = i+1:size(filtered_particles, 1)
        x1 = filtered_particles(i, 1);
        y1 = filtered_particles(i, 2);
        z1 = filtered_particles(i, 3);
        r1 = filtered_particles(i, 4);
        
        x2 = filtered_particles(j, 1);
        y2 = filtered_particles(j, 2);
        z2 = filtered_particles(j, 3);
        r2 = filtered_particles(j, 4);
        
        % 计算球心距离
        distance = sqrt((x2 - x1)^2 + (y2 - y1)^2 + (z2 - z1)^2);
        
        % 判断是否接触
        if distance <= (r1 + r2 + 0.01)
            line([x1, x2], [y1, y2], [z1, z2], 'Color', 'k', 'LineWidth', 1); % 绘制连线
        end
    end
end

% 绘制边界框
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
view(125,15);
hold off; % 释放图形