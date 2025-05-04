% 创建一个新的图形窗口
figure;

circles=xlsread("data.xlsx");

rm=1.4;

% 设置空间范围
axis([-rm 10+rm -rm 10+rm -rm 10+rm]); % 设置坐标轴范围为10x10x10的空间
hold on; % 保持当前图形

% 绘制深蓝色的球体 (半径 = 1.4)

for i=1:708
    x = circles(i, 1);
    y = circles(i, 2);
    z = circles(i, 3);
    r = circles(i, 4);

    if x<10 && y<10 && z<10 
        [X1, Y1, Z1] = sphere(200); % 生成球体的坐标
        if r<1.2
            surf(r*X1+x, r*Y1+y, r*Z1+z, 'FaceColor', [0.4 0.6 1.0], 'EdgeColor', 'none'); % 蓝色
        end

        if r>1.3
            surf(r*X1+x, r*Y1+y, r*Z1+z, 'FaceColor', [1.0 0.5 0.7], 'EdgeColor', 'none'); % 红色
        end
    end
end

line([0 10], [0 0], [0,0], Color = [0 0 0], LineStyle = "-.", LineWidth = 2)
line([0 10], [10 10], [0,0], Color = [0 0 0], LineStyle = "-.", LineWidth = 2)
line([0 10], [10 10], [10,10], Color = [0 0 0], LineStyle = "-.", LineWidth = 2)
line([0 10], [0 0], [10,10], Color = [0 0 0], LineStyle = "-.", LineWidth = 2)
line([0 0], [10 0], [0,0], Color = [0 0 0], LineStyle = "-.", LineWidth = 2)
line([0 0], [10 0], [10,10], Color = [0 0 0], LineStyle = "-.", LineWidth = 2)
line([10 10], [10 0], [10,10], Color = [0 0 0], LineStyle = "-.", LineWidth = 2)
line([10 10], [10 0], [0,0], Color = [0 0 0], LineStyle = "-.", LineWidth = 2)
line([0 0], [0 0], [0,10], Color = [0 0 0], LineStyle = "-.", LineWidth = 2)
line([10 10], [0 0], [0,10], Color = [0 0 0], LineStyle = "-.", LineWidth = 2)
line([10 10], [10 10], [0,10], Color = [0 0 0], LineStyle = "-.", LineWidth = 2)
line([0 0], [10 10], [0,10], Color = [0 0 0], LineStyle = "-.", LineWidth = 2)


% 添加光照
camlight; % 添加摄像机光源
lighting gouraud; % 使用 Gouraud 光照模型
material shiny; % 设置材质为光亮

% 添加图例和坐标轴标签

grid on; % 显示网格

xticks(0:2:10);
yticks(0:2:10);
zticks(0:2:10);

%axis equal; % 设置坐标轴等比例
view(3); % 以三维视角查看
hold off; % 释放图形
