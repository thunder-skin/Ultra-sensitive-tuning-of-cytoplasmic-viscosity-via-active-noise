% 定义参数
p = 1;
s = 0.04;
c = 0.5;

font_name = "Arial";

% 定义分段函数
f = @(x) (x < s) .* (p/2 * (x.^2 - c*s)) + ...
         (x >= s & x < c) .* (p/2 * s / (s - c) * (c - x).^2);

% 定义平移后的函数（向右平移2个单位）
f_shifted = @(x) f(x - 2);

% 定义绘图范围（x=1.4~2.8）
x = linspace(1.6, 2.8, 1000);
y = f_shifted(x);

% 设置全局字体
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');

% 绘制图形
figure;
plot(x, y, 'b', 'LineWidth', 2);

title('Interaction Potential Function', 'FontSize', 18);
grid on;
xlim([1.4, 2.8]);
ylim([-0.02, 0.1]);

% 添加垂直虚线
line([2 2], [-0.02, 0.1], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1);
line([2.5 2.5], [-0.02 0.1], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1);

% 定义x轴刻度值（按升序排列）
x_ticks = sort([1.5, 1.75, 2, 2.25, 2.5, 2.75]);

% 设置x轴刻度
xticks(x_ticks);

% 设置x轴刻度标签（水平显示）
xticklabels({'1.5', '1.75', '2', '2.25', '2.5', '2.75'});
ax = gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.XAxis.TickLabelRotation = 0; % 确保刻度标签水平

% 添加数据点连线
hold on;
scatter(2.04, f_shifted(2.04), 50, 'r', 'filled');
scatter(2.5, f_shifted(2.5), 50, 'r', 'filled');

% 在红点旁边添加 LaTeX 文字
text(2.04, f_shifted(2.04) + 0.005, '$r_s=0.04$', 'Interpreter', 'latex', 'FontSize', 16, 'HorizontalAlignment', 'left');
text(2.5, f_shifted(2.5) + 0.005, '$r_c=0.5$', 'Interpreter', 'latex', 'FontSize', 16, 'HorizontalAlignment', 'left');

% 添加区域标签
text(1.70, 0.09, 'Harmonic Repulsion', 'HorizontalAlignment', 'center', 'FontSize', 18);
text(2.25, 0.09, 'Attraction', 'HorizontalAlignment', 'center', 'FontSize', 18);
text(2.65, 0.09, 'Free', 'HorizontalAlignment', 'center', 'FontSize', 18);

% 调整图形布局
set(gcf, 'Position', [100, 100, 800, 600]);

% 设置 xlabel 和 ylabel
str = "$$\mathrm{Sphere\ Center\ Distance}\ r_{ij}$$";
xlabel(str, "Interpreter", "latex", 'FontSize', 20);

str = "$$\mathrm{Interaction\ Potential}\ U(r)$$";
ylabel(str, "Interpreter", "latex", 'FontSize', 20, "Rotation", 90);% 定义参数
p = 1;
s = 0.04;
c = 0.5;

font_name = "Arial";

% 定义分段函数
f = @(x) (x < s) .* (p/2 * (x.^2 - c*s)) + ...
         (x >= s & x < c) .* (p/2 * s / (s - c) * (c - x).^2);

% 定义平移后的函数（向右平移2个单位）
f_shifted = @(x) f(x - 2);

% 定义绘图范围（x=1.4~2.8）
x = linspace(1.575, 2.8, 1000);
y = f_shifted(x);

% 设置全局字体
set(groot, 'defaultAxesFontName', 'Arial');
set(groot, 'defaultTextFontName', 'Arial');

% 绘制图形
figure;
plot(x, y, 'b', 'LineWidth', 2);

title('Interaction Potential Function', 'FontSize', 18);
grid on;
xlim([1.5, 2.75]);
ylim([-0.02, 0.1]);

% 添加垂直虚线
line([2 2], [-0.02, 0.1], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1);
line([2.5 2.5], [-0.02 0.1], 'LineStyle', '--', 'Color', 'k', 'LineWidth', 1);

line([1.4 2.8], [0, 0], 'LineStyle', '-', 'Color', 'k', 'LineWidth', 1);

% 定义x轴刻度值（按升序排列）
x_ticks = sort([1.5, 1.75, 2, 2.25, 2.5, 2.75]);

% 设置x轴刻度
xticks(x_ticks);

% 设置x轴刻度标签（水平显示）
xticklabels({'1.5', '1.75', '2', '2.25', '2.5', '2.75'});
ax = gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;
ax.XAxis.TickLabelRotation = 0; % 确保刻度标签水平

% 添加数据点连线
hold on;
scatter(2.04, f_shifted(2.04), 50, 'r', 'filled');
scatter(2.5, f_shifted(2.5), 50, 'r', 'filled');

% 在红点旁边添加 LaTeX 文字
text(2.01,0.03, '$r_s=0.04$', 'Interpreter', 'latex', 'FontSize', 16, 'HorizontalAlignment', 'left',Color='r');
text(2.01,0.048, '$r_c=0.5$', 'Interpreter', 'latex', 'FontSize', 16, 'HorizontalAlignment', 'left',Color='r');

% 添加区域标签
text(1.75, 0.09, 'Harmonic Repulsion', 'HorizontalAlignment', 'center', 'FontSize', 18);
text(2.25, 0.09, 'Attraction', 'HorizontalAlignment', 'center', 'FontSize', 18);
text(2.625, 0.09, 'Free', 'HorizontalAlignment', 'center', 'FontSize', 18);

% 调整图形布局
set(gcf, 'Position', [100, 100, 1000, 600]);

% 添加竖直虚线及两端箭头
% 对于 x=2.04
line([2.04, 2.04], [f_shifted(2.04), 0.025], 'LineStyle', '--', 'Color', 'r', 'LineWidth', 3); % 竖直虚线

x_1=0.44
x_2=0.465
y_0=0.42

y_1=0.54

annotation('arrow', [0.45, 0.465], [y_0, y_0], ...
           'HeadLength', 7, 'HeadWidth', 7, 'LineWidth', 2, 'Color', 'r');

annotation('arrow', [0.455, 0.44], [y_0, y_0], ...
           'HeadLength', 7, 'HeadWidth', 7, 'LineWidth', 2, 'Color', 'r');

annotation('arrow', [0.45, 0.75], [y_1, y_1], ...
           'HeadLength', 7, 'HeadWidth', 7, 'LineWidth', 2, 'Color', 'r');

annotation('arrow', [0.455, 0.44], [y_1, y_1], ...
           'HeadLength', 7, 'HeadWidth', 7, 'LineWidth', 2, 'Color', 'r');


% 对于 x=2.5
line([2.5, 2.5], [f_shifted(2.5), 0.045], 'LineStyle', '--', 'Color', 'r', 'LineWidth', 3); % 竖直虚线


% 设置 xlabel 和 ylabel
str = "$$\mathrm{Sphere\ Center\ Distance}\ r_{ij}$$";
xlabel(str, "Interpreter", "latex", 'FontSize', 20);

str = "$$\mathrm{Interaction\ Potential}\ U(r)$$";
ylabel(str, "Interpreter", "latex", 'FontSize', 20, "Rotation", 90);