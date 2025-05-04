clc;clear
data=xlsread("large_active.xlsx")

%%画图格式
line_width_box=1;
line_width=0.5;
marker_width=0.5;

%%文字格式
font_name="Arial";
legend_font_size=16;
ax_font_size=10;
label_font_size=18;
t=40

fig=figure
set(fig, 'Position', [300, 300, 500, 500])

set(gca,"LineWidth",line_width_box,"FontSize",legend_font_size,"FontName",font_name)

colors = jet(size(data, 1));

for i = 1:size(data, 1) - 1
    plot(data(i:i+1, 1), data(i:i+1, 2),'LineWidth', line_width, 'Color', colors(i,:),"Marker","o","MarkerFaceColor",colors(i,:),"MarkerSize",marker_width);
    hold on;
end

colormap(jet);
colorbar;

ax=gca;
set(ax.XAxis, 'FontSize', ax_font_size);
set(ax.YAxis, 'FontSize', ax_font_size);
ax.YColor="k"

if t==40
    xticks([-40,-20,0,20,40]); 
    xticklabels({'-40','20','0','20','40'});
    
    yticks([-40,-20,0,20,40]); 
    yticklabels({'-40','20','0','20','40'});
end

if t==100
    xticks([-100,-50,0,50,100]); 
    xticklabels({'-100','50','0','50','100'});
    
    yticks([-100,-50,0,50,100]); 
    yticklabels({'-100','50','0','50','100'});
end

xlim([-t,t])
ylim([-t,t])

ax = gca;  % 获取当前坐标轴对象
ax.FontName = 'Arial';  % 设置字体类型为Arial
ax.FontSize = 12;       % 设置字体大小为12
axis square;

title('Large Radius, Active', 'FontSize', 16, 'FontName', 'Arial', 'Color', 'k');