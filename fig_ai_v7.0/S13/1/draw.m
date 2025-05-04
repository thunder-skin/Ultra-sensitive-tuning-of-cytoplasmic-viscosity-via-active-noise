clc;clear
data=xlsread("data_1.xlsx");

color=[[1 0 0];[0 0.8 0];[0 0 1];[0.1 0.8 0.8];[1 0 1]];

%%画图格式
line_width_box=1;
line_width=2;
symbol="s";
marker_size=7;
line_width_2=1.5;

%%文字格式
font_name="Arial";
legend_font_size=18;
ax_font_size=14;
label_font_size=18;

fig=figure;
set(fig, 'Position', [50, 50, 600, 500]);
hold on
box on

set(gca,"LineWidth",line_width_box,"FontSize",legend_font_size,"FontName",font_name)
set(gca,"LineWidth",line_width_box)
yyy=1.66
mk=6

errorbar(data(:,1)/0.8415,yyy*data(:,2)/2700000, yyy*data(:,27)/2700000,yyy*data(:,27)/2700000,"o","Marker","s","MarkerSize",mk,"MarkerFaceColor",color(1,:),"Color",color(1,:),"LineWidth",line_width); 
errorbar(data(:,1)/0.844,yyy*data(:,3)/12000000, yyy*data(:,28)/12000000,yyy*data(:,28)/12000000,"o","Marker","s","MarkerSize",mk,"MarkerFaceColor",color(2,:),"Color",color(2,:),"LineWidth",line_width); 
errorbar(data(:,4)/0.846,yyy*data(:,5)/21000000, yyy*data(:,29)/21000000,yyy*data(:,29)/21000000,"o","Marker","s","MarkerSize",mk,"MarkerFaceColor",color(3,:),"Color",color(3,:),"LineWidth",line_width); 
errorbar(data(:,6)/0.85,yyy*data(:,7)/42000000, yyy*data(:,30)/43000000,yyy*data(:,30)/43000000,"o","Marker","s","MarkerSize",mk,"MarkerFaceColor",color(4,:),"Color",color(4,:),"LineWidth",line_width); 
errorbar(data(:,6)/0.852,yyy*data(:,8)/48000000, yyy*data(:,31)/48000000,yyy*data(:,31)/48000000,"o","Marker","s","MarkerSize",mk,"MarkerFaceColor",color(5,:),"Color",color(5,:),"LineWidth",line_width); 


data=xlsread("data_2.xlsx");

plot(data(:,1)/0.8415,yyy*data(:,2)/2700000,"--","Color",color(1,:),"LineWidth",line_width)
plot(data(:,3)/0.844,yyy*data(:,4)/12000000,"--","Color",color(2,:),"LineWidth",line_width)
plot(data(:,5)/0.846,yyy*data(:,6)/21000000,"--","Color",color(3,:),"LineWidth",line_width)
plot(data(:,7)/0.85,yyy*data(:,8)/42000000,"--","Color",color(4,:),"LineWidth",line_width)
plot(data(:,9)/0.852,yyy*data(:,10)/48000000,"--","Color",color(5,:),"LineWidth",line_width)

plot([0.965 1],[1 1],"--","Color","k","LineWidth",line_width*0.5)
plot([1 1],[1e-3 1],"--","Color","k","LineWidth",line_width*0.5)

ax=gca;
set(ax.XAxis, 'FontSize', ax_font_size);
set(ax.YAxis, 'FontSize', ax_font_size);

h=legend("Athermal","$$1.0\times 10^{-5}$$","$$3.2 \times 10^{-5}$$","$$1.0\times 10^{-4}$$","$$1.8\times 10^{-4}$$","Location","eastoutside",'Interpreter', 'latex')
set(h, 'Box', 'off')

xlim([0.965,1.03])
ylim([3e-3,7e2])

str="$$\mathrm{Normalized\ Volume\ Fraction}\ \phi/\phi_c$$";
xlabel(str,"Interpreter","latex",'FontSize',18)

str="$$\mathrm{Normalized\ Viscosity}\ \eta(\phi)/\eta(\phi_c)$$";
ylabel(str,"Interpreter","latex",'FontSize',18)


set(gca, 'YScale', 'log');

