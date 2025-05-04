clc;clear
data=xlsread("data_1.xlsx");

color=[[1 0 0];[0 0.8 0];[0 0 1];[0.1 0.8 0.8];[1 0 1];[0.8 0.8 0.1];[0.6 0.4 0.2];[0.5 0 0.5]]

%%画图格式
line_width_box=1;
line_width=2;
symbol="s";
marker_size=5;
mk=5;
line_width_2=1.5;

%%文字格式
font_name="Arial";
legend_font_size=18;
ax_font_size=14;
label_font_size=26;

fig=figure;
set(fig, 'Position', [50, 50, 600, 600]);
hold on
box on

set(gca,"LineWidth",line_width_box,"FontWeight","normal","FontSize",legend_font_size,"FontName",font_name)

%描点
errorbar(data(:,1),data(:,2), data(:,3),data(:,3),"o","Marker","s","MarkerSize",mk,"MarkerFaceColor",color(1,:),"Color",color(1,:),"LineWidth",line_width); 
errorbar(data(:,1),data(:,4), data(:,5),data(:,5),"o","Marker","s","MarkerSize",mk,"MarkerFaceColor",color(2,:),"Color",color(2,:),"LineWidth",line_width); 
errorbar(data(:,1),data(:,6), data(:,7),data(:,7),"o","Marker","s","MarkerSize",mk,"MarkerFaceColor",color(3,:),"Color",color(3,:),"LineWidth",line_width); 
%errorbar(data(:,1),data(:,8), data(:,9),data(:,9),"o","Marker","s","MarkerSize",mk,"MarkerFaceColor",color(4,:),"Color",color(4,:),"LineWidth",line_width); 
errorbar(data(:,1),data(:,10), data(:,11),data(:,11),"o","Marker","s","MarkerSize",mk,"MarkerFaceColor",color(5,:),"Color",color(5,:),"LineWidth",line_width); 
errorbar(data(:,1),data(:,12), data(:,13),data(:,13),"o","Marker","s","MarkerSize",mk,"MarkerFaceColor",color(6,:),"Color",color(6,:),"LineWidth",line_width); 
errorbar(data(:,1),data(:,14), data(:,15),data(:,15),"o","Marker","s","MarkerSize",mk,"MarkerFaceColor",color(7,:),"Color",color(7,:),"LineWidth",line_width); 
errorbar(data(:,1),data(:,16), data(:,17),data(:,17),"o","Marker","s","MarkerSize",mk,"MarkerFaceColor",color(8,:),"Color",color(8,:),"LineWidth",line_width); 

data=xlsread("data_4.xlsx");

plot(data(:,1),data(:,2),"--","Color",color(1,:),"LineWidth",line_width);
plot(data(:,3),data(:,4),"--","Color",color(2,:),"LineWidth",line_width);
plot(data(:,5),data(:,6),"--","Color",color(3,:),"LineWidth",line_width);
%plot(data(:,7),data(:,8),"--","Color",color(4,:),"LineWidth",line_width);
plot(data(:,9),data(:,10),"--","Color",color(5,:),"LineWidth",line_width);
plot(data(:,11),data(:,12),"--","Color",color(6,:),"LineWidth",line_width);
plot(data(:,13),data(:,14),"--","Color",color(7,:),"LineWidth",line_width);
plot(data(:,15),data(:,16),"--","Color",color(8,:),"LineWidth",line_width);


plot([0.65 0.65],[6e3 1e10],"--","Color","k","LineWidth",line_width)
plot([0.675 0.675],[6e3 1e10],"--","Color","k","LineWidth",line_width)

ax=gca;
set(ax.XAxis, 'FontSize', ax_font_size);
set(ax.YAxis, 'FontSize', ax_font_size);

h = legend(...
    'Athermal', ...
    '$$1.0 \times 10^{-6}$$', ...
    '$$1.0 \times 10^{-5}$$', ...
    '$$1.0 \times 10^{-4}$$', ...
    '$$1.8 \times 10^{-4}$$', ...
    '$$3.2 \times 10^{-4}$$', ...
    '$$5.6 \times 10^{-4}$$', ...
    'Interpreter', 'latex');
set(h, 'Box', 'off')



xlim([0.615,0.705])
ylim([2e4,1e9])

str="$$\mathrm{Volume\ Fraction}\ \phi$$";
xlabel(str,"Interpreter","latex",'FontSize',22)

str="$$\mathrm{Viscosity}\ \eta(\phi)$$";
ylabel(str,"Interpreter","latex",'FontSize',22)

set(gca, 'YScale', 'log');

%annotation("textbox","string",'Fluid', "FontSize",18,"FontName",font_name,'FontWeight', 'normal','EdgeColor','none',"Color","b")
annotation("textbox","string",'Sensitive', "FontSize",18,"FontName",font_name,'FontWeight', 'normal','EdgeColor','none',"Color","k")
annotation("textbox","string",'Regime', "FontSize",18,"FontName",font_name,'FontWeight', 'normal','EdgeColor','none',"Color","k")
%annotation("textbox","string",'Solid', "FontSize",18,"FontName",font_name,'FontWeight', 'normal','EdgeColor','none',"Color","r")

