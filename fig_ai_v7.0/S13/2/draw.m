clc;clear

data=xlsread("data_2.xlsx");

color=[[1 0 0];[0 0.8 0];[0 0 1];[0.1 0.8 0.8];[1 0 1];[0.8 0.8 0.1]];

%%画图格式
line_width_box=1;
line_width=1.5;
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

plot(data(:,1),data(:,2),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",color(1,:),"Color",color(1,:),"LineWidth",line_width)
plot(data(:,3),data(:,4),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",color(2,:),"Color",color(2,:),"LineWidth",line_width)
plot(data(:,5),data(:,6),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",color(3,:),"Color",color(3,:),"LineWidth",line_width)
plot(data(:,7),data(:,8),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",color(4,:),"Color",color(4,:),"LineWidth",line_width)
plot(data(:,9),data(:,10),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",color(5,:),"Color",color(5,:),"LineWidth",line_width)


data=xlsread("data_7.xlsx");
plot(data(:,1),data(:,2),"--","Color",color(1,:),"LineWidth",line_width)
plot(data(:,3),data(:,4),"--","Color",color(2,:),"LineWidth",line_width)
plot(data(:,5),data(:,6),"--","Color",color(3,:),"LineWidth",line_width)
plot(data(:,7),data(:,8),"--","Color",color(4,:),"LineWidth",line_width)
plot(data(:,9),data(:,10),"--","Color",color(5,:),"LineWidth",line_width)


plot([0.841 0.841],[6e3 1e10],"--","Color","k","LineWidth",line_width*0.7)
plot([0.86 0.86],[6e3 1e10],"--","Color","k","LineWidth",line_width*0.7)

ax=gca;
set(ax.XAxis, 'FontSize', ax_font_size);
set(ax.YAxis, 'FontSize', ax_font_size);

h=legend("Athermal","$$1.0\times 10^{-5}$$","$$3.2 \times 10^{-5}$$","$$1.0\times 10^{-4}$$","$$1.8\times 10^{-4}$$","Location","eastoutside",'Interpreter', 'latex')
set(h, 'Box', 'off')

xlim([0.818,0.882])
ylim([6e3,1.8e9])

str="$$\mathrm{Volume\ Fraction}\ \phi$$";
xlabel(str,"Interpreter","latex",'FontSize',18)

str="$$\mathrm{Viscosity}\ \eta(\phi)$$";
ylabel(str,"Interpreter","latex",'FontSize',18)

set(gca, 'YScale', 'log');

%annotation("textbox","string",'Fluid', "FontSize",18,"FontName",font_name,'FontWeight', 'normal','EdgeColor','none',"Color","b")
annotation("textbox","string",{'Sensitive','Regime'}, "FontSize",16,"FontName",font_name,'FontWeight', 'normal','EdgeColor','none',"Color","k",'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle')
%annotation("textbox","string",'Jam', "FontSize",18,"FontName",font_name,'FontWeight', 'normal','EdgeColor','none',"Color","r")

