clc;clear
data=xlsread("data.xlsx")

color=[[1 0 0];[0 0.8 0];[0 0 1];[1 0 1];[0.6 0.4 0.2];[0.5 0 0.5]]


%%画图格式
line_width_box=1;
line_width=3;
symbol="s";
marker_size=7;
line_width_2=1.5;

%%文字格式
font_name="Arial";
legend_font_size=16;
ax_font_size=10;
label_font_size=16;

fig=figure;
set(fig, 'Position', [50, 50, 600, 600]);
hold on

set(gca,"LineWidth",line_width_box,"FontWeight","normal","FontSize",legend_font_size,"FontName",font_name)

a=errorbar(data(:,1),data(:,2), data(:,6),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",color(1,:),"Color",color(1,:),"LineWidth",line_width,"DisplayName","Relative Viscosity"); 
d=errorbar(data(:,1),data(:,10), data(:,11),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",color(2,:),"Color",color(2,:),"LineWidth",line_width,"DisplayName","Relative Viscosity"); 
b=errorbar(data(:,1),data(:,3), data(:,7),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",color(3,:),"Color",color(3,:),"LineWidth",line_width,"DisplayName","Relative Viscosity"); 
e=errorbar(data(:,1),data(:,12), data(:,13),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",color(4,:),"Color",color(4,:),"LineWidth",line_width,"DisplayName","Relative Viscosity"); 
c=errorbar(data(:,1),data(:,4), data(:,8),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",color(5,:),"Color",color(5,:),"LineWidth",line_width,"DisplayName","Relative Viscosity"); 



nnn=3;

x_test=data(1:13,1);
y_test=log10(data(1:13,2));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(0,0.0005,100)
y_fit = polyval(p, x_i)
plot(x_i,10.^y_fit,"--","LineWidth",2.5,"Color",color(1,:)); 

x_test=data(1:12,1);
y_test=log10(data(1:12,4));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(0,0.0005,100)
y_fit = polyval(p, x_i)
plot(x_i,10.^y_fit,"--","LineWidth",2.5,"Color",color(5,:)); 

x_test=data(1:11,1);
y_test=log10(data(1:11,10));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(0,0.00013,20)
y_fit = polyval(p, x_i)
plot(x_i,10.^(y_fit+0.07),"--","LineWidth",2.5,"Color",color(2,:)); 
plot([0.00013 0.0005],[643517,5e5],"--","LineWidth",2.5,"Color",color(2,:)); 

x_test=data(1:13,1);
y_test=log10(data(1:13,12));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(0,0.0005,20)
y_fit = polyval(p, x_i)
%plot(x_i,10.^(y_fit),"--","LineWidth",2.5,"Color",color(4,:)); 

data=xlsread("data_2.xlsx")
plot(data(:,1),data(:,2),"--","LineWidth",2.5,"Color",color(3,:)); 

data=xlsread("data_3.xlsx")
plot(data(:,1),data(:,2),"--","LineWidth",2.5,"Color",color(4,:)); 

xlim([0,0.0005]);
ylim([1.5e4,1.5e9])

plot([1.3e-4 1.3e-4],[3e6 4e7],"--","LineWidth",2.5,"Color","k");

xticks([0,0.0001,0.0002,0.0003,0.0004,0.0005]); 
xticklabels({'0','1','2','3','4','5'});

set(gca, 'YScale', 'log');
l = legend([a, d, b, e, c], '0.620', '0.650', '0.660', '0.675', '0.700', 'Interpreter', 'latex');
set(l, 'box', 'off', 'Orientation', 'horizontal', 'NumColumns', 5);  % 设置为5列

str="$$10^{-4}\ kT_{\rm{eff}}$$";
xlabel(str,"Interpreter","latex",'FontSize',22)

str="$$\mathrm{Viscosity}\ \eta(\phi)$$";
ylabel(str,"Interpreter","latex",'FontSize',22)

annotation('textbox','String','$$(\phi=0.66)$$','Interpreter','latex','FontSize',16,'FontName','Arial','EdgeColor','none', 'HorizontalAlignment', 'center');
annotation('textbox','String','$$T_{\rm eff,c}$$','Interpreter','latex','FontSize',20,'FontName','Arial','EdgeColor','none', 'HorizontalAlignment', 'center');