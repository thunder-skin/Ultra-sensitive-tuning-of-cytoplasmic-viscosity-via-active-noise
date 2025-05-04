clc;clear

color=[[072 198 235];[100, 220, 100];[000 082 164];[255 105 180];[251 205 017];[170 000 000];[171 130 255];[0, 150, 0];[255 165 0]]/255

%%画图格式
line_width_box=1;
line_width=1.5;
symbol="s";
marker_size=9;
line_width_2=1;

%%文字格式
font_name="Arial";
legend_font_size=14;
ax_font_size=14;
label_font_size=18;



fig=figure;
set(fig, 'Position', [50, 50, 600, 500]);

hold on

set(gca,"LineWidth",line_width_box,"FontSize",legend_font_size,"FontName",font_name)

xx=3e4;
yy1=3e3;
yy2=1.5e3
yy3=2.5e3
nnn=3;

data=xlsread("data_1.xlsx")
d=plot(data(:,1)*xx,data(:,4)*yy1,"o","Marker","x","MarkerSize",15,"MarkerFaceColor",color(2,:),"Color",color(2,:),"LineWidth",2); 

k1=plot(data(:,1)*xx,data(:,7)*yy3,"o","Marker","x","MarkerSize",15,"MarkerFaceColor",color(2,:),"Color",color(8,:),"LineWidth",2); 

e=plot(data(:,1)*xx,data(:,5)*yy2,"o","Marker","x","MarkerSize",15,"MarkerFaceColor",color(5,:),"Color",color(5,:),"LineWidth",2);

p1=plot(data(:,1)*xx,data(:,6)*yy2,"o","Marker","x","MarkerSize",15,"MarkerFaceColor",color(5,:),"Color",color(9,:),"LineWidth",2);

%%画散点
data=xlsread("data_2.xlsx")
a=errorbar(data(:,1),data(:,2), data(:,3),data(:,4),"o","Marker","s","MarkerSize",10,"MarkerFaceColor",color(4,:),"Color",color(4,:),"LineWidth",line_width); 
b=errorbar(data(:,1),data(:,6), data(:,7),data(:,8),"o","Marker","s","MarkerSize",10,"MarkerFaceColor",color(1,:),"Color",color(1,:),"LineWidth",line_width); 

pp=1.5;
data=xlsread("data_5.xlsx")
c=errorbar(data(:,1),data(:,8)*pp, data(:,12)*pp,data(:,13)*pp,"o","Marker","^","MarkerSize",10,"MarkerFaceColor",color(7,:),"Color",color(7,:),"LineWidth",line_width); 


data=xlsread("data_9.xlsx");
x_test=log10(data(:,1)*xx);
y_test=log10(data(:,4)*yy1);
p = polyfit(x_test, y_test, nnn);
x_i=linspace(-0.7,3.8,40)
y_fit = polyval(p, x_i)
f=plot(10.^x_i,10.^y_fit,"--","LineWidth",2.5,"Color",color(3,:)); 

x_test=log10(data(1:42,1)*xx);
y_test=log10(data(1:42,7)*yy3);
p = polyfit(x_test, y_test, nnn);
x_i=linspace(-0.6,3.8,40)
y_fit = polyval(p, x_i)
f=plot(10.^x_i,10.^y_fit,"--","LineWidth",2.5,"Color",color(3,:)); 

x_test=log10(data(1:36,1)*xx);
y_test=log10(data(1:36,5)*yy2);
p = polyfit(x_test, y_test, nnn);
x_i=linspace(-0.5,3.8,40)
y_fit = polyval(p, x_i)
f=plot(10.^x_i,10.^y_fit,"--","LineWidth",2.5,"Color",color(6,:)); 

x_test=log10(data(1:41,1)*xx);
y_test=log10(data(1:41,6)*yy2);
p = polyfit(x_test, y_test, nnn);
x_i=linspace(-0.5,3.8,40)
y_fit = polyval(p, x_i)
f=plot(10.^x_i,10.^y_fit,"--","LineWidth",2.5,"Color",color(6,:)); 



xlim([10.^(-0.9),8000])
ylim([2,2200])

xxx=0.7;
yyy=35;
g=plot([1,100]*xxx,[1 10]*yyy,"--","LineWidth",2.5,"Color",[0 0 0]);

ax=gca;
set(ax.XAxis, 'FontSize', ax_font_size);
set(ax.YAxis, 'FontSize', ax_font_size);

set(gca, 'YScale', 'log',"XScale","log");

l=legend([b,a,c,k1,d,e,p1], ...
    'ATP-depleted Cell (Ebata $$\textit{et al.}$$)','Live Cell (Ebata $$\textit{et al.}$$)','Live Cell (Muenker $$\textit{et al.}$$)','$$\phi=0.75\ kT_{\rm eff}=0$$','$$\phi=0.75\ kT_{\rm eff}=1.0 \times 10^{-3}$$','$$\phi=0.62\ kT_{\rm eff}=0$$','$$\phi=0.62\ kT_{\rm eff}=3.8\times 10^{-4}$$','Interpreter', 'latex');
set(l, 'box', 'off')

str="$$\mathrm{Frequency}\ f\ (\mathrm{Hz})$$";
xlabel(str,"Interpreter","latex",'FontSize',18)

str="$$\mathrm{Storage\ Modulus}\ G'\ (\mathrm{Pa})$$";
ylabel(str,"Interpreter","latex",'FontSize',18,"Rotation",90)