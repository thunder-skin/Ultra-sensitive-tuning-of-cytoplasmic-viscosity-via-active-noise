clc;clear

color=[[072 198 235];[129 201 152];[000 082 164];[255 105 180];[251 205 017];[170 000 000];[171 130 255]]/255

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

xx=3e5;
yy1=3e3;
yy2=4e3
nnn=3;

data=xlsread("data_8.xlsx")
d=plot(data(:,1),data(:,3),"o","Marker","x","MarkerSize",15,"MarkerFaceColor",color(2,:),"Color",color(2,:),"LineWidth",2); 
data=xlsread("data_7.xlsx")
e=plot(data(:,1),data(:,3),"o","Marker","x","MarkerSize",15,"MarkerFaceColor",color(5,:),"Color",color(5,:),"LineWidth",2);

%%画散点
data=xlsread("data_2.xlsx")
a=errorbar(data(:,1),data(:,2), data(:,3),data(:,4),"o","Marker","s","MarkerSize",10,"MarkerFaceColor",color(4,:),"Color",color(4,:),"LineWidth",line_width); 
b=errorbar(data(:,1),data(:,6), data(:,7),data(:,8),"o","Marker","s","MarkerSize",10,"MarkerFaceColor",color(1,:),"Color",color(1,:),"LineWidth",line_width); 

pp=1.5;
data=xlsread("data_5.xlsx")
c=errorbar(data(:,1),data(:,8)*pp, data(:,12)*pp,data(:,13)*pp,"o","Marker","^","MarkerSize",10,"MarkerFaceColor",color(7,:),"Color",color(7,:),"LineWidth",line_width); 

data=xlsread("data_8.xlsx")
x_test=log10(data(:,1));
y_test=log10(data(:,3));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(-0.5,3.8,40)
y_fit = polyval(p, x_i)
f=plot(10.^x_i,10.^y_fit,"--","LineWidth",2.5,"Color",color(3,:)); 

data=xlsread("data_7.xlsx")
x_test=log10(data(:,1));
y_test=log10(data(:,3));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(-0.5,3.8,40)
y_fit = polyval(p, x_i)
f=plot(10.^x_i,10.^y_fit,"--","LineWidth",2.5,"Color",color(6,:)); 



%%画拟合曲线
%data=xlsread("data_4.xlsx")
%f=plot(data(:,1),data(:,2),"--","LineWidth",2.5,"Color",color(6,:)); 
%g=plot(data(:,1),data(:,3),"--","LineWidth",2.5,"Color",color(3,:));

xlim([10.^(-1.2),4000])
ylim([5,1550])

xxx=2;
yyy=20;
g=plot([1,100]*xxx,[1 9]*yyy,"--","LineWidth",2.5,"Color",[0 0 0]);

ax=gca;
set(ax.XAxis, 'FontSize', ax_font_size);
set(ax.YAxis, 'FontSize', ax_font_size);

set(gca, 'YScale', 'log',"XScale","log");

l=legend([b,a,c,d,e], ...
    'ATP-depleted Cell (Ebata $$\textit{et al.}$$)','Live Cell (Ebata $$\textit{et al.}$$)','Live Cell (Muenker $$\textit{et al.}$$)','$$kT_{\rm{eff}}=2.5 \times 10^{-5}$$ (Simulation)','$$kT_{\rm{eff}}=1.3 \times 10^{-4}$$ (Simulation)','Interpreter', 'latex');
set(l, 'box', 'off')

str="$$\mathrm{Frequency}\ f\ (\mathrm{Hz})$$";
xlabel(str,"Interpreter","latex",'FontSize',18)

str="$$\mathrm{Storage\ Modulus}\ G'\ (\mathrm{Pa})$$";
ylabel(str,"Interpreter","latex",'FontSize',18,"Rotation",90)