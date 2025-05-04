clc;clear

color=[[255 105 180];[251 205 017];[170 000 000];[072 198 235];[129 201 152];[000 082 164];[171 130 255]]/255

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

xx=3e5;
yy=5e3;

fig=figure;
set(fig, 'Position', [50, 50, 600, 500]);
hold on

set(gca,"LineWidth",line_width_box,"FontSize",legend_font_size,"FontName",font_name)

ky=2.4;

p1=15
p2=9
yy2=4.5e3
nnn=3

data=xlsread("data_8.xlsx")
c=plot(data(:,5)*xx,data(:,2)*yy-p1,"o","Marker","x","MarkerSize",15,"MarkerFaceColor",color(5,:),"Color",color(5,:),"LineWidth",2); 
data=xlsread("data_7.xlsx")
d=plot(data(:,1)*xx,data(:,3)*yy2-p2,"o","Marker","x","MarkerSize",15,"MarkerFaceColor",color(2,:),"Color",color(2,:),"LineWidth",2);


%%画散点
data=xlsread("data_1.xlsx")
b=errorbar(data(:,1),data(:,3)*ky, data(:,4)*ky,data(:,5)*ky,"o","Marker","s","MarkerSize",10,"MarkerFaceColor",color(4,:),"Color",color(4,:),"LineWidth",line_width); 
a=errorbar(data(:,1),data(:,7)*ky, data(:,8)*ky,data(:,9)*ky,"o","Marker","s","MarkerSize",10,"MarkerFaceColor",color(1,:),"Color",color(1,:),"LineWidth",line_width); 

data=xlsread("data_8.xlsx")
x_test=log10(data(:,5)*xx);
y_test=log10(data(:,2)*yy-p1);
p = polyfit(x_test, y_test, 3);
x_i=linspace(-0.5,3.8,40)
y_fit=csaps(x_test,y_test,0.5,x_i);

f=plot(10.^x_i,10.^y_fit,"--","LineWidth",2.5,"Color",color(6,:)); 

data=xlsread("data_7.xlsx")
x_test=log10(data(:,1)*xx);
y_test=log10(data(:,3)*yy2-p2);
p = polyfit(x_test, y_test, 3);
x_i=linspace(-0.5,3.8,40)
y_fit = polyval(p, x_i)
f=plot(10.^x_i,10.^y_fit*0.97,"--","LineWidth",2.5,"Color",color(3,:)); 

pp=1.5;
data=xlsread("data_5.xlsx")
e=errorbar(data(:,1),data(:,3)*pp, data(:,6)*pp,data(:,7)*pp,"o","Marker","^","MarkerSize",10,"MarkerFaceColor",color(7,:),"Color",color(7,:),"LineWidth",line_width); 




ax=gca;
set(ax.XAxis, 'FontSize', ax_font_size);
set(ax.YAxis, 'FontSize', ax_font_size);

xlim([10.^(-1.2),4000])
ylim([4,1500])

xxx=2;
yyy=15;
g=plot([1,100]*xxx,[1 9]*yyy,"--","LineWidth",2.5,"Color",[0 0 0]);

set(gca, 'YScale', 'log',"XScale","log");

l=legend([b,a,e,c,d], ...
    'ATP-depleted Cell (Ebata $$\textit{et al.}$$)','Live Cell (Ebata $$\textit{et al.}$$)','Live Cell (Muenker $$\textit{et al.}$$)','$$T_{\rm{eff}}=2.2 \times 10^{-5}$$ (Simulation)','$$T_{\rm{eff}}=1.6 \times 10^{-4}$$ (Simulation)','Interpreter', 'latex');
set(l, 'box', 'off')

str="$$\mathrm{Frequency}\ f\ (\mathrm{Hz})$$";
xlabel(str,"Interpreter","latex",'FontSize',18)

str="$$\mathrm{Loss\ Modulus}\ G''\ (\mathrm{Pa})$$";
ylabel(str,"Interpreter","latex",'FontSize',18,"Rotation",90)