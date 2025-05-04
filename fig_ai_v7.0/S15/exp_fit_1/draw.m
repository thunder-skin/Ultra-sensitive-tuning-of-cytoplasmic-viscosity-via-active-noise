clc;clear
data=xlsread("data_2.xlsx")
color=[[072 198 235];[129 201 152];[000 082 164];[255 105 180];[251 205 017];[170 000 000];[171 130 255]]/255

%%画图格式
line_width_box=1;
line_width=1.5;
symbol="s";
marker_size=9;
line_width_2=1;

%%文字格式
font_name="Arial";
legend_font_size=15;
ax_font_size=14;
label_font_size=18;

xx=10;
yy=3;

fig=figure;
set(fig, 'Position', [50, 50, 600, 500]);

hold on

set(gca,"LineWidth",line_width_box,"FontSize",legend_font_size,"FontName",font_name)

%%画散点
a=errorbar(data(:,1),data(:,2), data(:,3),data(:,4),"o","Marker","s","MarkerSize",10,"MarkerFaceColor",color(4,:),"Color",color(4,:),"LineWidth",line_width); 
b=errorbar(data(:,1),data(:,6), data(:,7),data(:,8),"o","Marker","s","MarkerSize",10,"MarkerFaceColor",color(1,:),"Color",color(1,:),"LineWidth",line_width); 

pp=7.5;
data=xlsread("data_5.xlsx")
c=errorbar(data(:,1),data(:,2)*pp, (data(:,4)-data(:,2))*pp,(data(:,3)-data(:,2))*pp,"o","Marker","^","MarkerSize",10,"MarkerFaceColor",color(7,:),"Color",color(7,:),"LineWidth",line_width); 

data=xlsread("data_3.xlsx")
d=plot(data(:,3)*xx,data(:,4)*yy,"o","Marker","x","MarkerSize",15,"MarkerFaceColor",color(5,:),"Color",color(5,:),"LineWidth",2); 
e=plot(data(:,3)*xx,data(:,5)*yy*0.85,"o","Marker","x","MarkerSize",15,"MarkerFaceColor",color(2,:),"Color",color(2,:),"LineWidth",2);

data=xlsread("data_4.xlsx")
f=plot(data(:,1),data(:,2),"--","LineWidth",2.5,"Color",color(6,:)); 
g=plot(data(:,1),data(:,3),"--","LineWidth",2.5,"Color",color(3,:));

xlim([10.^(-1.2),4000])
ylim([5,1550])

ax=gca;
set(ax.XAxis, 'FontSize', ax_font_size);
set(ax.YAxis, 'FontSize', ax_font_size);

set(gca, 'YScale', 'log',"XScale","log");

xxx=2;
yyy=20;
g=plot([1,100]*xxx,[1 9]*yyy,"--","LineWidth",2.5,"Color",[0 0 0]);

l=legend([c,a,b,d,e],'ATP-depleted Cell (Ebata $$\textit{et al.}$$)','Live Cell (Ebata $$\textit{et al.}$$)','Live Cell (Muenker $$\textit{et al.}$$)','$$kT_{\rm{eff}}=3.2 \times 10^{-4}$$ (Simulation)','$$kT_{\rm{eff}}=1.8 \times 10^{-3}$$ (Simulation)','Interpreter', 'latex');
set(l, 'box', 'off')

str="$$\mathrm{Frequency}\ f\ (\mathrm{Hz})$$";
xlabel(str,"Interpreter","latex",'FontSize',18)

str="$$\mathrm{Storage\ Modulus}\ G'\ (\mathrm{Pa})$$";
ylabel(str,"Interpreter","latex",'FontSize',18,"Rotation",90)