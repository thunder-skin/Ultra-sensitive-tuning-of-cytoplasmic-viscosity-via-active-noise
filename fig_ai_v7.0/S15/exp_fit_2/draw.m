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
legend_font_size=15;
ax_font_size=14;
label_font_size=18;

xx=10;
yy=3;

fig=figure;
set(fig, 'Position', [50, 50, 600, 500]);
hold on

set(gca,"LineWidth",line_width_box,"FontSize",legend_font_size,"FontName",font_name)

ky=2.4;

%%画散点
data=xlsread("data_1.xlsx")
a=errorbar(data(:,1),data(:,3)*ky, data(:,4)*ky,data(:,5)*ky,"o","Marker","s","MarkerSize",10,"MarkerFaceColor",color(4,:),"Color",color(4,:),"LineWidth",line_width); 
b=errorbar(data(:,1),data(:,7)*ky, data(:,8)*ky,data(:,9)*ky,"o","Marker","s","MarkerSize",10,"MarkerFaceColor",color(1,:),"Color",color(1,:),"LineWidth",line_width); 

data=xlsread("data_2.xlsx")
e=plot(data(:,1)*xx,data(:,2)*yy,"o","Marker","x","MarkerSize",15,"MarkerFaceColor",color(2,:),"Color",color(2,:),"LineWidth",2); 
d=plot(data(:,1)*xx,data(:,3)*yy*0.85,"o","Marker","x","MarkerSize",15,"MarkerFaceColor",color(5,:),"Color",color(5,:),"LineWidth",2);

f=plot([0.07,3000],[14,900],"--","LineWidth",2.5,"Color",color(6,:)); 
g=plot([0.07,3000],[4.5,1000],"--","LineWidth",2.5,"Color",color(3,:));

ax=gca;
set(ax.XAxis, 'FontSize', ax_font_size);
set(ax.YAxis, 'FontSize', ax_font_size);

pp=1.5;
data=xlsread("data_5.xlsx")
c=errorbar(data(:,1),data(:,3)*pp, data(:,6)*pp,data(:,7)*pp,"o","Marker","^","MarkerSize",10,"MarkerFaceColor",color(7,:),"Color",color(7,:),"LineWidth",line_width); 

xxx=2;
yyy=15;
g=plot([1,100]*xxx,[1 9]*yyy,"--","LineWidth",2.5,"Color",[0 0 0]);

xlim([10.^(-1.2),4000])
ylim([3.75,1500])

set(gca, 'YScale', 'log',"XScale","log");

l=legend([c,a,b,d,e],'ATP-depleted Cell (Ebata $$\textit{et al.}$$)','Live Cell (Ebata $$\textit{et al.}$$)','Live Cell (Muenker $$\textit{et al.}$$)','$$kT_{\rm{eff}}=3.2 \times 10^{-4}$$ (Simulation)','$$kT_{\rm{eff}}=1.8 \times 10^{-3}$$ (Simulation)','Interpreter', 'latex');
set(l, 'box', 'off')

str="$$\mathrm{Frequency}\ f\ (\mathrm{Hz})$$";
xlabel(str,"Interpreter","latex",'FontSize',18)

str="$$\mathrm{Vicousity\ Modulus}\ G''\ (\mathrm{Pa})$$";
ylabel(str,"Interpreter","latex",'FontSize',18,"Rotation",90)