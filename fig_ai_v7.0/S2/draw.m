clc;clear
data=xlsread("calculation.xlsx")

%%画图格式
line_width_box=1;
line_width=2;
marker_size=12;

%%文字格式
font_name="Arial"
legend_font_size=18;
ax_font_size=10;
label_font_size=18;

fig=figure;
set(fig, 'Position', [50, 50, 900, 400]);

hold on
box on

rem=0.011

plot(data(:,1),(data(:,3)+rem)/(1+rem),"-","Color","k","LineWid",3)
plot(data(:,5)/3,data(:,6),"o","Color","r","LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor","r")
plot(data(:,8)/2,data(:,9),"o","Color","b","LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor","b")


set(gca,"LineWidth",line_width_box)
x1=1;
x2=2e4;
y1=8e-3;
y2=1.3;
a=plot([x1,x1],[y1,y2],"--","Color","k","LineWidth",line_width*0.8)

b=plot([x2,x2],[y1,0.3],"--","Color","k","LineWidth",line_width*0.8)

xlim([0.03,7e5])
ylim([y1,y2])

h=legend("Theory","3D Simulation","2D Simulation",'FontSize',14)
set(h, 'Box', 'off')

set(gca, 'XScale', 'log', 'YScale', 'log');
set(gca, 'XTickLabel', []);
%set(gca, 'YTickLabel', []);
ax=gca;
ax.XTick = []; 
%ax.YTick = []; 




str="$$\mathrm{log}\ \sigma(t)$$";
annotation('textbox',[0.13,0.72,0.1,0.1],"String",str,"Interpreter","latex",'FontSize',18, 'EdgeColor', 'none')

str="$$\sim t^{-1/2}$$";
annotation('textbox',[0.5,0.6,0.1,0.1],"String",str,"Interpreter","latex",'FontSize',18, 'EdgeColor', 'none')


%annotation('textbox',[0.11,0.07,0.1,0.1],"String","High Frequency",'FontSize',18, 'EdgeColor', 'none',"FontWeight","normal")
%annotation('textbox',[0.425,0.07,0.1,0.1],"String","Shear Thinning",'FontSize',18, 'EdgeColor', 'none',"FontWeight","normal")
%annotation('textbox',[0.755,0.07,0.1,0.1],"String","Quasistatic",'FontSize',18, 'EdgeColor', 'none',"FontWeight","normal")


str="$$t<1$$";
annotation('textbox',[0.1825,0.0,0.1,0.1],"String",str,"Interpreter","latex",'FontSize',18, 'EdgeColor', 'none')
str="$$1<t<\sigma_{\infty}^{-2}$$";
annotation('textbox',[0.455,0.0,0.1,0.1],"String",str,"Interpreter","latex",'FontSize',18, 'EdgeColor', 'none')
str="$$t>\sigma_{\infty}^{-2}$$";
annotation('textbox',[0.782,0.0,0.1,0.1],"String",str,"Interpreter","latex",'FontSize',18, 'EdgeColor', 'none')