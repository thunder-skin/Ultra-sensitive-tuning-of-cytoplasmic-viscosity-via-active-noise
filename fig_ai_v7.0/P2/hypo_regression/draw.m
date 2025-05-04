clc;clear

%%画图格式
line_width_box=1;
line_width=3;
symbol="s";
marker_size=10;
line_width_2=2;

%%文字格式
font_name="Arial";
legend_font_size=20;
ax_font_size=14;
label_font_size=20;

fig=figure;
set(fig, 'Position', [50, 50, 500, 600]);

hold on

set(gca,"LineWidth",line_width_box,"FontSize",legend_font_size,"FontName",font_name)

%a1=errorbar(data(:,1),data(:,2), data(:,3),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor","b","Color","b","LineWidth",line_width); 

data=xlsread("3D.xlsx")
a2=plot(data(:,2),data(:,1),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor","r","Color","r","LineWidth",line_width); 

%xi=linspace(-5.6,-2.9,30);
%yi=xi/2.04;
%aaa=1.43;
%plot(10.^xi,(10.^yi)/aaa,"--","Color","b","LineWidth",line_width_2)

yi=linspace(-6.1,-3.15,30);
xi=yi/2.03;
aaa=1.05;
plot(10.^xi/1.17,(10.^yi)/aaa,"--","Color","r","LineWidth",line_width_2)

ylim([10.^(-6.25),10.^(-3.05)])
xlim([0.0006,0.03])

xticks([0.002,0.005,0.01,0.02]); 
xticklabels({'0.002','0.005','0.01','0.02'});


ax=gca;
set(ax.XAxis, 'FontSize', ax_font_size);
set(ax.YAxis, 'FontSize', ax_font_size);

set(gca, 'YScale', 'log',"XScale","log");


str="$$\mathrm{Volume\ Fraction}\ \ \phi-\phi_J$$";
xlabel(str,"Interpreter","latex",'FontSize',22)

str="$$\mathrm{Critical\ Effective\ Force}\ \ T_{\mathrm{eff,c}}$$";
ylabel(str,"Interpreter","latex",'FontSize',22)

annotation('textbox','String','$T_{\mathrm{eff,c}} \sim (\phi-\phi_J)^{2.03}$','Interpreter','latex','FontSize',22,'FontName','Arial','EdgeColor','none');

%annotation("textbox","string",'k=0.490', "FontSize",18,"FontName",font_name,'EdgeColor','none',"Color","b","Rotation",51)