clc;clear
data=xlsread("hypo_proof.xlsx")
color1=[[0 24 134]/255;[9 102 255]/255;[17 214 255]/255;[74 255 185]/255;[186 252 70]/255;[254 215 7]/255;[254 101 4]/255;[246 1 1]/255;[133 4 2]/255;[234 63 247]/255];
color2=[[37 43 128]/255;[60 109 180]/255;[72 198 235]/255;[129 201 152]/255;[189 214 56]/255;[251 205 17]/255;[239 94 33]/255;[235 29 34]/255;[222 115 230]/255;[125 19 21]/255];
color=[[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0];[0 0 0]];

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

a2=plot(data(:,2),data(:,1),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor","r","Color","r","LineWidth",line_width); 

yi=linspace(-5.25,-3.05);
xi=yi/2.02;
aaa=0.7;
plot(10.^xi/1.17,(10.^yi)/aaa,"--","Color","r","LineWidth",line_width_2)

ylim([10.^(-5.2),10.^(-2.8)])
xlim([0.0016,0.03])

xticks([0.002,0.005,0.01,0.02]); 
xticklabels({'0.002','0.005','0.01','0.02'});

ax=gca;
set(ax.XAxis, 'FontSize', ax_font_size);
set(ax.YAxis, 'FontSize', ax_font_size);

set(gca, 'YScale', 'log',"XScale","log");

str="$$\mathrm{Volume\ Fraction}\ \ \phi-\phi_J$$";
xlabel(str,"Interpreter","latex",'FontSize',22)

str="$$\mathrm{Critical\ Effective\ Force}\ \ kT_{\mathrm{eff,c}}$$";
ylabel(str,"Interpreter","latex",'FontSize',22)

annotation('textbox','String','$kT_{\mathrm{eff,c}} \sim (\phi-\phi_J)^{2.02}$','Interpreter','latex','FontSize',22,'FontName','Arial','EdgeColor','none');

