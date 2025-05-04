clc;clear
data=xlsread("output_2.xlsx")
color1=[[0 24 134]/255;[9 102 255]/255;[17 214 255]/255;[74 255 185]/255;[186 252 70]/255;[254 215 7]/255;[254 101 4]/255;[246 1 1]/255;[133 4 2]/255;[234 63 247]/255];
color2=[[37 43 128]/255;[60 109 180]/255;[72 198 235]/255;[129 201 152]/255;[189 214 56]/255;[251 205 17]/255;[239 94 33]/255;[235 29 34]/255;[222 115 230]/255;[125 19 21]/255];
color=color2;

%%画图格式
line_width_box=1;
line_width=0.5;
symbol="x"
marker_size=4;

%%文字格式
font_name="Arial"
legend_font_size=14;
ax_font_size=10;
label_font_size=14;

fig=figure;
set(fig, 'Position', [50, 50, 800, 400]);
hold on
box on

set(gca,"LineWidth",line_width_box,"FontSize",legend_font_size,"FontName",font_name)
for i=1:7
    plot(data(:,1),data(:,i+1),"-","Color",color(i,:),"LineWidth",line_width,"Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",color(i,:))
end

xlim([0,1])
ylim([-1.1,1.1])

xticks([0,0.2,0.4,0.6,0.8,1]);
xticklabels(["0","0.2","0.4","0.6","0.8","1"]);


l=legend("$$1.0 \times 10^{-2}$$","$$3.3\times 10^{-3}$$","$$1.0 \times 10^{-3}$$","$$3.3 \times 10^{-4}$$","$$1.0 \times 10^{-4}$$","$$3.3 \times 10^{-5}$$","$$1.0 \times 10^{-5}$$","Location","eastoutside",'Interpreter', 'latex')
set(l, 'Box', 'off')
str="$$\mathrm{Time}\ \ t/T_{\mathrm{shear}}$$";
xlabel(str,"Interpreter","latex",'FontSize',20)

str="$$\mathrm{Shear\ Stress}\ \ \sigma(t)/\sigma_{\mathrm{max}}$$";
ylabel(str,"Interpreter","latex",'FontSize',20,"Rotation",90)