clc;clear
data=xlsread("relax_with_temp.xlsx")
color1=[[0 24 134]/255;[9 102 255]/255;[17 214 255]/255;[74 255 185]/255;[186 252 70]/255;[254 215 7]/255;[254 101 4]/255;[246 1 1]/255;[133 4 2]/255;[234 63 247]/255];
color2=[[37 43 128]/255;[60 109 180]/255;[72 198 235]/255;[129 201 152]/255;[189 214 56]/255;[251 205 17]/255;[239 94 33]/255;[235 29 34]/255;[222 115 230]/255;[125 19 21]/255];
color=color2;

%%画图格式
line_width_box=1;
line_width=2.5;
symbol="x"
marker_size=12;

%%文字格式
font_name="Arial"
legend_font_size=16;
ax_font_size=14;
label_font_size=18;

fig=figure;
set(fig, 'Position', [50, 50, 600, 600]);



hold on
box on



set(gca,"LineWidth",line_width_box,"FontSize",legend_font_size,"FontName",font_name,'Box', 'off')
for i=1:2
    plot(data(:,1),data(:,i+1),"o","Color",color(i,:),"LineWidth",line_width,"Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",color(i,:))
end

for i=3:9
    plot(data(:,1),data(:,i+1),"o","Color",color(i,:),"LineWidth",line_width,"Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",color(i,:))
end

for i=1:2
    x_test=log10(data(:,1))
    y_test=log10(data(:,i+1))
    xi=linspace(2,5.4,30)
    yi=csaps(x_test,y_test,0.9,xi)
    plot(10.^xi,10.^yi,"-","Color",color(i,:),"LineWidth",line_width)
end

for i=3:9
    x_test=log10(data(:,1))
    y_test=log10(data(:,i+1))
    xi=linspace(2,6.5,40)
    yi=csaps(x_test,y_test,0.95,xi)
    plot(10.^xi,10.^yi,"-","Color",color(i,:),"LineWidth",line_width)
end

a=5.5;
plot([10.^(2.5),1e5]*a,[10.^(-1.25),10.^(-2.5)],"--","Color","k","LineWidth",line_width)

xlim([8e1,2e6])
ylim([2e-3,2.5e-1])

ax=gca;
set(ax.XAxis, 'FontSize', ax_font_size);
set(ax.YAxis, 'FontSize', ax_font_size);

str="$$\mathrm{Time}\ \ t$$";
xlabel(str,"Interpreter","latex",'FontSize',18)

str="$$\mathrm{Shear\ Stress}\ \ \sigma(t)/\sigma(0)$$";
ylabel(str,"Interpreter","latex",'FontSize',18,"Rotation",90)


h=legend("$$3.2\times 10^{-3}$$","$$1.8\times 10^{-3}$$","$$1.0 \times 10^{-3}$$","$$5.6 \times 10^{-4}$$","$$3.2 \times 10^{-4}$$","$$1.8 \times 10^{-4}$$","$$1.0 \times 10^{-4}$$","$$3.2 \times 10^{-5}$$","$$1.0 \times 10^{-5}$$","Location","eastoutside",'Interpreter', 'latex')
set(h, 'Box', 'off')
set(gca, 'XScale', 'log', 'YScale', 'log');

