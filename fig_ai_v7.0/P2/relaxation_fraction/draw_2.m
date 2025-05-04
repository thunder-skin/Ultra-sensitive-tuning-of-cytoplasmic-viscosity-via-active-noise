clc;clear
data=xlsread("data.xlsx")
color=[[37 43 128]/255;[72 198 235]/255;[129 201 152]/255;[251 205 17]/255;[235 29 34]/255;[222 115 230]/255;[125 19 21]/255];
color = [[1 0 0];[0 0.8 0];[0 0 1];[1 0.5 0];[0.5 0 0.5];[0 1 1];[1 0 1];[0.8 0.8 0.1];[0.5 0.5 0.5]];
color2=[[37 43 128]/255;[60 109 180]/255;[72 198 235]/255;[129 201 152]/255;[189 214 56]/255;[251 205 17]/255;[239 94 33]/255;[235 29 34]/255;[222 115 230]/255;[125 19 21]/255];
color=color2;

%%画图格式
line_width_box=1;
line_width=3;
symbol="x"
marker_size=12;

%%文字格式
font_name="Arial"
legend_font_size=18;
ax_font_size=16;
label_font_size=20;



fig=figure
set(fig, 'Position', [50, 50, 700, 600]);
hold on
box on

set(gca,"LineWidth",line_width_box,"FontSize",legend_font_size,"FontName",font_name,'Box', 'off')

aa=1.2;
plot([1e3,10^5.5],[10.^(-1),10.^(-2.25)]*aa,"--","Color","k","LineWidth",line_width*1)

%jamming以上描点
plot(data(:,1),data(:,2),"o","Color",color(1,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(1,:));
plot(data(:,1),data(:,4),"o","Color",color(2,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(2,:));
plot(data(:,1),data(:,5),"o","Color",color(3,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(3,:));
plot(data(:,1),data(:,6),"o","Color",color(4,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(4,:));
plot(data(:,1),data(:,7),"o","Color",color(5,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(5,:));
plot(data(:,1),data(:,12),"o","Color",color(6,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(6,:));


%jamming以下描点
plot(data(:,1),data(:,8),"o","Color",color(7,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(7,:));
plot(data(:,1),data(:,9),"o","Color",color(8,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(8,:));
plot(data(:,1),data(:,10),"o","Color",color(9,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(9,:));
plot(data(:,1),data(:,11),"o","Color",color(10,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(10,:));






%绘制图标
xi=linspace(0,1e-5,1)
yi=linspace(0,1e-5,1)
b1=plot(xi,yi,"-","Color",color(1,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(1,:));
b3=plot(xi,yi,"-","Color",color(2,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(2,:));
b4=plot(xi,yi,"-","Color",color(3,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(3,:));
b5=plot(xi,yi,"-","Color",color(4,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(4,:));
b6=plot(xi,yi,"-","Color",color(5,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(5,:));
b7=plot(xi,yi,"-","Color",color(6,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(6,:));


b8=plot(xi,yi,"--","Color",color(7,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(7,:));
b9=plot(xi,yi,"--","Color",color(8,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(8,:));
b10=plot(xi,yi,"--","Color",color(9,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(9,:));
b11=plot(xi,yi,"--","Color",color(10,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(10,:));

nnn=3;



data=xlsread("data_2.xlsx")
x_test=log10(data(1:40,1));
y_test=log10(data(1:40,2));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(2,5.3,40);
y_fit = polyval(p, x_i);
plot(10.^(x_i),10.^(y_fit),"-","LineWidth",2.5,"Color",color(1,:))
plot([164204 1e6],[0.013383 0.013383],"-","LineWidth",2.5,"Color",color(1,:))

x_test=log10(data(1:40,1));
y_test=log10(data(1:40,4));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(2,5.4,40);
y_fit = polyval(p, x_i);
plot(10.^(x_i),10.^(y_fit),"-","LineWidth",2.5,"Color",color(2,:))
plot([251189 1e6],[0.0114745 0.0114745],"-","LineWidth",2.5,"Color",color(2,:))


x_test=log10(data(1:40,1));
y_test=log10(data(1:40,5));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(2,5.4,40);
y_fit = polyval(p, x_i);
plot(10.^(x_i),10.^(y_fit),"-","LineWidth",2.5,"Color",color(3,:))
plot([251189 1e6],[0.0102341 0.0102341],"-","LineWidth",2.5,"Color",color(3,:))

x_test=log10(data(1:40,1));
y_test=log10(data(1:40,6));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(2,5.5,40);
y_fit = polyval(p, x_i);
plot(10.^(x_i),10.^(y_fit),"-","LineWidth",2.5,"Color",color(4,:))
plot([316228 1e6],[0.0075 0.0073],"-","LineWidth",2.5,"Color",color(4,:))

x_test=log10(data(1:40,1));
y_test=log10(data(1:40,7));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(2,5.5,40);
y_fit = polyval(p, x_i);
plot(10.^(x_i),10.^(y_fit),"-","LineWidth",2.5,"Color",color(5,:))
plot([316228 1e6],[0.00547 0.0052],"-","LineWidth",2.5,"Color",color(5,:))

x_test=log10(data(1:39,1));
y_test=log10(data(1:39,12));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(2,6,40);
y_fit = polyval(p, x_i);
plot(10.^(x_i),10.^(y_fit),"-","LineWidth",2.5,"Color",color(6,:))



x_test=log10(data(1:38,1));
y_test=log10(data(1:38,8));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(2,6,40);
y_fit = polyval(p, x_i);
plot(10.^x_i,10.^y_fit,"--","LineWidth",2.5,"Color",color(7,:))

x_test=log10(data(1:36,1));
y_test=log10(data(1:36,9));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(2,6,40);
y_fit = polyval(p, x_i);
plot(10.^x_i,10.^y_fit,"--","LineWidth",2.5,"Color",color(8,:))

x_test=log10(data(1:35,1));
y_test=log10(data(1:35,10));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(2,6,40);
y_fit = polyval(p, x_i);
plot(10.^x_i,10.^y_fit,"--","LineWidth",2.5,"Color",color(9,:))

x_test=log10(data(1:33,1));
y_test=log10(data(1:33,11));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(2,6,40);
y_fit = polyval(p, x_i);
plot(10.^x_i,10.^y_fit,"--","LineWidth",2.5,"Color",color(10,:))






ax=gca;
set(ax.XAxis, 'FontSize', 14);
set(ax.YAxis, 'FontSize', 14);

ylim([1e-3,2e-1])
xlim([1e2,1e6])

l = legend([b11,b10, b9, b8, b7, b6, b5, b4, b3, b1], ...
    '$$6.3 \times 10^{-4}$$', ...
    '$$4.0 \times 10^{-4}$$', ...
    '$$2.5 \times 10^{-4}$$', ...
    '$$1.6 \times 10^{-4}$$', ...
    '$$1.3 \times 10^{-4}$$', ...
    '$$1.0 \times 10^{-4}$$', ...
    '$$4.7 \times 10^{-5}$$', ...
    '$$2.1 \times 10^{-5}$$', ...
    '$$1.0 \times 10^{-6}$$', ...
    'Athermal', ... % 使用 LaTeX 格式
    'Interpreter', 'latex');
set(l, 'Box', 'off')

set(gca, 'XScale', 'log', 'YScale', 'log');
str="$$\mathrm{Time}\ \ t/t_0$$";
xlabel(str,"Interpreter","latex",'FontSize',22)

annotation('textbox','String','$$T_{\rm eff,c}$$','Interpreter','latex','FontSize',22,'FontName','Arial','EdgeColor','none',"Color",color(6,:));

str="$$\mathrm{Shear\ Stress}\ \ \sigma(t)/\sigma_0$$";
ylabel(str,"Interpreter","latex",'FontSize',22,"Rotation",90)

