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
legend_font_size=20;
ax_font_size=16;
label_font_size=22;



fig=figure
set(fig, 'Position', [50, 50, 700, 600]);
hold on
box on

set(gca,"LineWidth",line_width_box,"FontSize",legend_font_size,"FontName",font_name,'Box', 'off')


%jamming以上描点
plot(data(:,1),data(:,2),"o","Color",color(1,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(1,:));
plot(data(:,1),data(:,4),"o","Color",color(2,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(2,:));
plot(data(:,1),data(:,5),"o","Color",color(3,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(3,:));
plot(data(:,1),data(:,6),"o","Color",color(4,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(4,:));
plot(data(:,1),data(:,7),"o","Color",color(5,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(5,:));


%jamming以下描点
plot(data(:,1),data(:,8),"o","Color",color(6,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(6,:));
plot(data(:,1),data(:,9),"o","Color",color(7,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(7,:));
plot(data(:,1),data(:,10),"o","Color",color(8,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(8,:));
plot(data(:,1),data(:,11),"o","Color",color(9,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(9,:));

%jamming以上连线拟合
%第一段
x_test=log10(data(:,1))
y_test=log10(data(:,2))
xi=linspace(-1.5,1.5,40)
yi=csaps(x_test,y_test,0.9,xi)
plot(10.^xi,10.^yi,"-","Color",color(1,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,2))
xi=linspace(1.5,4.5,40)
yi=csaps(x_test,y_test,0.99,xi)
plot(10.^xi,10.^yi,"-","Color",color(1,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,2))
xi=linspace(4.5,6,40)
yi=csaps(x_test,y_test,0.99,xi)
plot(10.^xi,10.^yi,"-","Color",color(1,:),"LineWidth",line_width)



%第二段
x_test=log10(data(:,1))
y_test=log10(data(:,4))
xi=linspace(-1.5,1.5,40)
yi=csaps(x_test,y_test,0.9,xi)
plot(10.^xi,10.^yi,"-","Color",color(2,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,4))
xi=linspace(1.5,4.5,40)
yi=csaps(x_test,y_test,0.99,xi)
plot(10.^xi,10.^yi,"-","Color",color(2,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,4))
xi=linspace(4.5,6,40)
yi=csaps(x_test,y_test,0.99,xi)
plot(10.^xi,10.^yi,"-","Color",color(2,:),"LineWidth",line_width)



%第三段
x_test=log10(data(:,1))
y_test=log10(data(:,5))
xi=linspace(-1.5,1.5,40)
yi=csaps(x_test,y_test,0.9,xi)
plot(10.^xi,10.^yi,"-","Color",color(3,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,5))
xi=linspace(1.5,4.5,40)
yi=csaps(x_test,y_test,0.99,xi)
plot(10.^xi,10.^yi,"-","Color",color(3,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,5))
xi=linspace(4.5,6,40)
yi=csaps(x_test,y_test,0.99,xi)
plot(10.^xi,10.^yi,"-","Color",color(3,:),"LineWidth",line_width)


%第四段
x_test=log10(data(:,1))
y_test=log10(data(:,6))
xi=linspace(-1.5,1.5,40)
yi=csaps(x_test,y_test,0.9,xi)
plot(10.^xi,10.^yi,"-","Color",color(4,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,6))
xi=linspace(1.5,4.5,40)
yi=csaps(x_test,y_test,0.99,xi)
plot(10.^xi,10.^yi,"-","Color",color(4,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,6))
xi=linspace(4.5,6,40)
yi=csaps(x_test,y_test,0.99,xi)
plot(10.^xi,10.^yi,"-","Color",color(4,:),"LineWidth",line_width)


%第五段
x_test=log10(data(:,1))
y_test=log10(data(:,7))
xi=linspace(-1.5,1.5,40)
yi=csaps(x_test,y_test,0.9,xi)
plot(10.^xi,10.^yi,"-","Color",color(5,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,7))
xi=linspace(1.5,4.5,40)
yi=csaps(x_test,y_test,0.99,xi)
plot(10.^xi,10.^yi,"-","Color",color(5,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,7))
xi=linspace(4.5,6,40)
yi=csaps(x_test,y_test,0.99,xi)
plot(10.^xi,10.^yi,"-","Color",color(5,:),"LineWidth",line_width)


x_test=log10(data(:,1))
y_test=log10(data(:,8))
xi=linspace(-1.5,6,40)
yi=csaps(x_test,y_test,0.8,xi)
plot(10.^xi,10.^yi,"--","Color",color(6,:),"LineWidth",line_width)


x_test=log10(data(:,1))
y_test=log10(data(:,9))
xi=linspace(-1.5,6,40)
yi=csaps(x_test,y_test,0.8,xi)
plot(10.^xi,10.^yi,"--","Color",color(7,:),"LineWidth",line_width)


x_test=log10(data(:,1))
y_test=log10(data(:,10))
xi=linspace(-1.5,6,40)
yi=csaps(x_test,y_test,0.8,xi)
plot(10.^xi,10.^yi,"--","Color",color(8,:),"LineWidth",line_width)


x_test=log10(data(:,1))
y_test=log10(data(:,11))
xi=linspace(-1.5,6,40)
yi=csaps(x_test,y_test,0.8,xi)
plot(10.^xi,10.^yi,"--","Color",color(9,:),"LineWidth",line_width)


aa=0.5
plot([1e3,1e6],[10.^(-1),10.^(-2.5)]*aa,"--","Color","k","LineWidth",line_width*1.5)

%绘制图标
xi=linspace(0,1e-5,1)
yi=linspace(0,1e-5,1)
b1=plot(xi,yi,"-","Color",color(1,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(1,:));
b3=plot(xi,yi,"-","Color",color(2,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(2,:));
b4=plot(xi,yi,"-","Color",color(3,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(3,:));
b5=plot(xi,yi,"-","Color",color(4,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(4,:));
b6=plot(xi,yi,"-","Color",color(5,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(5,:));

b7=plot(xi,yi,"--","Color",color(6,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(6,:));
b8=plot(xi,yi,"--","Color",color(7,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(7,:));
b9=plot(xi,yi,"--","Color",color(8,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(8,:));
b10=plot(xi,yi,"--","Color",color(9,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(9,:));

ax=gca;
set(ax.XAxis, 'FontSize', 14);
set(ax.YAxis, 'FontSize', 14);

ylim([1e-3,2e-1])
xlim([1e2,1e6])

l = legend([b10, b9, b8, b7, b6, b5, b4, b3, b1], ...
    '$$6.3 \times 10^{-4}$$', ...
    '$$4.0 \times 10^{-4}$$', ...
    '$$2.5 \times 10^{-4}$$', ...
    '$$1.6 \times 10^{-4}$$', ...
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

str="$$\mathrm{Shear\ Stress}\ \ \sigma(t)/\sigma_0$$";
ylabel(str,"Interpreter","latex",'FontSize',22,"Rotation",90)

