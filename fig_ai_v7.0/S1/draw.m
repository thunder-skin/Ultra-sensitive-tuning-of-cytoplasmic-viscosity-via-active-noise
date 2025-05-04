clc;clear
data=xlsread("data.xlsx")
color=[[37 43 128]/255;[60 109 180]/255;[72 198 235]/255;[129 201 152]/255;[251 205 17]/255;[239 94 33]/255;[235 29 34]/255;[222 115 230]/255;[125 19 21]/255];


%%画图格式
line_width_box=1;
line_width=2;
symbol="x"
marker_size=8;

%%文字格式
font_name="Arial"
legend_font_size=18;
ax_font_size=16;
label_font_size=18;



fig=figure
set(fig, 'Position', [50, 50, 800, 600]);
hold on
box on

set(gca,"LineWidth",line_width_box,"FontSize",legend_font_size,"FontName",font_name,'Box', 'off')

%jamming以下描点
plot(data(:,1),data(:,9),"--","Color",color(1,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(1,:));
plot(data(:,1),data(:,10),"--","Color",color(2,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(2,:));
plot(data(:,1),data(:,11),"--","Color",color(3,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(3,:));
plot(data(:,1),data(:,12),"--","Color",color(4,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(4,:));
plot(data(:,1),data(:,13),"--","Color",color(5,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(5,:));
plot(data(:,1),data(:,14),"--","Color",color(6,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(6,:));

%jamming以上描点
plot(data(:,1),data(:,2),"o","Color",color(1,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(1,:));
plot(data(:,1),data(:,3),"o","Color",color(2,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(2,:));
plot(data(:,1),data(:,4),"o","Color",color(3,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(3,:));
plot(data(:,1),data(:,5),"o","Color",color(4,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(4,:));
plot(data(:,1),data(:,6),"o","Color",color(5,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(5,:));
plot(data(:,1),data(:,7),"o","Color",color(6,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(6,:));
plot(data(:,1),data(:,8),"o","Color",color(7,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(7,:));

%jamming以上连线拟合
%第一段
x_test=log10(data(:,1))
y_test=log10(data(:,2))
xi=linspace(-1.5,1.5,40)
yi=csaps(x_test,y_test,0.9,xi)
a1=plot(10.^xi,10.^yi,"-","Color",color(1,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,2))
xi=linspace(1.5,3.5,10)
yi=csaps(x_test,y_test,0.6,xi)
plot(10.^xi,10.^yi,"-","Color",color(1,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,2))
xi=linspace(3.5,5,40)
yi=csaps(x_test,y_test,0.99,xi)
plot(10.^xi,10.^yi+0.0013,"-","Color",color(1,:),"LineWidth",line_width)

%第二段
x_test=log10(data(:,1))
y_test=log10(data(:,3))
xi=linspace(-1.5,1.5,40)
yi=csaps(x_test,y_test,0.9,xi)
a2=plot(10.^xi,10.^yi,"-","Color",color(2,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,3))
xi=linspace(1.5,3.5,10)
yi=csaps(x_test,y_test,0.6,xi)
plot(10.^xi,10.^yi,"-","Color",color(2,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,3))
xi=linspace(3.5,5,40)
yi=csaps(x_test,y_test,0.99,xi)
plot(10.^xi,10.^yi+0.001,"-","Color",color(2,:),"LineWidth",line_width)

%第三段
x_test=log10(data(:,1))
y_test=log10(data(:,4))
xi=linspace(-1.5,1.5,40)
yi=csaps(x_test,y_test,0.9,xi)
a3=plot(10.^xi,10.^yi,"-","Color",color(3,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,4))
xi=linspace(1.5,3.1,10)
yi=csaps(x_test,y_test,0.6,xi)
plot(10.^xi,10.^yi,"-","Color",color(3,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,4))
xi=linspace(3.1,5,40)
yi=csaps(x_test,y_test,0.99,xi)
plot(10.^xi,10.^yi+0.001,"-","Color",color(3,:),"LineWidth",line_width)

%第四段
x_test=log10(data(:,1))
y_test=log10(data(:,5))
xi=linspace(-1.5,1.5,40)
yi=csaps(x_test,y_test,0.9,xi)
a4=plot(10.^xi,10.^yi,"-","Color",color(4,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,5))
xi=linspace(1.5,3.3,10)
yi=csaps(x_test,y_test,0.6,xi)
plot(10.^xi,10.^yi,"-","Color",color(4,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,5))
xi=linspace(3.3,5,40)
yi=csaps(x_test,y_test,0.99,xi)
plot(10.^xi,10.^yi*1.03,"-","Color",color(4,:),"LineWidth",line_width)

%第五段
x_test=log10(data(:,1))
y_test=log10(data(:,6))
xi=linspace(-1.5,1.5,40)
yi=csaps(x_test,y_test,0.9,xi)
a5=plot(10.^xi,10.^yi,"-","Color",color(5,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,6))
xi=linspace(1.5,3.3,10)
yi=csaps(x_test,y_test,0.6,xi)
plot(10.^xi,10.^yi,"-","Color",color(5,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,6))
xi=linspace(3.3,5,40)
yi=csaps(x_test,y_test,0.99,xi)
plot(10.^xi,10.^yi*1.06-0.001,"-","Color",color(5,:),"LineWidth",line_width)

%第六段
x_test=log10(data(:,1))
y_test=log10(data(:,7))
xi=linspace(-1.5,1.5,40)
yi=csaps(x_test,y_test,0.9,xi)
a6=plot(10.^xi,10.^yi,"-","Color",color(6,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,7))
xi=linspace(1.5,3.3,10)
yi=csaps(x_test,y_test,0.6,xi)
plot(10.^xi,10.^yi,"-","Color",color(6,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,7))
xi=linspace(3.3,5,40)
yi=csaps(x_test,y_test,0.99,xi)
plot(10.^xi,(10.^yi)*1.07,"-","Color",color(6,:),"LineWidth",line_width)

%第七段
x_test=log10(data(:,1))
y_test=log10(data(:,8))
xi=linspace(-1.5,5,40)
yi=csaps(x_test,y_test,0.7,xi)
a7=plot(10.^xi,10.^yi,"-","Color",color(7,:),"LineWidth",line_width)

%jamming以下连线
x_test=log10(data(:,1))
y_test=log10(data(:,9))
xi=linspace(-1.5,5,40)
yi=csaps(x_test,y_test,0.99,xi)
a8=plot(10.^xi,10.^yi,"--","Color",color(1,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,10))
xi=linspace(-1.5,5,40)
yi=csaps(x_test,y_test,0.99,xi)
a9=plot(10.^xi,10.^yi,"--","Color",color(2,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,11))
xi=linspace(-1.5,5,40)
yi=csaps(x_test,y_test,0.99,xi)
a10=plot(10.^xi,10.^yi,"--","Color",color(3,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,12))
xi=linspace(-1.5,3.5,40)
yi=csaps(x_test,y_test,0.99,xi)
a11=plot(10.^xi,10.^yi,"--","Color",color(4,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,13))
xi=linspace(-1.5,3.5,40)
yi=csaps(x_test,y_test,0.99,xi)
a12=plot(10.^xi,10.^yi,"--","Color",color(5,:),"LineWidth",line_width)

x_test=log10(data(:,1))
y_test=log10(data(:,14))
xi=linspace(-1.5,3.5,40)
yi=csaps(x_test,y_test,0.99,xi)
a13=plot(10.^xi,10.^yi,"--","Color",color(6,:),"LineWidth",line_width)

aa=0.3

plot([1e3,1e6],[10.^(-1),10.^(-2.6)]*aa,"--","Color","k","LineWidth",line_width*1.5)

xi=linspace(0,1e-5,1)
yi=linspace(0,1e-5,1)

b1=plot(xi,yi,"-","Color",color(1,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(1,:));
b2=plot(xi,yi,"-","Color",color(2,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(2,:));
b3=plot(xi,yi,"-","Color",color(3,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(3,:));
b4=plot(xi,yi,"-","Color",color(4,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(4,:));
b5=plot(xi,yi,"-","Color",color(5,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(5,:));
b6=plot(xi,yi,"-","Color",color(6,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(6,:));
b7=plot(xi,yi,"-","Color",color(7,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(7,:));

b8=plot(xi,yi,"--","Color",color(1,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(1,:));
b9=plot(xi,yi,"--","Color",color(2,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(2,:));
b10=plot(xi,yi,"--","Color",color(3,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(3,:));
b11=plot(xi,yi,"--","Color",color(4,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(4,:));
b12=plot(xi,yi,"--","Color",color(5,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(5,:));
b13=plot(xi,yi,"--","Color",color(6,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(6,:));


ax=gca;
set(ax.XAxis, 'FontSize', 14);
set(ax.YAxis, 'FontSize', 14);

ylim([3e-3,1])
xlim([3e-2,1e5])

set(gca, 'xtick', [0.1 1 1e1 1e2 1e3 1e4 1e5]);

l=legend([b13,b12,b11,b10,b9,b8,b7,b6,b5,b4,b3,b2,b1],"$$0.620$$","$$0.625$$","$$0.630$$","$$0.635$$","$$0.640$$","$$0.645$$","$$0.650$$","$$0.655$$","$$0.660$$","$$0.670$$","$$0.680$$","$$0.690$$","$$0.700$$",'Interpreter', 'latex');
set(l, 'Box', 'off')
set(gca, 'XScale', 'log', 'YScale', 'log');
str="$$\mathrm{Time}\ \ t$$";
xlabel(str,"Interpreter","latex",'FontSize',20)

str="$$\mathrm{Shear\ Stress}\ \ \sigma(t)/\sigma_0$$";
ylabel(str,"Interpreter","latex",'FontSize',20,"Rotation",90)

