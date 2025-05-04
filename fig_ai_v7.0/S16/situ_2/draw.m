clc;clear
data=xlsread("data.xlsx");
color=[[37 43 128]/255;[60 109 180]/255;[72 198 235]/255;[129 201 152]/255;[251 205 17]/255;[239 94 33]/255;[235 29 34]/255;[222 115 230]/255;[125 19 21]/255];


%%画图格式
line_width_box=1;
line_width=2.5;
symbol="x";
marker_size=12;

%%文字格式
font_name="Arial";
legend_font_size=18;
ax_font_size=16;
label_font_size=18;



fig=figure;
set(fig, 'Position', [50, 50, 800, 600]);
hold on
box on

set(gca,"LineWidth",line_width_box,"FontSize",legend_font_size,"FontName",font_name,'Box', 'off')

%jamming以下描点
plot(data(:,1),data(:,2),"o","Color",color(5,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(5,:));
plot(data(:,3),data(:,4),"o","Color",color(6,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(6,:));
plot(data(:,5),data(:,6),"o","Color",color(7,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(7,:));

%jamming以上描点
plot(data(:,7),data(:,8),"o","Color",color(8,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(8,:));
plot(data(:,9),data(:,10),"o","Color",color(1,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(1,:));
plot(data(:,11),data(:,12),"o","Color",color(2,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(2,:));
plot(data(:,13),data(:,14),"o","Color",color(3,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(3,:));
plot(data(:,15),data(:,16),"o","Color",color(4,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(4,:));
plot(data(:,17),data(:,18),"o","Color",color(5,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(5,:));
plot(data(:,19),data(:,20),"o","Color",color(6,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(6,:));
plot(data(:,21),data(:,22),"o","Color",color(7,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(7,:));

%jamming以下连线拟合
x_test=log10(data(:,1));
y_test=log10(data(:,2));
xi=linspace(-1.5,3.2,100);
yi=csaps(x_test,y_test,0.9,xi);
a1=plot(10.^xi,10.^yi,"--","Color",color(5,:),"LineWidth",line_width);

x_test=log10(data(:,3));
y_test=log10(data(:,4));
xi=linspace(-1.5,3.5,100);
yi=csaps(x_test,y_test,0.9,xi);
a1=plot(10.^xi,10.^yi,"--","Color",color(6,:),"LineWidth",line_width);

x_test=log10(data(:,5));
y_test=log10(data(:,6));
xi=linspace(-1.5,5.5,100);
yi=csaps(x_test,y_test,0.9,xi);
a1=plot(10.^xi,10.^yi,"--","Color",color(7,:),"LineWidth",line_width);

%jamming以上连线拟合

x_test=log10(data(:,7));
y_test=log10(data(:,8));
xi=linspace(-1.5,6,100);
yi=csaps(x_test,y_test,0.9,xi);
a1=plot(10.^xi,10.^yi,"-","Color",color(8,:),"LineWidth",line_width);

x_test=log10(data(:,9));
y_test=log10(data(:,10));
xi=linspace(-1.5,6,100);
yi=csaps(x_test,y_test,0.9,xi);
a1=plot(10.^xi,10.^yi,"-","Color",color(1,:),"LineWidth",line_width);

x_test=log10(data(:,11));
y_test=log10(data(:,12));
xi=linspace(-1.5,6,100);
yi=csaps(x_test,y_test,0.9,xi);
a1=plot(10.^xi,10.^yi,"-","Color",color(2,:),"LineWidth",line_width);

x_test=log10(data(:,13));
y_test=log10(data(:,14));
xi=linspace(-1.5,6,100);
yi=csaps(x_test,y_test,0.9,xi);
a1=plot(10.^xi,10.^yi,"-","Color",color(3,:),"LineWidth",line_width);

x_test=log10(data(:,15));
y_test=log10(data(:,16));
xi=linspace(-1.5,6,100);
yi=csaps(x_test,y_test,0.9,xi);
a1=plot(10.^xi,10.^yi,"-","Color",color(4,:),"LineWidth",line_width);

x_test=log10(data(:,17));
y_test=log10(data(:,18));
xi=linspace(-1.5,6,100);
yi=csaps(x_test,y_test,0.9,xi);
a1=plot(10.^xi,10.^yi,"-","Color",color(5,:),"LineWidth",line_width);

x_test=log10(data(:,19));
y_test=log10(data(:,20));
xi=linspace(-1.5,6,100);
yi=csaps(x_test,y_test,0.9,xi);
a1=plot(10.^xi,10.^yi,"-","Color",color(6,:),"LineWidth",line_width);

x_test=log10(data(:,21));
y_test=log10(data(:,22));
xi=linspace(-1.5,6,100);
yi=csaps(x_test,y_test,0.9,xi);
a1=plot(10.^xi,10.^yi,"-","Color",color(7,:),"LineWidth",line_width);


xi=linspace(0,1e-5,1);
yi=linspace(0,1e-5,1);

b8=plot(xi,yi,"--","Color",color(5,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(5,:));
b9=plot(xi,yi,"--","Color",color(6,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(6,:));
b10=plot(xi,yi,"--","Color",color(7,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(7,:));

b11=plot(xi,yi,"-","Color",color(8,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(8,:));
b1=plot(xi,yi,"-","Color",color(1,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(1,:));
b2=plot(xi,yi,"-","Color",color(2,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(2,:));
b3=plot(xi,yi,"-","Color",color(3,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(3,:));
b4=plot(xi,yi,"-","Color",color(4,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(4,:));
b5=plot(xi,yi,"-","Color",color(5,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(5,:));
b6=plot(xi,yi,"-","Color",color(6,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(6,:));
b7=plot(xi,yi,"-","Color",color(7,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(7,:));

aa=2;
plot([1e3,10^5.5],[10.^(-1),10.^(-2.2)]*aa,"--","Color","k","LineWidth",line_width*1.5)

ax=gca;
set(ax.XAxis, 'FontSize', 14);
set(ax.YAxis, 'FontSize', 14);

ylim([4e-3,1])
xlim([3e-2,1e6])

set(gca, 'xtick', [0.1 1 1e1 1e2 1e3 1e4 1e5]);

l=legend([b8,b9,b10,b11,b1,b2,b3,b4,b5,b6,b7],"$$0.640$$","$$0.650$$","$$0.655$$","$$0.660$$","$$0.670$$","$$0.680$$","$$0.690$$","$$0.700$$","$$0.710$$","$$0.720$$","$$0.730$$",'Interpreter', 'latex');
set(l, 'Box', 'off')
set(gca, 'XScale', 'log', 'YScale', 'log');
str="$$\mathrm{Time}\ \ t$$";
xlabel(str,"Interpreter","latex",'FontSize',20)

str="$$\mathrm{Shear\ Stress}\ \ \sigma(t)/\sigma_0$$";
ylabel(str,"Interpreter","latex",'FontSize',20,"Rotation",90)

