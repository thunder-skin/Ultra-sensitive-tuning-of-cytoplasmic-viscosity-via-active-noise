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


plot(data(:,1),data(:,2),"x","Color",color(9,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(9,:));
plot(data(:,1),data(:,3),"x","Color",color(8,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(8,:));
plot(data(:,1),data(:,4),"x","Color",color(7,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(7,:));
plot(data(:,1),data(:,5),"x","Color",color(6,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(6,:));
plot(data(:,1),data(:,6),"x","Color",color(5,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(5,:));
plot(data(:,1),data(:,7),"x","Color",color(4,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(4,:));
plot(data(:,1),data(:,8),"x","Color",color(3,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(3,:));
plot(data(:,1),data(:,9),"x","Color",color(2,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(2,:));
plot(data(:,1),data(:,10),"x","Color",color(1,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(1,:));

x_test=log10(data(:,1));
y_test=log10(data(:,2));
xi=linspace(-2,7,100);
yi=csaps(x_test,y_test,0.9,xi);
a1=plot(10.^xi,10.^yi,"--","Color",color(9,:),"LineWidth",line_width);

x_test=log10(data(:,1));
y_test=log10(data(:,3));
xi=linspace(-2,7,100);
yi=csaps(x_test,y_test,0.9,xi);
a1=plot(10.^xi,10.^yi,"--","Color",color(8,:),"LineWidth",line_width);

x_test=log10(data(:,1));
y_test=log10(data(:,4));
xi=linspace(-2,7,100);
yi=csaps(x_test,y_test,0.9,xi);
a1=plot(10.^xi,10.^yi,"--","Color",color(7,:),"LineWidth",line_width);

x_test=log10(data(:,1));
y_test=log10(data(:,5));
xi=linspace(-2,7,100);
yi=csaps(x_test,y_test,0.9,xi);
a1=plot(10.^xi,10.^yi,"-","Color",color(6,:),"LineWidth",line_width);

x_test=log10(data(:,1));
y_test=log10(data(:,6));
xi=linspace(-2,7,100);
yi=csaps(x_test,y_test,0.9,xi);
a1=plot(10.^xi,10.^yi,"-","Color",color(5,:),"LineWidth",line_width);

x_test=log10(data(:,1));
y_test=log10(data(:,7));
xi=linspace(-2,7,100);
yi=csaps(x_test,y_test,0.9,xi);
a1=plot(10.^xi,10.^yi,"-","Color",color(4,:),"LineWidth",line_width);

x_test=log10(data(:,1));
y_test=log10(data(:,8));
xi=linspace(-2,7,100);
yi=csaps(x_test,y_test,0.9,xi);
a1=plot(10.^xi,10.^yi,"-","Color",color(3,:),"LineWidth",line_width);

x_test=log10(data(:,1));
y_test=log10(data(:,9));
xi=linspace(-2,7,100);
yi=csaps(x_test,y_test,0.9,xi);
a1=plot(10.^xi,10.^yi,"-","Color",color(2,:),"LineWidth",line_width);

x_test=log10(data(:,1));
y_test=log10(data(:,10));
xi=linspace(-2,7,100);
yi=csaps(x_test,y_test,0.9,xi);
a1=plot(10.^xi,10.^yi,"-","Color",color(1,:),"LineWidth",line_width);






xi=linspace(0,1e-5,1);
yi=linspace(0,1e-5,1);


b1=plot(xi,yi,"--","Color",color(9,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(9,:));
b2=plot(xi,yi,"--","Color",color(8,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(8,:));
b3=plot(xi,yi,"--","Color",color(7,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(7,:));
b4=plot(xi,yi,"-","Color",color(6,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(6,:));
b5=plot(xi,yi,"-","Color",color(5,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(5,:));
b6=plot(xi,yi,"-","Color",color(4,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(4,:));
b7=plot(xi,yi,"-","Color",color(3,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(3,:));
b8=plot(xi,yi,"-","Color",color(2,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(2,:));
b9=plot(xi,yi,"-","Color",color(1,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(1,:));

aa=4;
plot([1e3,10^5.5],[10.^(-1),10.^(-2.2)]*aa,"--","Color","k","LineWidth",line_width*1.5)

ax=gca;
set(ax.XAxis, 'FontSize', 14);
set(ax.YAxis, 'FontSize', 14);

ylim([1e-3,1.1])
xlim([0.01,7e6])

set(gca, 'xtick', [0.1 1 1e1 1e2 1e3 1e4 1e5]);

l = legend([b1, b2, b3, b4, b5, b6, b7, b8, b9], ...
    '$$2.5 \times 10^{-3}$$', ...
    '$$1.6 \times 10^{-3}$$', ...
    '$$1.0 \times 10^{-3}$$', ...
    '$$5.6 \times 10^{-4}$$', ...
    '$$3.2 \times 10^{-4}$$', ...
    '$$1.8 \times 10^{-4}$$', ...
    '$$1.0 \times 10^{-4}$$', ...
    '$$3.2 \times 10^{-5}$$', ...
    '$$1.0 \times 10^{-5}$$', ...
    'Interpreter', 'latex');set(l, 'Box', 'off')
set(gca, 'XScale', 'log', 'YScale', 'log');
str="$$\mathrm{Time}\ \ t$$";
xlabel(str,"Interpreter","latex",'FontSize',20)

str="$$\mathrm{Shear\ Stress}\ \ \sigma(t)/\sigma_0$$";
ylabel(str,"Interpreter","latex",'FontSize',20,"Rotation",90)

