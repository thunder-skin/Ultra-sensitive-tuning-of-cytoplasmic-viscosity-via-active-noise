clc;clear
data=xlsread("relax_normal.xlsx")
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

a1=plot(data(:,16),data(:,17),"--","Color",color(1,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(1,:));
b1=plot(data(:,16),data(:,18),"--","Color",color(2,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(2,:));
c1=plot(data(:,16),data(:,19),"--","Color",color(3,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(3,:));
d1=plot(data(:,16),data(:,20),"--","Color",color(4,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(4,:));
e1=plot(data(:,16),data(:,21),"--","Color",color(5,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(5,:));

f1=plot(data(:,23),data(:,24),"--","Color",color(6,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(6,:));
g1=plot([1e1,1e1],[1e-4,1e-5],"--","Color",color(7,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(7,:));

h1=plot([1e1,1e1],[1e-4,1e-5],"-","Color",color(1,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(1,:));
i1=plot([1e1,1e1],[1e-4,1e-5],"-","Color",color(2,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(2,:));

plot(data(:,11),data(:,14),"o","Color",color(7,:),"LineWidth",line_width,"Marker","+","MarkerSize",marker_size,"MarkerFaceColor",color(7,:))
x_test=log10(data(:,11))
y_test=log10(data(:,14))
xi=linspace(1,5.9,40)
yi=csaps(x_test,y_test,0.9,xi)
plot(10.^xi,10.^yi,"--","Color",color(7,:),"LineWidth",line_width)

plot(data(:,11),data(:,13),"o","Color",color(1,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(1,:))
x_test=log10(data(:,11))
y_test=log10(data(:,13))
xi=linspace(1,5.9,40)
yi=csaps(x_test,y_test,0.9,xi)
plot(10.^xi,10.^yi,"--","Color",color(1,:),"LineWidth",line_width)

plot(data(:,11),data(:,12),"o","Color",color(2,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(2,:))
x_test=log10(data(:,11))
y_test=log10(data(:,12))
xi=linspace(1,5.9,40)
yi=csaps(x_test,y_test,0.9,xi)
plot(10.^xi,10.^yi,"--","Color",color(2,:),"LineWidth",line_width)

j1=plot(data(:,1),data(:,3),"-","Color",color(3,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(3,:));
k1=plot(data(:,1),data(:,4),"-","Color",color(4,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(4,:));
l1=plot(data(:,1),data(:,5),"-","Color",color(5,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(5,:));
m1=plot(data(:,1),data(:,6),"-","Color",color(6,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(6,:));
n1=plot(data(:,1),data(:,7),"-","Color",color(7,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(7,:));
o1=plot(data(:,1),data(:,8),"-","Color",color(8,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(8,:));
p1=plot(data(:,1),data(:,9),"-","Color",color(9,:),"LineWidth",line_width,"Marker","x","MarkerSize",marker_size,"MarkerFaceColor",color(9,:));

aa=0.27;
plot([1e2,1e6],[10.^(-0.5),10.^(-2.2)]*aa,"--","Color","k","LineWidth",line_width*1.5)


ax=gca;
set(ax.XAxis, 'FontSize', 14);
set(ax.YAxis, 'FontSize', 14);

ylim([3e-3,1])
xlim([1e-2,1e6])

l=legend([a1,b1,c1,d1,e1,f1,g1,h1,i1,j1,k1,l1,m1,n1,o1,p1],"0.820","0.825","0.830","0.835","0.840","0.841","0.842","0.843","0.844","0.845","0.850","0.860","0.870","0.880","0.890","0.900");
set(l, 'Box', 'off')
set(gca, 'XScale', 'log', 'YScale', 'log');
str="$$\mathrm{Time}\ \ t$$";
xlabel(str,"Interpreter","latex",'FontSize',20)

str="$$\mathrm{Shear\ Stress}\ \ \sigma(t)/\sigma_0$$";
ylabel(str,"Interpreter","latex",'FontSize',20,"Rotation",90)

