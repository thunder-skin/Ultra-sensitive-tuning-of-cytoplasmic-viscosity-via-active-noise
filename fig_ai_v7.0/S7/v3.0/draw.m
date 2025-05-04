clc;clear
color1=[[0 24 134]/255;[9 102 255]/255;[17 214 255]/255;[74 255 185]/255;[186 252 70]/255;[254 215 7]/255;[254 101 4]/255;[246 1 1]/255;[133 4 2]/255;[234 63 247]/255];
color2=[[37 43 128]/255;[60 109 180]/255;[72 198 235]/255;[129 201 152]/255;[189 214 56]/255;[251 205 17]/255;[239 94 33]/255;[235 29 34]/255;[222 115 230]/255;[125 19 21]/255];
color=color2;

%%画图格式
line_width_box=1;
line_width=3;
symbol="s";
marker_size=7;
line_width_2=1.5;

%%文字格式
font_name="Arial";
legend_font_size=16;
ax_font_size=16;
label_font_size=18;

fig=figure;
set(fig, 'Position', [50, 50, 600, 600]);
hold on

set(gca,"LineWidth",line_width_box,"FontSize",20,"FontName",font_name)

data=xlsread("data_1.xlsx")
a=plot(data(:,1),data(:,2),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",color(1,:),"Color",color(1,:),"LineWidth",line_width); 
b=plot(data(:,1),data(:,3),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",color(2,:),"Color",color(2,:),"LineWidth",line_width); 
c=plot(data(:,1),data(:,4),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",color(3,:),"Color",color(3,:),"LineWidth",line_width); 
d=plot(data(:,1),data(:,5),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",color(4,:),"Color",color(4,:),"LineWidth",line_width); 
e=plot(data(:,1),data(:,6),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",color(5,:),"Color",color(5,:),"LineWidth",line_width); 

nnn=3;

x_test=data(:,1);
y_test=log10(data(:,2));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(0.6,0.80,40);
y_fit = polyval(p, x_i);
f=plot(x_i,10.^y_fit,"--","LineWidth",2.5,"Color",color(1,:)); 

x_test=data(:,1);
y_test=log10(data(:,3));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(0.6,0.80,40);
y_fit = polyval(p, x_i);
f=plot(x_i,10.^y_fit,"--","LineWidth",2.5,"Color",color(2,:)); 

x_test=data(:,1);
y_test=log10(data(:,4));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(0.6,0.80,40);
y_fit = polyval(p, x_i);
f=plot(x_i,10.^y_fit,"--","LineWidth",2.5,"Color",color(3,:)); 

x_test=data(:,1);
y_test=log10(data(:,5));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(0.6,0.80,40);
y_fit = polyval(p, x_i);
f=plot(x_i,10.^y_fit,"--","LineWidth",2.5,"Color",color(4,:)); 

x_test=data(:,1);
y_test=log10(data(:,6));
p = polyfit(x_test, y_test, nnn);
x_i=linspace(0.6,0.80,40);
y_fit = polyval(p, x_i);
f=plot(x_i,10.^y_fit,"--","LineWidth",2.5,"Color",color(5,:)); 


ax=gca;
set(ax.XAxis, 'FontSize', ax_font_size);
set(ax.YAxis, 'FontSize', ax_font_size);



xlim([0.58,0.82])

xticks(0.6:0.05:0.8);

ylim([0.4,600])


set(gca, 'YScale', 'log');



legend=legend([a,b,c,d,e],"r=0.3","r=0.5","r=1.0","r=1.4","r=2.0","Interpreter","latex")
set(legend, 'box', 'off',"FontSize",20)

str="$$\mathrm{Volume\ Fraction}\ \ \phi$$";
xlabel(str,"Interpreter","latex",'FontSize',20)

str="$$\mathrm{Diffusion\ Coefficient}\ \ 10^{-5}\ D$$";
ylabel(str,"Interpreter","latex",'FontSize',20,"Rotation",90)

