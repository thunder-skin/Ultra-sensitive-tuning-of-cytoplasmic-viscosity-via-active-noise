clc;clear
data=xlsread("data.xlsx")
color1=[[0 24 134]/255;[9 102 255]/255;[17 214 255]/255;[74 255 185]/255;[186 252 70]/255;[254 215 7]/255;[254 101 4]/255;[246 1 1]/255;[133 4 2]/255;[234 63 247]/255];
color2=[[37 43 128]/255;[60 109 180]/255;[72 198 235]/255;[129 201 152]/255;[189 214 56]/255;[251 205 17]/255;[239 94 33]/255;[235 29 34]/255;[222 115 230]/255;[125 19 21]/255];
color3=[[1 0 0];[0 0 1];[0 0 0]]
color=color3;

%%画图格式
line_width_box=1;
line_width=3;
symbol="s";
marker_size=7;
line_width_2=1.5;

%%文字格式
font_name="Arial";
legend_font_size=14;
ax_font_size=16;
label_font_size=18;

fig=figure;
set(fig, 'Position', [50, 50, 600, 600]);
hold on

set(gca,"LineWidth",line_width_box,"FontSize",16,"FontName",font_name)

yyaxis left;

a=errorbar(data(:,1),data(:,2), data(:,3),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",color(1,:),"Color",color(1,:),"LineWidth",line_width,"DisplayName","Relative Viscosity"); 
b=errorbar(data(:,1),data(:,4), data(:,5),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",color(2,:),"Color",color(2,:),"LineWidth",line_width,"DisplayName","Controling Ability"); 

yyaxis right;

c=plot(data(:,1),data(:,6),"o","Marker","^","MarkerSize",5,"MarkerFaceColor",color(3,:),"Color",color(3,:),"LineWidth",line_width,"DisplayName","Relative Viscosity")

yyaxis left;
plot(data(:,7),data(:,8),"--","LineWidth",line_width_2,"Color",color(1,:))
plot(data(:,9),data(:,10),"--","LineWidth",line_width_2,"Color",color(2,:))


ax=gca;
set(ax.XAxis, 'FontSize', ax_font_size);
yyaxis left;
set(ax.YAxis, 'FontSize', ax_font_size);
ax.YColor="k"
yyaxis right;
set(ax.YAxis, 'FontSize', ax_font_size);
ax.YColor="k"

xlim([-0.1,0.9])
xticks(0:0.2:0.8);

yyaxis left;
ylim([0.003,8])
yticks([0.01, 0.1, 1, 5]); 
yticklabels({'0.01', '0.1', '1', '5'});

yyaxis right;
ylim([1,55]);
yticks([1,2, 5, 10, 20, 50]); 
yticklabels({'1','2', '5', '10', '20', '50'});


yyaxis left;
set(gca, 'YScale', 'log');
yyaxis right;
set(gca, 'YScale', 'log');

legen=legend([a,b,c],"Active Diffusion Coefficient","Passive Diffusion Coefficient","Fluctuation Enhancement")
set(legen, 'box', 'off')

yyaxis left;
str="$$\mathrm{Volume\ Fraction}\ \ \phi$$";
xlabel(str,"Interpreter","latex",'FontSize',20)

str="$$\mathrm{Diffusion\ Coefficient}\ \ D/D_0$$";
ylabel(str,"Interpreter","latex",'FontSize',20,"Rotation",90)

yyaxis right;

str="$$\mathrm{Fluctuation\ Enhancement}\ \ D_{\mathrm{active}}/D_{\mathrm{passive}}$$";
ylabel(str,"Interpreter","latex",'FontSize',20,"Rotation",90)