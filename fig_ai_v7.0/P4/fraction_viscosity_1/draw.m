clc;clear
data=xlsread("visc.xlsx")
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
legend_font_size=18;
ax_font_size=10;
label_font_size=18;

fig=figure;
set(fig, 'Position', [50, 50, 600, 500]);
hold on

set(gca,"LineWidth",line_width_box,"FontWeight","normal","FontSize",legend_font_size,"FontName",font_name)

yyaxis left;
a=errorbar(data(:,1),data(:,3), data(:,4),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",'r',"Color",'r',"LineWidth",line_width,"DisplayName","Relative Viscosity"); 

yyaxis right;
b=errorbar(data(:,5),data(:,8), data(:,11),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor",'b',"Color",'b',"LineWidth",line_width,"DisplayName","Controling Ability"); 

data=xlsread("line.xlsx")
yyaxis left;
plot(data(:,1),data(:,2),"--","Color","r","LineWidth",line_width)

yyaxis right;
%plot(data(:,3),data(:,4),"--","Color","b","LineWidth",line_width)

plot([0.000135 0.000135],[0 10],"--","Color",'k',"LineWidth",line_width_2)


ax=gca;
set(ax.XAxis, 'FontSize', ax_font_size);
yyaxis left;
set(ax.YAxis, 'FontSize', ax_font_size);
ax.YColor="k"
yyaxis right;
set(ax.YAxis, 'FontSize', ax_font_size);
ax.YColor="k"


xlim([10.^(-7.25),10.^(-3)])


yyaxis left;
ylim([4e6,7e7])
yticks([1e7,2e7,3e7,5e7]); 
yticklabels({'1','2','3','5'});

yyaxis left;
%set(gca, 'YScale', 'log',"XScale","log");

l=legend([a,b],'Viscosity','Sensitivity');
set(l, 'box', 'off')

str="$$\mathrm{Active\ Force}\ \ kT_{\mathrm{eff}}/E_0$$";
xlabel(str,"Interpreter","latex",'FontSize',18)

str="$$\mathrm{Viscosity}\ \ 10^7\ \eta/\eta_0$$";
ylabel(str,"Interpreter","latex",'FontSize',18)

yyaxis right;
ylim([0,8])
str="$$\mathrm{Sensitivity}\ S$$";
ylabel(str,"Interpreter","latex",'FontSize',18)

%annotation("textbox","string",{'Glass Transition', 'Temperature'}, "FontSize",14,"FontName",font_name,'FontWeight', 'normal','EdgeColor','none',"Color","k",'HorizontalAlignment', 'center')
