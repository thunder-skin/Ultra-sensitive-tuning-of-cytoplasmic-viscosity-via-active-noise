clc;clear
data=xlsread("robust.xlsx")
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
ax_font_size=14;
label_font_size=18;

fig=figure;
set(fig, 'Position', [50, 50, 600, 500]);

hold on

set(gca,"LineWidth",line_width_box,"FontSize",legend_font_size,"FontName",font_name)

yyaxis right;
for i=0:0
    a=plot(data(:,1),data(:,2),"o","Color","r","LineWidth",line_width,"Marker","o","MarkerFaceColor","r")
end

yyaxis left;
for i=1:1
    b=errorbar(data(:,1),data(:,3), data(:,4),"o","Marker",symbol,"MarkerSize",marker_size,"MarkerFaceColor","b","Color","b","LineWidth",line_width); 
end



yyaxis right;
for i=0:0
    xi=linspace(-5.3,-1.7,30)
    yi=linspace(0.248,0.248,30)
    plot(10.^xi,yi,"--","Color","r","LineWidth",line_width_2)
end



ax=gca;
set(ax.XAxis, 'FontSize', ax_font_size);
yyaxis left;
set(ax.YAxis, 'FontSize', ax_font_size);
ax.YColor="k"
yyaxis right;
set(ax.YAxis, 'FontSize', ax_font_size);
ax.YColor="k"


xlim([10.^(-5.5),10.^(-1.5)])
xticks([1e-5,1e-4,1e-3,1e-2]); 
yyaxis left;
ylim([100,600])
yticks([100,200,300,400,500,600]); 
yticklabels({'100','200','300','400','500','600'});
yyaxis right;
ylim([0,0.5])
yticks(0:0.1:0.5);

ax=gca;
set(ax.XAxis, 'FontSize', ax_font_size);
set(ax.YAxis, 'FontSize', ax_font_size);

yyaxis left;
set(gca,"XScale","log");

l=legend([a,b],'Shear Frequency','Stress Amplitude',"Interpreter","latex");
set(l, 'box', 'off')

yyaxis left;
str="$$\mathrm{Frequency}\ f/t^{-1}_0\ $$";
xlabel(str,"Interpreter","latex",'FontSize',18)

str="$$\mathrm{Stress\ Amplitude}\ \sigma_{\rm max}$$";
ylabel(str,"Interpreter","latex",'FontSize',18,"Rotation",90)

yyaxis right;
str="$$\mathrm{Phase\ Shift}\ \Delta \theta/\pi $$";
ylabel(str,"Interpreter","latex",'FontSize',18,"Rotation",90)