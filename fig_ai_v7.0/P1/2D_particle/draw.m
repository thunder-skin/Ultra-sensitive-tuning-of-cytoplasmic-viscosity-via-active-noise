clc;clear;
circles=xlsread("data.xlsx");
color=[[1 0 0];[0 0.7 0];[0 0 1];[0 0 0]];
p=[[0 0];[1 0];[1 1];[1 -1];[0 -1];[0 1];[-1 -1];[-1 0];[-1 1]];
color=[[255 165 0]]/255;

marker_size=1;

fig=figure;
set(fig, 'Position', [50, 50, 600, 600]);
hold on;

for j = 1:9
    xi=p(j,1);
    yi=p(j,2);
    for i = 1:size(circles, 1)
        x = circles(i, 1);
        y = circles(i, 2);
        r = circles(i, 3);
        a = circles(i, 5);
        b = circles(i, 6);
        theta = linspace(0, 2*pi, 100);
        cx = x + r * cos(theta)+xi*20;
        cy = y + r * sin(theta)+yi*20;
        plot(cx, cy, "k","MarkerSize",marker_size);  % 绘制圆形
        if j<2
            plot([x,x+a*1500], [y,y+b*1500], 'b', 'LineWidth', 2);
            %plot(x, y, "k","MarkerSize",2,"Marker","o","MarkerFaceColor","k");
            if a<0
                theta=atan(b/a)+pi;
            end
    
            if a>0
                theta=atan(b/a);
            end
    
            L=0.25;
    
            x=x+a*1500;
            y=y+b*1500;
    
            x1 = x + L*cos(theta);
            y1 = y + L*sin(theta);
    
            x2 = x + L*cos(theta + 2*pi/3);
            y2 = y + L*sin(theta + 2*pi/3);
    
            x3 = x + L*cos(theta - 2*pi/3);
            y3 = y + L*sin(theta - 2*pi/3);
            
            fill([x1, x2, x3], [y1, y2, y3], 'b', 'EdgeColor', 'none');
    
    
            sigma=0.00015;
            a = normrnd(0, sigma)*r*r;
            b = normrnd(0, sigma)*r*r;
            x = circles(i, 1);
            y = circles(i, 2);
            plot([x,x+a*1500], [y,y+b*1500], 'r', 'LineWidth', 2);
            %plot(x, y, "k","MarkerSize",2,"Marker","o","MarkerFaceColor","k");
            if a<0
                theta=atan(b/a)+pi;
            end    
            if a>0
                theta=atan(b/a);
            end    
            L=0.2;
            x=x+a*1500;
            y=y+b*1500;
            x1 = x + L*cos(theta);
            y1 = y + L*sin(theta);    
            x2 = x + L*cos(theta + 2*pi/3);
            y2 = y + L*sin(theta + 2*pi/3);    
            x3 = x + L*cos(theta - 2*pi/3);
            y3 = y + L*sin(theta - 2*pi/3);            
            fill([x1, x2, x3], [y1, y2, y3], 'r', 'EdgeColor', 'none');

            sigma=0.0003;
            a = normrnd(0, sigma)/r;
            b = normrnd(0, sigma)/r;
            x = circles(i, 1);
            y = circles(i, 2);
            plot([x,x+a*1500], [y,y+b*1500],"Color", color(1,:), 'LineWidth', 2);
            %plot(x, y, "k","MarkerSize",2,"Marker","o","MarkerFaceColor","k");
            if a<0
                theta=atan(b/a)+pi;
            end    
            if a>0
                theta=atan(b/a);
            end    
            L=0.2;
            x=x+a*1500;
            y=y+b*1500;
            x1 = x + L*cos(theta);
            y1 = y + L*sin(theta);    
            x2 = x + L*cos(theta + 2*pi/3);
            y2 = y + L*sin(theta + 2*pi/3);    
            x3 = x + L*cos(theta - 2*pi/3);
            y3 = y + L*sin(theta - 2*pi/3);            
            fill([x1, x2, x3], [y1, y2, y3], color(1,:), 'EdgeColor', 'none');
        end      
    end
end

plot([7,13], [20.5,20.5],"Color", "k", 'LineWidth', 6);
theta=0;

L=0.5;
x=13;
y=20.5;
x1 = x + L*cos(theta);
y1 = y + L*sin(theta);    
x2 = x + L*cos(theta + 2*pi/3);
y2 = y + L*sin(theta + 2*pi/3);    
x3 = x + L*cos(theta - 2*pi/3);
y3 = y + L*sin(theta - 2*pi/3);            
fill([x1, x2, x3], [y1, y2, y3], "k", 'EdgeColor', 'none');

plot([13,7], [-0.5,-0.5],"Color", "k", 'LineWidth', 6);
theta=-pi;

L=0.5;
x=7;
y=-0.5;
x1 = x + L*cos(theta);
y1 = y + L*sin(theta);    
x2 = x + L*cos(theta + 2*pi/3);
y2 = y + L*sin(theta + 2*pi/3);    
x3 = x + L*cos(theta - 2*pi/3);
y3 = y + L*sin(theta - 2*pi/3);            
fill([x1, x2, x3], [y1, y2, y3], "k", 'EdgeColor', 'none');

yy=0.06
annotation('textbox', [0.13, yy, 0.27, 0.05], 'String', "Elastic Force", 'EdgeColor', 'none', 'BackgroundColor', 'w','FontSize',14);
annotation('textbox', [0.38, yy, 0.30, 0.05], 'String', 'Thermal Fluctuation', 'EdgeColor', 'none', 'BackgroundColor', 'w','FontSize',14);
annotation('textbox', [0.71, yy, 0.20, 0.05], 'String', 'Active Force', 'EdgeColor', 'none', 'BackgroundColor', 'w','FontSize',14);

yyy=0.115
annotation('arrow', [0.185, 0.285], [yyy, yyy], 'LineWidth', 3, 'Color', "b");
annotation('arrow', [0.47, 0.57], [yyy, yyy], 'LineWidth', 3, 'Color', color(1,:));
annotation('arrow', [0.76, 0.86], [yyy, yyy], 'LineWidth', 3, 'Color', "r");

axis square
xlim([-1,21])
ylim([-1,21])

x = [0, 20, 20, 0, 0];
y = [0, 0, 20, 20, 0];

plot(x, y,'--k', 'LineWidth', 1.5);

axis off

