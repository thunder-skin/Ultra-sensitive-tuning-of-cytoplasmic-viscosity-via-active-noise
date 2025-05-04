clear, clc

data = readmatrix("data_5.xlsx", NumHeaderLines = 1);
T = readcell("data_5.xlsx", TextType = "string");
T = regexprep(string(T(1, 2 : 2 : end)), "(\d+)e((+|-)?\d+\.\d+)", "$1*10^($2)");
T = arrayfun(@str2num, T);
a = num2cell(reshape(data, size(data, 1), 2, ''), 1 : 2);
a = squeeze(cellfun(@(x)x(all(x == x, 2), :), a, UniformOutput = false));
n = 101;
phi = zeros(n, numel(a));
logz = zeros(n, numel(a));
logT = log(repmat(T, n, 1));
figure
for i = 1 : numel(T)
phit = a{i}(:, 1);
logzt = log(a{i}(:, 2));
phi(:, i) = linspace(min(phit), max(phit), n);
logz(:, i) = interp1(phit, logzt, phi(:, i), "makima");
end
b = interp1(cat(3, phi', logT', logz'), linspace(1, numel(a), n), "makima");
phi = b(:, :, 1);
T = exp(b(:, :, 2));
z = exp(b(:, :, 3));
surf(phi, T, z, CData = log(z))
hl = [line(phi(1 : 20 : end, :)', T(1 : 20 : end, :)', z(1 : 20 : end, :)')
line(phi(:, 1 : 20 : end), T(:, 1 : 20 : end), z(:, 1 : 20 : end))];
set(hl, Color = "k", LineWidth = 1)
yscale log, zscale log, shading interp

axis square manual, grid minor

xlim([0.62 0.70])
ylim([10^-7 10^-3.25])
zlim([2e4 5e8])
zticks(10.^(4:9))
yticks(10.^(-7:-3))
xticks([0.62,0.66,0.70])

%line([.66 .66], [10^-3.25 10^-3.25], [5.8e6,10^4], Color = [0 0 0], LineStyle = "-.", LineWidth = 2)


line(data(:,1),[1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7 1e-7],data(:,2), Color = [0 0 0], LineStyle = "-", LineWidth = 2.5)
line(data(:,1),[10^-3.25 10^-3.25 10^-3.25 10^-3.25 10^-3.25 10^-3.25 10^-3.25 10^-3.25 10^-3.25 10^-3.25 10^-3.25 10^-3.25 10^-3.25 10^-3.25 10^-3.25 10^-3.25 10^-3.25],data(:,16)*1.01, Color = [0 0 0], LineStyle = "-", LineWidth = 2.5)

line([.66 .66], [1 1] * 10^-7, [6.5e7 1e4], Color = [0 0 0], LineStyle = "-.", LineWidth = 2)

Fz = scatteredInterpolant(phi(:), T(:), z(:));
T85 = interp1(phi(1, :)', T', 0.66);
Z85 = interp1(phi(1, :)', z', 0.66);
line(repmat(.66, size(T85)), T85, Z85, Color = [0 0 0], LineWidth = 2, LineStyle="-.")



x=[0.648 0.649 0.652 0.654 0.657 0.660 0.663 0.669]
y=[1e-7 1e-6 1e-5 10^-4.5 1e-4 10^-3.75 10^-3.5 10^-3.25]
z=[1700000 3300000 8000000 12000000 19000000 23000000 26000000 31000000]
line(x, y, z, Color = [0 0 0], LineStyle = "-", LineWidth = 2)

yticks([1e-6,1e-5,1e-4]); 

figure1 = gcf;
set(figure1, 'Position', [50, 50, 600, 600]);


% 创建 textbox
annotation(figure1,'textbox',...
[0.292617785367786 0.883112632336975 0.0456957340759088 0.0369020494621816],...
'String','$\phi_c(T)$',...
'Interpreter','latex',...
'FontSize',18,...
'FontName','Arial',...
'EdgeColor','none');

% 创建 textbox
annotation(figure1,'textbox',...
[0.189092144342145 0.814092131276089 0.0426282040488262 0.0358769924365036],...
'Color',[0 0 0],...
'String','Boundary',...
'Interpreter','none',...
'FontSize',16,...
'FontName','Arial',...
'EdgeColor','none');

% 创建 textbox
annotation(figure1,'textbox',...
[0.189092144342145 0.814092131276089 0.0426282040488262 0.0358769924365036],...
'Color',[0 0 1],...
'String','Fluid',...
'Interpreter','none',...
'FontSize',22,...
'FontName','Arial',...
'EdgeColor','none');

% 创建 textbox
annotation(figure1,'textbox',...
[0.18172034947035 0.746324477517545 0.0394230759392182 0.0358769924365036],...
'Color',[1 0 0],...
'String','Solid',...
'Interpreter','none',...
'FontSize',22,...
'FontName','Arial',...
'EdgeColor','none');

set(gca, fontsize=14)
zlabel("\eta", Rotation=0,FontSize=20)
ylabel("$$kT_{\mathrm{eff}}$$","Interpreter","latex",FontSize=18)
xlabel("\phi",FontSize=18)
set(gca, 'SortMethod', 'childorder');