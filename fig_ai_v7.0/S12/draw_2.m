clear, clc

data = readmatrix("data_2.xlsx", NumHeaderLines = 1);
T = readcell("data_2.xlsx", TextType = "string");
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

Fz = scatteredInterpolant(phi(:), T(:), z(:));
ylim([10 ^ -5.5 1e-3])
T85 = interp1(phi(1, :)', T', 0.85);
Z85 = interp1(phi(1, :)', z', 0.85);
line(repmat(.85, size(T85)), T85, Z85, Color = [0 0 0], LineWidth = 2, LineStyle="-.")

v = reshape([10 ^ -5.5, 0.842, 10 ^ -5, 0.844, 10 ^ -4.5, 0.846, ...
10 ^ -4, 0.851, 10 ^ -3.5, 0.853, 10 ^ -3, 0.855], 2, '');
vi = interp1(v', linspace(1, size(v, 2), n));
line(vi(:, 2), vi(:, 1), Fz(vi(:, 2), vi(:, 1)), ...
Color = "k", LineStyle = "--", LineWidth = 2)
zlim([1e4 1e9])
zticks(10.^(4:9))
yticks(10.^(-5:-3))
xticks([0.83,0.85,0.87])
line([.85 .85], [1 1] * 10^-5.5, [min(z(:)) Fz(.85, 10^-5.5)], ...
Color = [0 0 0], LineStyle = "-.", LineWidth = 2)


figure1 = gcf;
set(figure1, 'Position', [50, 50, 600, 500]);
% 创建 textbox
annotation(figure1,'textbox',...
[0.292617785367786 0.883112632336975 0.0456957340759088 0.0369020494621816],...
'String','Boundary',...
'FontSize',18,...
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
'String','Jam',...
'Interpreter','none',...
'FontSize',22,...
'FontName','Arial',...
'EdgeColor','none');
set(gca, fontsize=14)
zlabel("\eta", Rotation=0,FontSize=20)
ylabel("$$kT_{\mathrm{eff}}$$","Interpreter","latex",FontSize=18)
xlabel("\phi",FontSize=18)
set(gca, 'SortMethod', 'childorder');