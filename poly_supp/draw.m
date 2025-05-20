% 读取粒子数据
function particles = read_particles(file_path)
    fid = fopen(file_path, 'r');
    if fid == -1
        error('File not found');
    end
    particles = [];
    while ~feof(fid)
        line = fgetl(fid);
        data = str2num(line);
        if length(data) == 5
            particles = [particles; data];
        end
    end
    fclose(fid);
end

% 绘制球体
function draw_sphere(ax, center, radius, color)
    [u, v] = meshgrid(linspace(0, 2*pi, 100), linspace(0, pi, 100));
    x = center(1) + radius * cos(u) .* sin(v);
    y = center(2) + radius * sin(u) .* sin(v);
    z = center(3) + radius * cos(v);
    surf(ax, x, y, z, 'FaceColor', color, 'EdgeColor', 'none');
end

% 主函数
function main()
    file_path = '.txt'; % 替换为你的txt文件路径
    particles = read_particles(file_path);

    % 指定立方空间范围
    min_x = 3;
    max_x = 7;
    min_y = 3;
    max_y = 7;
    min_z = 3;
    max_z = 7;

    % 创建3D图形
    fig = figure;
    ax = axes('NextPlot', 'add');

    % 定义颜色映射
    cmap = jet(256); % 使用jet颜色映射，生成256个颜色
    radius_min = min(particles(:, 4));
    radius_max = max(particles(:, 4));

    for i = 1:size(particles, 1)
        x = particles(i, 1);
        y = particles(i, 2);
        z = particles(i, 3);
        radius = particles(i, 4);
        id = particles(i, 5);

        % 检查粒子是否在指定立方空间内
        if x >= min_x && x <= max_x && y >= min_y && y <= max_y && z >= min_z && z <= max_z
            % 根据半径选择颜色
            color_index = round((radius - radius_min) / (radius_max - radius_min) * (size(cmap, 1) - 1)) + 1;
            color = cmap(color_index, :);
            draw_sphere(ax, [x, y, z], radius, color);
        end
    end

    % 设置坐标轴范围
    xlim(ax, [min_x, max_x]);
    ylim(ax, [min_y, max_y]);
    zlim(ax, [min_z, max_z]);

    % 添加坐标轴标签
    xlabel(ax, 'X');
    ylabel(ax, 'Y');
    zlabel(ax, 'Z');

    % 添加颜色条
    colorbar; % 直接调用 colorbar 函数
    colormap(fig, cmap); % 设置颜色映射
    caxis(fig, [radius_min, radius_max]); % 设置颜色条的范围
    ylabel(colorbar, 'Radius'); % 设置颜色条的标签

    % 设置视图
    view(ax, 3);
end

% 调用主函数
main();