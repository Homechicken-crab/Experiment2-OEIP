%% 读取图像并设置保存位置
clear;
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root, 'results', 'turbulence');
if ~exist(out, 'dir')
    mkdir(out);
end
set(groot, 'DefaultFigureVisible', 'off');
I = imread(fullfile(root, 'images', 'aerial_view_no_turb.tif'));
if size(I, 3) == 3
    I = rgb2gray(I);
end
I = im2double(I); % 转为 [0,1] 双精度
[M, N] = size(I);

%% 设置补零尺寸与中心化频率坐标
% 原图放在左上角；完整网格参与运算，最后取 1:M、1:N。
% 奇数网格的带符号频率坐标关于零频对称。
P = 2*M + 1;
Q = 2*N + 1;
u = -floor(Q/2):floor(Q/2); % 水平频率索引，对应列
v = -floor(P/2):floor(P/2); % 垂直频率索引，对应行
[U, V] = meshgrid(u, v);
F = fftshift(fft2(I, P, Q));

%% 按指导书公式生成大气湍流 OTF 并退化
k_values = [0.0025 0.001 0.00025]; % k 越大，退化越严重
for index = 1:numel(k_values)
    k = k_values(index);
    H = exp(-k * (U.^2 + V.^2).^(5/6));
    g_complex = ifft2(ifftshift(F .* H));
    g = real(g_complex);
    J = g(1:M, 1:N);
    mse_value = mean((J(:)-I(:)).^2);

    %% 显示原图、模型和退化图
    f = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1560 620], ...
        'DefaultAxesFontName', 'Microsoft YaHei', 'DefaultTextFontName', 'Microsoft YaHei', ...
        'DefaultAxesFontSize', 12, 'DefaultTextFontSize', 13);
    subplot(1, 3, 1);
    imshow(I, [0 1]);
    title('原始图像，显示 [0,1]');
    subplot(1, 3, 2);
    imagesc(u, v, abs(H), [0 1]);
    axis image;
    colorbar;
    xlabel('水平频率索引 u');
    ylabel('垂直频率索引 v');
    title(sprintf('湍流 OTF 幅度，k=%g', k));
    subplot(1, 3, 3);
    imshow(J, [0 1]);
    title(sprintf('退化图像，MSE=%.5f', mse_value));
    colormap(gray(256));
    sgtitle({sprintf('大气湍流退化：k=%g', k), ...
        sprintf('补零 %d×%d；统一显示 [0,1]；未作单幅归一化', P, Q)});
    drawnow;
    exportgraphics(f, fullfile(out, sprintf('turbulence_k_%g.png', k)), 'Resolution', 180);
    close(f);
end
fprintf('大气湍流退化实验完成。\n');
