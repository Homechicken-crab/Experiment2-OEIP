%% 读取图像
clear;
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root, 'results', 'motion');
if ~exist(out, 'dir')
    mkdir(out);
end
set(groot, 'DefaultFigureVisible', 'off');
I = imread(fullfile(root, 'images', 'book.tif'));
if size(I, 3) == 3
    I = rgb2gray(I);
end
I = im2double(I);
[M, N] = size(I);

%% 设置频率坐标
P = 2*M + 1;
Q = 2*N + 1;
u = -floor(Q/2):floor(Q/2);
v = -floor(P/2):floor(P/2);
[U, V] = meshgrid(u, v);
F = fftshift(fft2(I, P, Q));

%% 运动退化模型
T = 1;
ab_values = [0.1 0.1; 0.05 0; 0.02 0.06];
for index = 1:size(ab_values, 1)
    a = ab_values(index, 1);
    b = ab_values(index, 2);
    S = a*U + b*V;
    H = T * ones(P, Q); % S=0 时 H=T
    nonzero = (S ~= 0);
    H(nonzero) = T * sin(pi*S(nonzero)) ./ (pi*S(nonzero)) ...
        .* exp(-1i*pi*S(nonzero));

    %% 运动模糊
    g_complex = ifft2(ifftshift(F .* H));
    g = real(g_complex);
    J = g(1:M, 1:N);
    mse_value = mean((J(:)-I(:)).^2);

    %% 显示并保存结果
    f = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1560 620], ...
        'DefaultAxesFontName', 'Microsoft YaHei', 'DefaultTextFontName', 'Microsoft YaHei', ...
        'DefaultAxesFontSize', 12, 'DefaultTextFontSize', 13);
    subplot(1, 3, 1);
    imshow(I, [0 1]);
    title('原始图像，显示 [0,1]');
    subplot(1, 3, 2);
    imagesc(u, v, abs(H), [0 T]);
    axis image;
    colorbar;
    xlabel('水平频率索引 u');
    ylabel('垂直频率索引 v');
    title(sprintf('运动 OTF 幅度，a=%g，b=%g', a, b));
    subplot(1, 3, 3);
    imshow(J, [0 1]);
    title(sprintf('运动退化图像，MSE=%.5f', mse_value));
    colormap(gray(256));
    sgtitle({sprintf('匀速直线运动：T=%g，a=%g，b=%g', T, a, b), ...
        sprintf('公式实现；补零 %d×%d；图像统一显示 [0,1]', P, Q)});
    drawnow;
    exportgraphics(f, fullfile(out, sprintf('motion_a_%g_b_%g.png', a, b)), 'Resolution', 180);
    close(f);
end
fprintf('运动模糊退化实验完成。\n');
