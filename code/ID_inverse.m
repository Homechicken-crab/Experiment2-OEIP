%% 读取图像
clear;
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root, 'results', 'inverse');
if ~exist(out, 'dir')
    mkdir(out);
end
set(groot, 'DefaultFigureVisible', 'off');
I = imread(fullfile(root, 'images', 'aerial_view_no_turb.tif'));
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

%% 模糊与加噪
k = 0.00025;
seed = 20261030;
noise_mean = 0;
noise_variance = 0.0001; % 高斯噪声方差
H = exp(-k * (U.^2 + V.^2).^(5/6));
blurred_complex = ifft2(ifftshift(F .* H));
blurred = real(blurred_complex);
rng(seed, 'twister');
noise = noise_mean + sqrt(noise_variance)*randn(P, Q);
g = blurred + noise;
G = fftshift(fft2(g));
degraded = g(1:M, 1:N);
mse_degraded = mean((degraded(:)-I(:)).^2);

%% 全逆滤波
F_hat = G ./ H;
restored_complex = ifft2(ifftshift(F_hat));
J_full = real(restored_complex);
J = J_full(1:M, 1:N); % 复原后再裁剪

%% 计算误差
difference = J-I;
mse_restored = mean(difference(:).^2);
outside_fraction = mean(J(:)<0 | J(:)>1);
value_range = [min(J(:)) max(J(:))];

%% 显示并保存结果
f = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1560 620], ...
    'DefaultAxesFontName', 'Microsoft YaHei', 'DefaultTextFontName', 'Microsoft YaHei', ...
    'DefaultAxesFontSize', 12, 'DefaultTextFontSize', 13);
subplot(1, 3, 1);
imshow(degraded, [0 1]);
title(sprintf('同一退化图，MSE=%.5f', mse_degraded));
subplot(1, 3, 2);
imshow(J, [0 1]);
title({sprintf('复原结果，显示 [0,1]，MSE=%.4g', mse_restored), ...
    sprintf('原始值域 [%.3g, %.3g]；越界 %.2f%%', value_range, 100*outside_fraction)});
subplot(1, 3, 3);
limit = max(abs(difference(:)));
imagesc(difference, [-limit limit]);
axis image off;
colorbar;
title({'有符号差值：复原结果 - 理想图像', ...
    sprintf('原始差值色标 [-%.3g, %.3g]', limit, limit)});
colormap(gray(256));
sgtitle({'全逆滤波：Fhat=G./H（无阈值、无正则项）', ...
    '湍流 k=0.00025；高斯噪声方差 0.0001；显示限幅不参与误差计算'});
drawnow;
exportgraphics(f, fullfile(out, 'inverse_comparison.png'), 'Resolution', 180);
close(f);
fprintf('全逆滤波实验完成。\n');
