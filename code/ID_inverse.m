%% 读取图像并设置保存位置
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

%% 构造三种复原共用的模糊与噪声模型
k = 0.00025; % 湍流参数，小写 k
seed = 20261030;
noise_mean = 0;
noise_variance = 0.0001; % 输入灰度以 [0,1] 为单位，标准差 0.01
H = exp(-k * (U.^2 + V.^2).^(5/6));
blurred_complex = ifft2(ifftshift(F .* H));
blurred = real(blurred_complex);
rng(seed, 'twister');
noise = noise_mean + sqrt(noise_variance)*randn(P, Q);
g = blurred + noise; % 对整个补零网格加噪，保留原始浮点数据，不截断
G = fftshift(fft2(g)); % 不能先裁剪 g 再进行复原
degraded = g(1:M, 1:N); % 此裁剪仅用于显示，不用于构造 G
mse_degraded = mean((degraded(:)-I(:)).^2);

%% 全逆滤波：所有频率直接相除，不加阈值、eps 或正则项
% 所选湍流 H 严格大于零，直接求逆，不修改 H。
F_hat = G ./ H;
restored_complex = ifft2(ifftshift(F_hat));
J_full = real(restored_complex);
J = J_full(1:M, 1:N);

%% 使用未归一化、未截断的浮点数据计算误差
difference = J-I; % 有符号差值：复原结果减去理想图像
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
