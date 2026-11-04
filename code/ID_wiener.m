%% 读取图像并设置保存位置
clear;
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root, 'results', 'wiener');
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

%% 自编 Wiener 滤波，比较三组 K
K_values = [0.0001 0.001 0.01]; % 大写 K 为 Wiener 参数，与湍流 k 区分
for index = 1:numel(K_values)
    K = K_values(index);
    F_hat = conj(H) ./ (abs(H).^2 + K) .* G;
    restored_complex = ifft2(ifftshift(F_hat));
    J_full = real(restored_complex);
    J = J_full(1:M, 1:N);
    % 使用原始浮点结果计算，不能先截断或分别归一化。
    mse_value = mean((J(:)-I(:)).^2);
    psnr_value = 10*log10(1/mse_value);

    %% 每组同时展示理想图、同一退化图和复原结果
    f = figure('Visible', 'off', 'Color', 'w', 'Position', [80 80 1560 620], ...
        'DefaultAxesFontName', 'Microsoft YaHei', 'DefaultTextFontName', 'Microsoft YaHei', ...
        'DefaultAxesFontSize', 12, 'DefaultTextFontSize', 13);
    subplot(1, 3, 1);
    imshow(I, [0 1]);
    title('理想图像（退化前）');
    subplot(1, 3, 2);
    imshow(degraded, [0 1]);
    title(sprintf('同一退化图，MSE=%.5f', mse_degraded));
    subplot(1, 3, 3);
    imshow(J, [0 1]);
    title({sprintf('Wiener 复原，K=%g', K), ...
        sprintf('MSE=%.5f；PSNR=%.2f dB', mse_value, psnr_value)});
    colormap(gray(256));
    sgtitle({sprintf('Wiener 滤波：K=%g；湍流参数 k=0.00025', K), ...
        '高斯噪声方差 0.0001；三图统一显示 [0,1]；误差采用原始浮点值'});
    drawnow;
    exportgraphics(f, fullfile(out, sprintf('wiener_K_%g.png', K)), 'Resolution', 180);
    close(f);
end
fprintf('Wiener 滤波实验完成。\n');
