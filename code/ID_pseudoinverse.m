%% 读取图像并设置保存位置
clear;
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root, 'results', 'pseudoinverse');
if ~exist(out, 'dir')
    mkdir(out);
end
diary(fullfile(out, 'run.log'));
fprintf('开始：%s\nMATLAB：%s\n', ...
    char(datetime('now', 'TimeZone', 'Asia/Shanghai')), version);
set(groot, 'DefaultFigureVisible', 'off');
I = imread(fullfile(root, 'images', 'aerial_view_no_turb.tif'));
if size(I, 3) == 3
    I = rgb2gray(I);
end
I = im2double(I); % 转为 [0,1] 双精度
[M, N] = size(I);

%% 设置补零尺寸与中心化频率坐标
% 原图放在左上角；完整网格参与运算，最后取 1:M、1:N。
% 使用奇数网格，运动 OTF 无单独的 Nyquist 行列，便于保持共轭对称。
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
degraded = g(1:M, 1:N); % 此裁剪仅用于显示与保存，不用于构造 G
mse_degraded = mean((degraded(:)-I(:)).^2);
assert(all(H(:) > 0) && all(isfinite(G(:))));
assert(H(floor(P/2)+1, floor(Q/2)+1) == 1);
fprintf('模型：大气湍流 k=%g；原图 %d×%d；补零 %d×%d；裁剪左上角。\n', k, M, N, P, Q);
fprintf('噪声：高斯，均值=%g，方差=%g，rng(%d, twister)。\n', ...
    noise_mean, noise_variance, seed);
fprintf('min(H)=%.8g，最大逆增益=%.8g，退化图 MSE=%.8g。\n', ...
    min(H(:)), max(1./H(:)), mse_degraded);

%% 伪逆滤波：严格采用指导书式（5）的幅度阈值
threshold = 0.08; % |H| >= threshold 时求逆，否则置零
mask = abs(H) >= threshold;
inverse_filter = zeros(P, Q);
inverse_filter(mask) = 1 ./ H(mask);
F_hat = inverse_filter .* G;
restored_complex = ifft2(ifftshift(F_hat));
J_full = real(restored_complex);
J = J_full(1:M, 1:N);
assert(all(inverse_filter(~mask) == 0));
fprintf('幅度阈值=%.4g，保留频率比例=%.4f%%；没有使用半径截断。\n', ...
    threshold, 100*mean(mask(:)));

%% 使用未归一化、未截断的浮点数据计算误差
assert(all(isfinite(J(:))));
difference = J-I; % 有符号差值：复原结果减去理想图像
mse_restored = mean(difference(:).^2);
psnr_db = 10*log10(1/mse_restored);
outside_fraction = mean(J(:)<0 | J(:)>1);
value_range = [min(J(:)) max(J(:))];
fprintf('复原值域=[%.8g, %.8g]；超出 [0,1] 比例=%.4f%%。\n', ...
    value_range, 100*outside_fraction);
fprintf('原始浮点 MSE=%.8g；PSNR=%.6f dB。\n', mse_restored, psnr_db);

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
sgtitle({sprintf('伪逆滤波（式5）：|H| >= %.2f 时求逆，其余置零', threshold), ...
    '湍流 k=0.00025；高斯噪声方差 0.0001；显示限幅不参与误差计算'});
drawnow;
exportgraphics(f, fullfile(out, 'pseudoinverse_comparison.png'), 'Resolution', 180);
close(f);
save(fullfile(out, 'pseudoinverse_data.mat'), 'I', 'H', 'g', 'noise', 'G', ...
    'J', 'J_full', 'difference', 'P', 'Q', 'k', 'seed', 'noise_mean', 'noise_variance', ...
    'threshold', 'mask', 'inverse_filter', 'mse_degraded', 'mse_restored', 'psnr_db', 'outside_fraction', 'value_range');
metrics = table(mse_degraded, mse_restored, psnr_db, outside_fraction);
writetable(metrics, fullfile(out, 'metrics.csv'));
fprintf('PASS：计算、检查及结果导出完成。\n');
diary off;
