%% 读取图像并设置保存位置
clear;
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root, 'results', 'motion');
if ~exist(out, 'dir')
    mkdir(out);
end
diary(fullfile(out, 'run.log'));
fprintf('开始：%s\nMATLAB：%s\n', ...
    char(datetime('now', 'TimeZone', 'Asia/Shanghai')), version);
set(groot, 'DefaultFigureVisible', 'off');
I = imread(fullfile(root, 'images', 'book.tif'));
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

%% 按匀速直线运动公式生成 OTF
T = 1;
ab_values = [0.1 0.1; 0.05 0; 0.02 0.06]; % 对角、水平、偏垂直运动
mse_values = zeros(size(ab_values, 1), 1);
for index = 1:size(ab_values, 1)
    a = ab_values(index, 1);
    b = ab_values(index, 2);
    S = a*U + b*V;
    H = T * ones(P, Q);
    nonzero = (S ~= 0);
    H(nonzero) = T * sin(pi*S(nonzero)) ./ (pi*S(nonzero)) ...
        .* exp(-1i*pi*S(nonzero));
    % S=0 时按极限 H=T 处理，不对整幅 H 取实部或幅度来退化。
    assert(all(H(~nonzero) == T) && all(isfinite(H(:))));
    conjugate_error = max(abs(H(:)-reshape(conj(rot90(H, 2)), [], 1)));
    assert(conjugate_error < 1e-12);

    %% 频域相乘，反中心化后逆变换
    g_complex = ifft2(ifftshift(F .* H));
    imaginary_max = max(abs(imag(g_complex(:))));
    assert(imaginary_max < 1e-12);
    g = real(g_complex);
    J = g(1:M, 1:N);
    mse_values(index) = mean((J(:)-I(:)).^2);

    %% 显示复数 OTF 的幅度和退化结果
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
    title(sprintf('运动退化图像，MSE=%.5f', mse_values(index)));
    colormap(gray(256));
    sgtitle({sprintf('匀速直线运动：T=%g，a=%g，b=%g', T, a, b), ...
        sprintf('公式实现；补零 %d×%d；图像统一显示 [0,1]', P, Q)});
    drawnow;
    exportgraphics(f, fullfile(out, sprintf('motion_a_%g_b_%g.png', a, b)), 'Resolution', 180);
    close(f);
    save(fullfile(out, sprintf('motion_a_%g_b_%g.mat', a, b)), ...
        'I', 'H', 'J', 'T', 'a', 'b', 'P', 'Q', 'imaginary_max', 'conjugate_error');
    fprintf('T=%g，a=%g，b=%g，虚部最大值=%.3g，MSE=%.8g\n', ...
        T, a, b, imaginary_max, mse_values(index));
end
metrics = table(ab_values(:, 1), ab_values(:, 2), mse_values, ...
    'VariableNames', {'a', 'b', 'mse'});
writetable(metrics, fullfile(out, 'metrics.csv'));
fprintf('PASS：计算、检查及结果导出完成。\n');
diary off;
