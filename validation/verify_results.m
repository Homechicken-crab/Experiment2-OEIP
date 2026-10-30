%% 从已保存数据核查模型、复原公式和导出结果
clear;
root = fileparts(fileparts(mfilename('fullpath')));
out = fullfile(root, 'results');
diary(fullfile(out, 'logs', 'verification.log'));
fprintf('MATLAB %s；核查时间 %s\n', version, ...
    char(datetime('now', 'TimeZone', 'Asia/Shanghai')));
inverse = load(fullfile(out, 'inverse', 'inverse_data.mat'));
pseudo = load(fullfile(out, 'pseudoinverse', 'pseudoinverse_data.mat'));
wiener = load(fullfile(out, 'wiener', 'wiener_data.mat'));

%% 三个独立进程必须生成逐元素相同的退化数据
names = {'same_ideal_image'; 'same_H'; 'same_noise'; 'same_g'; 'same_G'};
fields = {'I', 'H', 'noise', 'g', 'G'};
errors = zeros(numel(fields), 1);
tolerances = zeros(numel(fields), 1);
for index = 1:numel(fields)
    field = fields{index};
    assert(isequal(inverse.(field), pseudo.(field), wiener.(field)), ...
        '三种独立脚本的 %s 不一致。', field);
end

%% 独立重建坐标、噪声和完整退化模型
[M, N] = size(inverse.I);
P = inverse.P;
Q = inverse.Q;
[U, V] = meshgrid(-floor(Q/2):floor(Q/2), -floor(P/2):floor(P/2));
H = exp(-inverse.k*(U.^2+V.^2).^(5/6));
assert(isequal(H, inverse.H));
assert(H(floor(P/2)+1, floor(Q/2)+1) == 1);
rng(inverse.seed, 'twister');
noise = inverse.noise_mean + sqrt(inverse.noise_variance)*randn(P, Q);
assert(isequal(noise, inverse.noise));
expected_g = real(ifft2(fft2(inverse.I, P, Q).*ifftshift(H))) + noise;
names{end+1, 1} = 'degradation_full_grid';
errors(end+1, 1) = max(abs(expected_g(:)-inverse.g(:)));
tolerances(end+1, 1) = 1e-12;
names{end+1, 1} = 'noiseless_inverse';
errors(end+1, 1) = inverse.noiseless_error;
tolerances(end+1, 1) = 1e-8;

%% 直接核查全逆、幅度伪逆以及有符号差值
reference = real(ifft2(ifftshift(inverse.G./H)));
names{end+1, 1} = 'full_inverse_formula';
errors(end+1, 1) = max(abs(reference(:)-inverse.J_full(:)));
tolerances(end+1, 1) = 1e-12;
assert(isequal(inverse.J, inverse.J_full(1:M, 1:N)));
assert(isequal(inverse.difference, inverse.J-inverse.I));
assert(inverse.mse_restored == mean(inverse.difference(:).^2));
assert(inverse.outside_fraction > 0.5, '应如实保留全逆滤波的噪声放大。');

mask = abs(H) >= pseudo.threshold;
assert(isequal(mask, pseudo.mask));
assert(all(pseudo.inverse_filter(~mask) == 0));
assert(all(pseudo.inverse_filter(mask) == 1./H(mask)));
reference = real(ifft2(ifftshift(pseudo.G.*pseudo.inverse_filter)));
names{end+1, 1} = 'amplitude_pseudoinverse_formula';
errors(end+1, 1) = max(abs(reference(:)-pseudo.J_full(:)));
tolerances(end+1, 1) = 1e-12;
assert(isequal(pseudo.difference, pseudo.J-pseudo.I));
assert(pseudo.mse_restored == mean(pseudo.difference(:).^2));

%% 三个 Wiener K 均按原理计算，指标使用原始浮点结果
for index = 1:numel(wiener.K_values)
    K = wiener.K_values(index);
    reference = real(ifft2(ifftshift(conj(H)./(abs(H).^2+K).*wiener.G)));
    reference = reference(1:M, 1:N);
    J = wiener.J_all(:, :, index);
    names{end+1, 1} = sprintf('wiener_K_%g', K);
    errors(end+1, 1) = max(abs(reference(:)-J(:)));
    tolerances(end+1, 1) = 1e-12;
    assert(wiener.mse_values(index) == mean((J(:)-wiener.I(:)).^2));
end

%% 湍流三组及运动三组均重新检查公式和虚部
files = dir(fullfile(out, 'turbulence', '*.mat'));
assert(numel(files) == 3);
for index = 1:numel(files)
    item = load(fullfile(files(index).folder, files(index).name));
    [U, V] = meshgrid(-floor(item.Q/2):floor(item.Q/2), -floor(item.P/2):floor(item.P/2));
    expected_H = exp(-item.k*(U.^2+V.^2).^(5/6));
    assert(isequal(expected_H, item.H));
    reference = real(ifft2(fft2(item.I, item.P, item.Q).*ifftshift(expected_H)));
    reference = reference(1:size(item.I,1), 1:size(item.I,2));
    names{end+1, 1} = sprintf('turbulence_k_%g', item.k);
    errors(end+1, 1) = max(abs(reference(:)-item.J(:)));
    tolerances(end+1, 1) = 1e-12;
end
files = dir(fullfile(out, 'motion', '*.mat'));
assert(numel(files) == 3);
for index = 1:numel(files)
    item = load(fullfile(files(index).folder, files(index).name));
    [U, V] = meshgrid(-floor(item.Q/2):floor(item.Q/2), -floor(item.P/2):floor(item.P/2));
    S = item.a*U+item.b*V;
    expected_H = item.T*ones(item.P, item.Q);
    mask = S ~= 0;
    expected_H(mask) = item.T*sin(pi*S(mask))./(pi*S(mask)).*exp(-1i*pi*S(mask));
    assert(isequal(expected_H, item.H) && all(item.H(~mask) == item.T));
    assert(item.imaginary_max < 1e-12 && item.conjugate_error < 1e-12);
    reference = real(ifft2(fft2(item.I, item.P, item.Q).*ifftshift(expected_H)));
    reference = reference(1:size(item.I,1), 1:size(item.I,2));
    names{end+1, 1} = sprintf('motion_a_%g_b_%g', item.a, item.b);
    errors(end+1, 1) = max(abs(reference(:)-item.J(:)));
    tolerances(end+1, 1) = 1e-12;
end

%% 检查组合图确实保存且没有残留图窗
files = dir(fullfile(out, '**', '*.png'));
assert(numel(files) == 11, '应有 11 张组合结果图。');
for index = 1:numel(files)
    info = imfinfo(fullfile(files(index).folder, files(index).name));
    assert(info.Width >= 2000 && info.Height >= 700);
end
assert(isempty(findall(groot, 'Type', 'figure')));
passed = errors <= tolerances;
checks = table(names, errors, tolerances, passed);
writetable(checks, fullfile(out, 'verification_checks.csv'));
disp(checks);
assert(all(passed));
fprintf('PASS：%d 项数值核查，11 张组合图；三种复原的 I、H、noise、g、G 完全一致。\n', numel(passed));
diary off;
