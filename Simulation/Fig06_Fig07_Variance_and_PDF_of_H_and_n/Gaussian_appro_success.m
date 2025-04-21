close all
clear all
clc
tic

seed = 1;
rng(seed)

data_size = 1000000;
gca_size = 22;

%% H_FA的方差估计准确性
var_h_arr = 0.001:0.001:0.5;
var_h_true = zeros(1,length(var_h_arr));
var_h_est  = zeros(1,length(var_h_arr));
parfor var_h_idx = 1:length(var_h_arr)

    var_h = var_h_arr(var_h_idx);
    
    hab1 = gen_complex_h(var_h, data_size);
    hba1 = hab1;
    hab2 = gen_complex_h(var_h, data_size);
    hba2 = hab2;
    ha = hab1 .* hba1 .* hab2 .* hba2;

    var_h_true(var_h_idx) = var(real(ha)) + var(imag(ha));
    var_h_est(var_h_idx) = 4 * ((2 * var_h) ^ 4);
end

figure
plot(2*var_h_arr, var_h_true,'-',LineWidth=2)
hold on
plot(2*var_h_arr, var_h_est,':',LineWidth=3)
xlim([-0.01,1.01])
ylim([-0.05,4.2])
xlabel({'\itr'; '\rm(a)'},'FontSize',gca_size) 
ylabel('Variance of \itH_{FA}','FontSize',gca_size) 
set(gca,'FontSize',gca_size);
legend('Numerical Simulation', 'Gaussian Approximation', 'Fontsize', 20)

% save('var_h_true.mat','var_h_true')
% save('var_h_est.mat','var_h_est')

%% H_FA的分布估计准确性
var_h = 0.05;
hab1 = gen_complex_h(var_h, data_size);
hba1 = hab1;
hab2 = gen_complex_h(var_h, data_size);
hba2 = hab2;
ha = hab1 .* hba1 .* hab2 .* hba2;
var_h_est = 4 * ((2 * var_h) ^ 4);
var_real_part_h = var_h_est/2;

x_limits = [-3 3];

figure
sample_hist = histogram(real(ha), 'Normalization', 'probability', 'BinLimits', x_limits);
bin_egd = sample_hist.BinEdges;
bin_val = sample_hist.Values;
mul_par = 1 / sum(bin_val*(bin_egd(2) - bin_egd(1)));
x_val = (bin_egd(2:end)+bin_egd(1:end-1))/2;

y = (1/sqrt(2*pi*var_real_part_h))*exp(- x_val.^2 / (2*(var_real_part_h)));
pdf_h_true = mul_par*bin_val;
pdf_h_est = y;

figure
plot(x_val, pdf_h_true,'-', LineWidth=2)
hold on
plot(x_val, pdf_h_est, ':',LineWidth=3)
xlim([-0.05, 0.05])
ylim([0, 310])
xlabel({'Real Part of \itH_{FA}';'\rm(a)'})
ylabel('PDF')
set(gca,'FontSize',gca_size);
legend('Numerical Simulation','Gaussian Approximation', 'Fontsize', 20)

% save('pdf_h_true.mat','pdf_h_true')
% save('pdf_h_est.mat','pdf_h_est')

%% n的方差估计准确性
snr_arr = 0:25;
var_h = 0.05;
parfor snr_idx = 1:length(snr_arr)

    snr_val = snr_arr(snr_idx);
    var_n = (2 * var_h / (10 ^ (snr_val / 10)))/2;
    
    hab1 = gen_complex_h(var_h, data_size);
    hba1 = hab1;
    hab2 = gen_complex_h(var_h, data_size);
    hba2 = hab2;
    nb11 = gen_complex_n(var_n, data_size);
    na12 = gen_complex_n(var_n, data_size);
    nb23 = gen_complex_n(var_n, data_size);
    na24 = gen_complex_n(var_n, data_size);

    na = hba1 .* hab2 .* hba2 .* nb11 + hab2 .* hba2 .* na12 + hba2 .* nb23 + na24;

    var_n_true(snr_idx) = var(real(na)) + var(imag(na));
    var_n_est(snr_idx) = 2 * ((2 * var_h) ^ 3) * 2 * var_n + 2 * ((2 * var_h) ^ 2) * 2 * var_n + (2 * var_h) * 2 * var_n + 2 * var_n;
end

figure 
plot(snr_arr, var_n_true, '-', LineWidth=2)
hold on
plot(snr_arr, var_n_est, ':', LineWidth=3)
xlim([0, 25])
ylim([-0.002,0.12])
xlabel({'SNR'; '(b)'},'FontSize',gca_size) 
ylabel('Variance of \itn_{FA}','FontSize',gca_size) 
set(gca,'FontSize',gca_size);
legend('Numerical Simulation', 'Gaussian Approximation', 'Fontsize', 20)

% save('var_n_true.mat','var_n_true')
% save('var_n_est.mat','var_n_est')

%% n的分布估计准确性
% Low SNR
var_h = 0.05;
snr_val = 0;
var_n = (2 * var_h / (10 ^ (snr_val / 10)))/2;

hab1 = gen_complex_h(var_h, data_size);
hba1 = hab1;
hab2 = gen_complex_h(var_h, data_size);
hba2 = hab2;
nb11 = gen_complex_n(var_n, data_size);
na12 = gen_complex_n(var_n, data_size);
nb23 = gen_complex_n(var_n, data_size);
na24 = gen_complex_n(var_n, data_size);
na21 = gen_complex_n(var_n, data_size);
nb22 = gen_complex_n(var_n, data_size);
na13 = gen_complex_n(var_n, data_size);
nb14 = gen_complex_n(var_n, data_size);

na = hba1 .* hab2 .* hba2 .* nb11 + hab2 .* hba2 .* na12 + hba2 .* nb23 + na24;
var_n_est = 2 * ((2 * var_h) ^ 3) * 2 * var_n + 2 * ((2 * var_h) ^ 2) * 2 * var_n + (2 * var_h) * 2 * var_n + 2 * var_n;
var_real_part_n = var_n_est / 2;

x_limits = [-3 3];

figure
sample_hist = histogram(real(na), 'Normalization', 'probability', 'BinLimits', x_limits);
bin_egd = sample_hist.BinEdges;
bin_val = sample_hist.Values;
mul_par = 1 / sum(bin_val*(bin_egd(2) - bin_egd(1)));
x_val = (bin_egd(2:end)+bin_egd(1:end-1))/2;

y = (1/sqrt(2*pi*var_real_part_n))*exp(- x_val.^2 / (2*(var_real_part_n)));

pdf_n_true_low_snr = mul_par*bin_val;
pdf_n_est_low_snr = y;

figure
plot(x_val, pdf_n_true_low_snr,'-', LineWidth=2)
hold on 
plot(x_val, pdf_n_est_low_snr,':',LineWidth=3)
xlim([-1, 1])
ylim([0,2.1])
xlabel({'Real Part of \itn_{FA}'; '\rm(b)'})
ylabel('PDF')
set(gca,'FontSize',gca_size);
legend('Numerical Simulation', 'Gaussian Approximation', 'Fontsize', 20)

% save('pdf_n_true_low_snr.mat','pdf_n_true_low_snr')
% save('pdf_n_est_low_snr.mat','pdf_n_est_low_snr')

% High SNR
var_h = 0.05;
snr_val = 25;
var_n = (2 * var_h / (10 ^ (snr_val / 10)))/2;

hab1 = gen_complex_h(var_h, data_size);
hba1 = hab1;
hab2 = gen_complex_h(var_h, data_size);
hba2 = hab2;
nb11 = gen_complex_n(var_n, data_size);
na12 = gen_complex_n(var_n, data_size);
nb23 = gen_complex_n(var_n, data_size);
na24 = gen_complex_n(var_n, data_size);
na21 = gen_complex_n(var_n, data_size);
nb22 = gen_complex_n(var_n, data_size);
na13 = gen_complex_n(var_n, data_size);
nb14 = gen_complex_n(var_n, data_size);

na = hba1 .* hab2 .* hba2 .* nb11 + hab2 .* hba2 .* na12 + hba2 .* nb23 + na24;
var_n_est = 2 * ((2 * var_h) ^ 3) * 2 * var_n + 2 * ((2 * var_h) ^ 2) * 2 * var_n + (2 * var_h) * 2 * var_n + 2 * var_n;
var_real_part_n = var_n_est / 2;

x_limits = [-3 3];

figure
sample_hist = histogram(real(na), 'Normalization', 'probability', 'BinLimits', x_limits);
bin_egd = sample_hist.BinEdges;
bin_val = sample_hist.Values;
mul_par = 1 / sum(bin_val*(bin_egd(2) - bin_egd(1)));
x_val = (bin_egd(2:end)+bin_egd(1:end-1))/2;

y = (1/sqrt(2*pi*var_real_part_n))*exp(- x_val.^2 / (2*(var_real_part_n)));

pdf_n_true_high_snr = mul_par*bin_val;
pdf_n_est_high_snr = y;

figure
plot(x_val, pdf_n_true_high_snr, '-', LineWidth=2)
hold on 
plot(x_val, pdf_n_est_high_snr, ':', LineWidth=3)
xlim([-0.07, 0.07])
ylim([0,37])
xlabel({'Real Part of \itn_{FA}'; '\rm(c)'})
ylabel('PDF')
set(gca,'FontSize',gca_size);
legend('Numerical Simulation', 'Gaussian Approximation', 'Fontsize', 20)

% save('pdf_n_true_high_snr.mat','pdf_n_true_high_snr')
% save('pdf_n_est_high_snr.mat','pdf_n_est_high_snr')

toc

function h = gen_complex_h(var_h, data_size)
    real_part = sqrt(var_h)*randn(1, data_size);
    imag_part = sqrt(var_h)*randn(1, data_size);
    h = real_part + 1j*imag_part;
end
% 
function n = gen_complex_n(var_n, data_size)
    real_part = sqrt(var_n)*randn(1, data_size);
    imag_part = sqrt(var_n)*randn(1, data_size);
    n = real_part + 1j*imag_part;
end