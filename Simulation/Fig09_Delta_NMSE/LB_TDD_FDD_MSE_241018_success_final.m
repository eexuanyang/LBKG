clear all;
% close all;
clc
r = 1;
snr_arr = -30:5:30;
d_arr = 0:0.05:0.2;

mse_fdd = zeros(length(snr_arr), length(d_arr));
mse_tdd = zeros(length(snr_arr), length(d_arr));
mse_jrnso = zeros(length(snr_arr), length(d_arr));
mse_arr2 = zeros(length(snr_arr), length(d_arr));
mse_arr = zeros(length(snr_arr), length(d_arr));
a1 = 1; b2 = 1;
for snr_idx = 1:length(snr_arr)
    for d_idx = 1:length(d_arr)
        n = r/(10^(snr_arr(snr_idx)/10));
        d = d_arr(d_idx);
        a2 = a1 - d; b1 = b2 - d;

        mse_fdd(snr_idx, d_idx) = (2*a2^2*b1^2*(a1^2+b2^2)*r^3*n + 2*(a1^2*b1^2 + a2^2*b2^2)*r^2*n + (a1^2+b2^2)*r*n + 2*n)...
                                 /(4*a1^2*a2^2*b1^2*b2^2*r^4 + 2*a2^2*b1^2*b2^2*r^3*n + 2*a2^2*b2^2*r^2*n + b2^2*r*n + n);
        mse_tdd(snr_idx, d_idx) = ((a1*b2 - a2*b1)^2*r^2 + (a2^2 + b2^2)*r*n + 2*n)...
                                 /(a1^2*b2^2*r^2 + b2^2*r*n + n);
        mse_jrnso(snr_idx, d_idx) = (a1^2*r*n + b2^2*r*n + 2*n) / (a1^2*b2^2*r^2 + b2^2*r*n + n);
        mse_arr(snr_idx, d_idx) = mse_fdd(snr_idx, d_idx) - mse_tdd(snr_idx, d_idx);
        mse_arr2(snr_idx, d_idx) = mse_fdd(snr_idx, d_idx) - mse_jrnso(snr_idx, d_idx);

    end
end

LineWidth_val = 2; MarkerSize_val = 8;

figure
plot(snr_arr,mse_arr(:,1),'k-^',LineWidth=LineWidth_val, MarkerSize=MarkerSize_val);
hold on
plot(snr_arr,mse_arr(:,2),'b-s',LineWidth=LineWidth_val, MarkerSize=MarkerSize_val);
hold on
plot(snr_arr,mse_arr(:,3),'c-o',LineWidth=LineWidth_val, MarkerSize=MarkerSize_val);
hold on
plot(snr_arr,mse_arr(:,4),'g-*',LineWidth=LineWidth_val, MarkerSize=MarkerSize_val);
hold on
plot(snr_arr,mse_arr(:,5),'r-v',LineWidth=LineWidth_val, MarkerSize=MarkerSize_val);
xlim([snr_arr(1), snr_arr(end)])
ylim([-0.14, 0.19])
xticks([-30, -15, 0, 15, 30])
xlabel('SNR (dB)', 'Fontsize', 18);
ylabel('\DeltaNMSE', 'Fontsize', 18);
set(gca,'FontSize',18,'FontName', 'Times New Roman');
legend('{\it{d}} = 0','{\it{d}} = 0.05','{\it{d}} = 0.10','{\it{d}} = 0.15','{\it{d}} = 0.2', 'Fontsize', 18);
