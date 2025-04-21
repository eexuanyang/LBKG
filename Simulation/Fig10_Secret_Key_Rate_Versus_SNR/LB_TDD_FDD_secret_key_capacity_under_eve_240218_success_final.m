clear all;
close all;
clc
r = 1;
snr_arr = 0:2:30;
d_arr = 0:0.1:0.2;
I_f = zeros(length(snr_arr), length(d_arr));
I_t = zeros(length(snr_arr), length(d_arr));
I_diff = zeros(length(snr_arr), length(d_arr));
a1 = 1; b2 = 1; 
c = 1;
lambda_1 = 0.2;
lambda_2 = 0.2;
for snr_idx = 1:length(snr_arr)
    for d_idx = 1:length(d_arr)
        n = r/(10^(snr_arr(snr_idx)/10));
        d = d_arr(d_idx);
        a2 = a1 - d; b1 = b2 - d;

        r_taa = a1^2*b2^2*r^2 + b2^2*r*n + n;
        r_tbb = a2^2*b1^2*r^2 + a2^2*r*n + n;
        r_tee = a1^2*c^2*r^2 + c^2*r*n + n;
        r_tab = a1*a2*b1*b2*r^2;
        r_tae = lambda_1*a1^2*b2*c*r^2 + lambda_1*b2*c*r*n;
        r_tbe = lambda_1*a1*a2*b1*c*r^2;

        r_hat_tae = r_taa*r_tee - r_tae^2;
        r_hat_tbe = r_tbb*r_tee - r_tbe^2;
        r_hat_tee = r_tee;
        r_hat_tabe= det([r_taa r_tab r_tae; r_tab r_tbb r_tbe; r_tae r_tbe r_tee]);

        I_t(snr_idx, d_idx) = log2((r_hat_tae*r_hat_tbe)/(r_hat_tee*r_hat_tabe));
        
        r_faa = 4*a1^2*a2^2*b1^2*b2^2*r^4 + 2*b1^2*a2^2*b2^2*r^3*n + 2*a2^2*b2^2*r^2*n + b2^2*r*n + n;
        r_fbb = 4*a1^2*a2^2*b1^2*b2^2*r^4 + 2*a1^2*b1^2*a2^2*r^3*n + 2*a1^2*b1^2*r^2*n + a1^2*r*n + n;
        r_fee = 2*(1+lambda_2^2)*a1^2*a2^2*b1^2*c^2*r^4 + (1+lambda_2^2)*a2^2*b1^2*c^2*r^3*n + (1+lambda_2^2)*a2^2*c^2*r^2*n + c^2*r*n + n;
        r_fab = 4*a1^2*a2^2*b1^2*b2^2*r^4;
        r_fae = 4*lambda_2*a1^2*a2^2*b1^2*b2*c*r^4 + 2*lambda_2*a2^2*b1^2*b2*c*r^3*n + 2*lambda_2*a2^2*b2*c*r^2*n + lambda_2*b2*c*r*n;
        r_fbe = 4*lambda_2*a1^2*a2^2*b1^2*b2*c*r^4;

        r_hat_fae = r_faa*r_fee - r_fae^2;
        r_hat_fbe = r_fbb*r_fee - r_fbe^2;
        r_hat_fee = r_fee;
        fabe = [r_faa r_fab r_fae; r_fab r_fbb r_fbe; r_fae r_fbe r_fee];
        r_hat_fabe= det(fabe);

        I_f(snr_idx, d_idx) = log2((r_hat_fae*r_hat_fbe)/(r_hat_fee*r_hat_fabe));

        I_diff(snr_idx, d_idx) = I_f(snr_idx, d_idx) - I_t(snr_idx, d_idx);
    end
end

LineWidth_val = 4; MarkerSize_val = 10;
figure
plot(snr_arr,I_f(:,1),'r-^',LineWidth=LineWidth_val, MarkerSize=MarkerSize_val);
hold on
plot(snr_arr,I_f(:,2),'b-s',LineWidth=LineWidth_val, MarkerSize=MarkerSize_val);
hold on
plot(snr_arr,I_f(:,3),'c-o',LineWidth=LineWidth_val, MarkerSize=MarkerSize_val);
hold on
plot(snr_arr,I_t(:,1),'g-v',LineWidth=LineWidth_val, MarkerSize=MarkerSize_val);
hold on
plot(snr_arr,I_t(:,2),'m-*',LineWidth=LineWidth_val, MarkerSize=MarkerSize_val);
hold on
plot(snr_arr,I_t(:,3),'k-p',LineWidth=LineWidth_val, MarkerSize=MarkerSize_val);
hold on
ylim([0,max(I_f(:,1))])
xlabel({'SNR (dB)';'(a)'}, 'Fontsize', 30);
ylabel('{\it{I}}_{\it{F,sk}} and {\it{I}}_{\it{T,sk}} (bit/symbol)', 'Fontsize', 30);
set(gca,'FontSize',30,'FontName', 'Times New Roman');
legend('LBKG {\it{d}} = 0','LBKG, {\it{d}} = 0.1','LBKG, {\it{d}} = 0.2','LB-TDD, {\it{d}} = 0','LB-TDD, {\it{d}} = 0.1','LB-TDD, {\it{d}} = 0.2', 'Fontsize', 26);

LineWidth_val = 3;
figure
plot(snr_arr,I_diff(:,1),'k-^',LineWidth=LineWidth_val, MarkerSize=MarkerSize_val);
hold on
plot(snr_arr,I_diff(:,2),'b-s',LineWidth=LineWidth_val, MarkerSize=MarkerSize_val);
hold on
plot(snr_arr,I_diff(:,3),'r-v',LineWidth=LineWidth_val, MarkerSize=MarkerSize_val);
hold on
yticks([0,0.1,0.2,0.3,0.4,0.5])
xlabel({'SNR (dB)';'(b)'}, 'Fontsize', 30);
ylabel('\Delta {\it{I}}_{\it{sk}} (bit/symbol)', 'Fontsize', 30);
set(gca,'FontSize',30,'FontName', 'Times New Roman');
legend('{\it{d}} = 0','{\it{d}} = 0.1','{\it{d}} = 0.2', 'Fontsize', 26);
