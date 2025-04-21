clear all
close all
clc

load("key_fdd_norm_sort_kdr_zz_mv25.mat")
A_B_mv_25 = key_fdd_norm_sort_kdr_zz;
load("key_fdd_norm_sort_kdr_zz_mv75.mat")
A_B_mv_75 = key_fdd_norm_sort_kdr_zz;
load("key_fdd_norm_sort_kdr_zz_mv150.mat")
A_B_mv_150 = key_fdd_norm_sort_kdr_zz;
load("key_fdd_norm_sort_kdr_zz_mv200.mat")
A_B_mv_200 = key_fdd_norm_sort_kdr_zz;
load("key_fdd_norm_sort_kdr_zz_mv300.mat")
A_B_mv_300 = key_fdd_norm_sort_kdr_zz;

A_B_mv_25_mean   = mean(A_B_mv_25);
A_B_mv_75_mean  = mean(A_B_mv_75);
A_B_mv_150_mean  = mean(A_B_mv_150);
A_B_mv_200_mean  = mean(A_B_mv_200);
A_B_mv_300_mean = mean(A_B_mv_300);

A_B_mv_25_mean_arr   = [A_B_mv_25_mean 0];
A_B_mv_75_mean_arr  = [A_B_mv_75_mean 0];
A_B_mv_150_mean_arr  = [A_B_mv_150_mean 0];
A_B_mv_200_mean_arr  = [A_B_mv_200_mean 0];
A_B_mv_300_mean_arr = [A_B_mv_300_mean 0];

kdr = [A_B_mv_25_mean_arr; A_B_mv_75_mean_arr; A_B_mv_150_mean_arr; A_B_mv_200_mean_arr; A_B_mv_300_mean_arr];

figure
yyaxis left
b=bar(kdr);
ch = get(b,'children');
set(gca,'XTickLabel',{'25', '75', '150', '200', '300'},'FontSize',50)

% x = 1:size(kdr,1);
% for i = 1:size(kdr,1)
%     text(x(i)-0.22, kdr(i,1), sprintf('%.4f', kdr(i,1)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 24);
% end
yticks([0, 0.05, 0.1, 0.15, 0.2, 0.25])
ylim([0, 0.27]);
ylabel('Alice-Bob KDR','FontSize',50);

load("key_fdd_norm_sort_A_E_kdr_zz_mv25.mat")
A_E_mv_25 = key_fdd_norm_sort_A_E_kdr_zz;
load("key_fdd_norm_sort_A_E_kdr_zz_mv75.mat")
A_E_mv_75 = key_fdd_norm_sort_A_E_kdr_zz;
load("key_fdd_norm_sort_A_E_kdr_zz_mv150.mat")
A_E_mv_150 = key_fdd_norm_sort_A_E_kdr_zz;
load("key_fdd_norm_sort_A_E_kdr_zz_mv200.mat")
A_E_mv_200 = key_fdd_norm_sort_A_E_kdr_zz;
load("key_fdd_norm_sort_A_E_kdr_zz_mv300.mat")
A_E_mv_300 = key_fdd_norm_sort_A_E_kdr_zz;

A_E_mv_25_mean   = mean(A_E_mv_25);
A_E_mv_75_mean  = mean(A_E_mv_75);
A_E_mv_150_mean  = mean(A_E_mv_150);
A_E_mv_200_mean  = mean(A_E_mv_200);
A_E_mv_300_mean = mean(A_E_mv_300);

A_E_mv_25_mean_arr   = [0 A_E_mv_25_mean];
A_E_mv_75_mean_arr  = [0 A_E_mv_75_mean];
A_E_mv_150_mean_arr  = [0 A_E_mv_150_mean];
A_E_mv_200_mean_arr  = [0 A_E_mv_200_mean];
A_E_mv_300_mean_arr = [0 A_E_mv_300_mean];

kdr_A_E = [A_E_mv_25_mean_arr; A_E_mv_75_mean_arr; A_E_mv_150_mean_arr; A_E_mv_200_mean_arr; A_E_mv_300_mean_arr];
yyaxis right
b=bar(kdr_A_E);
ch = get(b,'children');
% x = 1:size(kdr_A_E,1);
% for i = 1:size(kdr_A_E,1)
%     text(x(i)+0.22, kdr_A_E(i,2), sprintf('%.4f', kdr_A_E(i,2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 24);
% end
ylim([0, 0.5]);
yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
xlabel('(b) Moving average filter length','FontSize',52);
ylabel('Alice-Eve KDR','FontSize',50);
% set(gca,'FontSize',30)
