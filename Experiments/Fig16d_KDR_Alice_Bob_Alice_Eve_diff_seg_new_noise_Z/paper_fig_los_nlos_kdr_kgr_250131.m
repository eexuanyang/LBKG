clear all
close all
clc

load("key_fdd_norm_sort_kdr_zz_seg10.mat")
A_B_seg_10 = key_fdd_norm_sort_kdr_zz;
load("key_fdd_norm_sort_kdr_zz_seg20.mat")
A_B_seg_20 = key_fdd_norm_sort_kdr_zz;
load("key_fdd_norm_sort_kdr_zz_seg50.mat")
A_B_seg_50 = key_fdd_norm_sort_kdr_zz;
load("key_fdd_norm_sort_kdr_zz_seg100.mat")
A_B_seg_100 = key_fdd_norm_sort_kdr_zz;
load("key_fdd_norm_sort_kdr_zz_seg150.mat")
A_B_seg_150 = key_fdd_norm_sort_kdr_zz;

A_B_seg_10_mean   = mean(A_B_seg_10);
A_B_seg_20_mean  = mean(A_B_seg_20);
A_B_seg_50_mean  = mean(A_B_seg_50);
A_B_seg_100_mean  = mean(A_B_seg_100);
A_B_seg_150_mean = mean(A_B_seg_150);

A_B_seg_10_mean_arr   = [A_B_seg_10_mean 0];
A_B_seg_20_mean_arr  = [A_B_seg_20_mean 0];
A_B_seg_50_mean_arr  = [A_B_seg_50_mean 0];
A_B_seg_100_mean_arr  = [A_B_seg_100_mean 0];
A_B_seg_150_mean_arr = [A_B_seg_150_mean 0];

kdr = [A_B_seg_10_mean_arr; A_B_seg_20_mean_arr; A_B_seg_50_mean_arr; A_B_seg_100_mean_arr; A_B_seg_150_mean_arr];

figure
yyaxis left
b=bar(kdr);
ch = get(b,'children');
set(gca,'XTickLabel',{'10', '20', '50', '100', '150'},'FontSize', 50)

% x = 1:size(kdr,1);
% for i = 1:size(kdr,1)
%     text(x(i)-0.22, kdr(i,1), sprintf('%.4f', kdr(i,1)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 24);
% end
yticks([0, 0.05, 0.1, 0.15, 0.2, 0.25])
ylim([0, 0.25]);
ylabel('Alice-Bob KDR','FontSize', 50);

load("key_fdd_norm_sort_A_E_kdr_zz_seg10.mat")
A_E_seg_10 = key_fdd_norm_sort_A_E_kdr_zz;
load("key_fdd_norm_sort_A_E_kdr_zz_seg20.mat")
A_E_seg_20 = key_fdd_norm_sort_A_E_kdr_zz;
load("key_fdd_norm_sort_A_E_kdr_zz_seg50.mat")
A_E_seg_50 = key_fdd_norm_sort_A_E_kdr_zz;
load("key_fdd_norm_sort_A_E_kdr_zz_seg100.mat")
A_E_seg_100 = key_fdd_norm_sort_A_E_kdr_zz;
load("key_fdd_norm_sort_A_E_kdr_zz_seg150.mat")
A_E_seg_150 = key_fdd_norm_sort_A_E_kdr_zz;

A_E_seg_10_mean   = mean(A_E_seg_10);
A_E_seg_20_mean  = mean(A_E_seg_20);
A_E_seg_50_mean  = mean(A_E_seg_50);
A_E_seg_100_mean  = mean(A_E_seg_100);
A_E_seg_150_mean = mean(A_E_seg_150);

A_E_seg_10_mean_arr   = [0 A_E_seg_10_mean];
A_E_seg_20_mean_arr  = [0 A_E_seg_20_mean];
A_E_seg_50_mean_arr  = [0 A_E_seg_50_mean];
A_E_seg_100_mean_arr  = [0 A_E_seg_100_mean];
A_E_seg_150_mean_arr = [0 A_E_seg_150_mean];

kdr_A_E = [A_E_seg_10_mean_arr; A_E_seg_20_mean_arr; A_E_seg_50_mean_arr; A_E_seg_100_mean_arr; A_E_seg_150_mean_arr];

yyaxis right
b=bar(kdr_A_E);
ch = get(b,'children');
% x = 1:size(kdr_A_E,1);
% for i = 1:size(kdr_A_E,1)
%     text(x(i)+0.22, kdr_A_E(i,2), sprintf('%.4f', kdr_A_E(i,2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 24);
% end
% ylim([0, 2.2]);
yticks([0, 0.1, 0.2, 0.3, 0.4])
xlabel('(d) Segmentation number','FontSize', 52);
ylabel('Alice-Eve KDR','FontSize', 50);
% set(gca,'FontSize',30)
