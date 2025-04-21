clear all
close all
clc

load("key_fdd_norm_sort_kdr_zz_range1.mat")
A_B_rg_1 = key_fdd_norm_sort_kdr_zz;
load("key_fdd_norm_sort_kdr_zz_range2.mat")
A_B_rg_2 = key_fdd_norm_sort_kdr_zz;
load("key_fdd_norm_sort_kdr_zz_range3.mat")
A_B_rg_3 = key_fdd_norm_sort_kdr_zz;
load("key_fdd_norm_sort_kdr_zz_range4.mat")
A_B_rg_4 = key_fdd_norm_sort_kdr_zz;
load("key_fdd_norm_sort_kdr_zz_range5.mat")
A_B_rg_5 = key_fdd_norm_sort_kdr_zz;

A_B_rg_1_mean = mean(A_B_rg_1);
A_B_rg_2_mean = mean(A_B_rg_2);
A_B_rg_3_mean = mean(A_B_rg_3);
A_B_rg_4_mean = mean(A_B_rg_4);
A_B_rg_5_mean = mean(A_B_rg_5);

% A_B_rg_25_mean_arr   = [A_B_rg_25_mean 0];
A_B_rg_1_mean_arr  = [A_B_rg_1_mean 0];
A_B_rg_2_mean_arr  = [A_B_rg_2_mean 0];
A_B_rg_3_mean_arr  = [A_B_rg_3_mean 0];
A_B_rg_4_mean_arr = [A_B_rg_4_mean 0];
A_B_rg_5_mean_arr = [A_B_rg_5_mean 0];

kdr = [A_B_rg_1_mean_arr; A_B_rg_2_mean_arr; A_B_rg_3_mean_arr; A_B_rg_4_mean_arr; A_B_rg_5_mean_arr];

figure
yyaxis left
b=bar(kdr);
ch = get(b,'children');
% x = 1:size(kdr,1);
% for i = 1:size(kdr,1)
%     text(x(i)-0.22, kdr(i,1), sprintf('%.4f', kdr(i,1)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 24);
% end
ylim([0, 0.1]);
% set(gca,'XTickLabel',{'{\itRg} = 0.04', '{\itRg} = 0.06', '{\itRg} = 0.08', '{\itRg} = 0.1', '{\itRg} = 0.12'},'FontSize',30)
set(gca,'XTickLabel',{'0.04', '0.06', '0.08', '0.1', '0.12'},'FontSize',50)
ylabel('Alice-Bob KDR','FontSize',50);
yticks([0, 0.025, 0.05, 0.075, 0.1])
% set(gca,'FontSize',38)

load("key_fdd_norm_sort_A_E_kdr_zz_range1.mat")
A_E_rg_1 = key_fdd_norm_sort_A_E_kdr_zz;
load("key_fdd_norm_sort_A_E_kdr_zz_range2.mat")
A_E_rg_2 = key_fdd_norm_sort_A_E_kdr_zz;
load("key_fdd_norm_sort_A_E_kdr_zz_range3.mat")
A_E_rg_3 = key_fdd_norm_sort_A_E_kdr_zz;
load("key_fdd_norm_sort_A_E_kdr_zz_range4.mat")
A_E_rg_4 = key_fdd_norm_sort_A_E_kdr_zz;
load("key_fdd_norm_sort_A_E_kdr_zz_range5.mat")
A_E_rg_5 = key_fdd_norm_sort_A_E_kdr_zz;

A_E_rg_1_mean = mean(A_E_rg_1);
A_E_rg_2_mean = mean(A_E_rg_2);
A_E_rg_3_mean = mean(A_E_rg_3);
A_E_rg_4_mean = mean(A_E_rg_4);
A_E_rg_5_mean = mean(A_E_rg_5);

% A_E_rg_25_mean_arr   = [0 A_E_rg_25_mean];
A_E_rg_1_mean_arr = [0 A_E_rg_1_mean];
A_E_rg_2_mean_arr = [0 A_E_rg_2_mean];
A_E_rg_3_mean_arr = [0 A_E_rg_3_mean];
A_E_rg_4_mean_arr = [0 A_E_rg_4_mean];
A_E_rg_5_mean_arr = [0 A_E_rg_5_mean];

kdr_A_E = [A_E_rg_1_mean_arr; A_E_rg_2_mean_arr; A_E_rg_3_mean_arr; A_E_rg_4_mean_arr; A_E_rg_5_mean_arr];

yyaxis right
b=bar(kdr_A_E);
ch = get(b,'children');

% x = 1:size(kdr_A_E,1);
% for i = 1:size(kdr_A_E,1)
%     text(x(i)+0.22, kdr_A_E(i,2), sprintf('%.4f', kdr_A_E(i,2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 24);
% end
ylim([0, 0.5]);
yticks([0, 0.1, 0.2, 0.3, 0.4, 0.5])
xlabel('(c) Amplitude range of \itZ','FontSize',52);
ylabel('Alice-Eve KDR','FontSize',50);
% set(gca,'FontSize',30)
