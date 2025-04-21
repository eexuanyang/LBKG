clear all
close all
clc

load("key_fdd_norm_sort_kdr_zz_seg100.mat")
load("key_fdd_norm_sort_kdr_wo_Z.mat")
load("key_tdd_norm_sort_kdr.mat")
load("key_sar_norm_sort_kdr.mat")
load("kdr_MV_300_DL_1.mat")
load("kdr_MV_300_DL_2.mat")

fdd_kdr_term_z = mean(key_fdd_norm_sort_kdr_zz);
fdd_kdr = mean(key_fdd_norm_sort_kdr_wo_Z);
tdd_kdr = key_tdd_norm_sort_kdr;
sar_kdr = key_sar_norm_sort_kdr;
DL_kdr_1 = kdr_MV_300_DL_1(:);
DL_kdr_2 = kdr_MV_300_DL_2(:);
DL_kdr = mean([DL_kdr_1(1:5); DL_kdr_1(1:5)]);

a = [fdd_kdr_term_z 0];
b = [fdd_kdr 0];
c = [tdd_kdr 0];
d = [sar_kdr 0];
f = [DL_kdr 0];

kdr = [a;b;c;d;f];
% kgr = [fdd_kgr_LOS  tdd_kgr_LOS;
%        fdd_kgr_NLOS tdd_kgr_NLOS];
figure
yyaxis left
b=bar(kdr);
% ch = get(b,'children');
% set(gca,'XTickLabel',{'LBKG w/\newline  term \itZ', 'LBKG w/o\newline   term \itZ','LB-TDD','    SAR, \newline amplitude','     DL, \newline amplitude'},'FontSize',28)
set(gca,'XTickLabel',{' ', ' ', ' ', ' ', ' '},'FontSize',28,'FontName', 'Times New Roman');
title(' ', 'FontSize',48,'FontName', 'Times New Roman');
% % 调整刻度标签的位置
% ax = gca; % 获取当前坐标轴
% ax.XAxis.TickLabelRotation = 0; % 设置标签旋转角度（0 表示不旋转）
% 
% % 将第 2 和第 4 个标签显示在上方
% for i = 1:length(ax.XTickLabel)
%     if mod(i, 2) == 0 % 偶数位置（2, 4）
%         ax.XTickLabel{i} = ['\color{red}' ax.XTickLabel{i}]; % 设置为红色并显示在上方
%         ax.XRuler.TickLabelGapOffset = -10; % 调整标签位置（负值表示上方）
%     end
% end

% set(gca,'FontName','Times New Roman','FontSize',45);
% set(gca,'XTickLabel',{'LBKG w/ term z', 'LBKG w/o term z','LB-TDD','SAR','DL'},'FontSize',28)

% x = 1:size(kdr,1);
% for i = 1:size(kdr,1)
%     text(x(i)-0.22, kdr(i,1), sprintf('%.4f', kdr(i,1)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 28);
% end
ylim([0, 0.6]);
yticks([0.1:0.1:0.6]);
ax = gca;
set(ax.YAxis, 'FontSize', 50,'FontName', 'Times New Roman');

% legend('LB-FDD','LB-TDD','FontSize',14);
ylabel('KDR','FontSize', 50,'FontName', 'Times New Roman');
% set(gca,'FontSize',22)
 
load("key_fdd_norm_sort_kgr_zz_seg100.mat")
load("key_fdd_norm_sort_kgr_wo_Z.mat")
load("key_tdd_norm_sort_kgr.mat")
load("key_sar_norm_sort_kgr.mat")
load("kgr_MV_300_DL_1.mat")
load("kgr_MV_300_DL_2.mat")

fdd_kgr_term_z = mean(key_fdd_norm_sort_kgr_zz) * 128 / 150;
fdd_kgr = mean(key_fdd_norm_sort_kgr_wo_Z) * 128 / 150;
tdd_kgr = key_tdd_norm_sort_kgr * 128 / 150;
sar_kgr = key_sar_norm_sort_kgr * 128 / 150;
DL_kgr_1 = kgr_MV_300_DL_1(:);
DL_kgr_2 = kgr_MV_300_DL_2(:);
DL_kgr = mean([DL_kgr_1(1:5); DL_kgr_1(1:5)]) * 128 / 150;

a = [0 fdd_kgr_term_z];
b = [0 fdd_kgr];
c = [0 tdd_kgr];
d = [0 sar_kgr];
f = [0 DL_kgr];

kgr = [a;b;c;d;f];

yyaxis right
b=bar(kgr);
ch = get(b,'children');
% x = 1:size(kgr,1);
% for i = 1:size(kgr,1)
%     text(x(i)+0.22, kgr(i,2), sprintf('%.4f', kgr(i,2)), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 24);
% end
ylim([0, 2]);
yticks([0, 0.5, 1, 1.5, 2]);
% legend('LB-FDD','LB-TDD','FontSize',14);
xlabel('(a) Different schemes','FontSize', 50,'FontName', 'Times New Roman');
ylabel('KGR (bit/subcarrier)','FontSize', 50,'FontName', 'Times New Roman');
% set(gca,'FontSize',28)
ax = gca;
set(ax.YAxis, 'FontSize', 50,'FontName', 'Times New Roman');