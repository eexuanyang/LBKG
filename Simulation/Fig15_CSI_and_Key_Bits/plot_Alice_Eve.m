clear all; close all; clc

load('CSI_FDD_Alice.mat')
CSI_A = DMRS_arr_proc;
load('CSI_FDD_Bob.mat')
CSI_B = DMRS_arr_proc;
load('CSI_FDD_Eve.mat')
CSI_E = DMRS_arr_proc;

subcarr_idx = 426:1625;

figure
plot(subcarr_idx, CSI_A,'r',LineWidth=1.5)
hold on
plot(subcarr_idx, CSI_B,':k',LineWidth=2.5)
hold on
plot(subcarr_idx, CSI_E,'--b',LineWidth=3)
legend('Alice', 'Bob', 'Eve', 'Fontsize', 22)
xlabel({'Subcarrier Index','(a) \itZ\rm-multiplied CSI'}, 'Fontsize', 22)
ylabel('CSI Amplitude', 'Fontsize', 22)
xlim([426,1625])
set(gca,'FontSize',22)

load('CSI_Seg_FDD_Alice.mat')
CSI_A_seg = temp_arr_norm;
load('CSI_Seg_FDD_Bob.mat')
CSI_B_seg = temp_arr_norm;
load('CSI_Seg_FDD_Eve.mat')
CSI_E_seg = temp_arr_norm;

x_label_low_bound = round((426+1625)/2) - 30;
x_label_high_bound = round((426+1625)/2) + 29;

CSI_idx_low_bound = round((1+1200)/2) - 30;
CSI_idx_high_bound = round((1+1200)/2) + 29;

x_label_index = x_label_low_bound : x_label_high_bound;
CSI_idx = CSI_idx_low_bound : CSI_idx_high_bound;

figure
plot(x_label_index, CSI_A_seg(CSI_idx),'r',LineWidth=2)
hold on
plot(x_label_index, CSI_B_seg(CSI_idx),':k',LineWidth=3)
legend('Alice', 'Bob', 'Fontsize', 18)
xlabel({'Subcarrier Index','(b) Amplitude-normalized CSI w/ \itZ'}, 'Fontsize', 22)
ylabel('CSI Amplitude', 'Fontsize', 22)
set(gca,'FontSize',22)
xlim([x_label_low_bound, x_label_high_bound])
ylim([-0.05,1.05])

load('CSI_Seg_FDD_Alice_wo_Z.mat')
CSI_A_seg_wo_Z = temp_arr_norm;

figure
plot(x_label_index, CSI_A_seg(CSI_idx),'r',LineWidth=2)
hold on
plot(x_label_index, CSI_E_seg(CSI_idx),'--b',LineWidth=3)
legend('Alice', 'Eve', 'Fontsize', 18)
xlabel({'Subcarrier Index','(c) Amplitude-normalized CSI w/ \itZ'}, 'Fontsize', 22)
ylabel('CSI Amplitude', 'Fontsize', 22)
set(gca,'FontSize',22)
xlim([x_label_low_bound, x_label_high_bound])
ylim([-0.05,1.05])

figure
plot(x_label_index, CSI_A_seg_wo_Z(CSI_idx),'r',LineWidth=2)
hold on
plot(x_label_index, CSI_E_seg(CSI_idx),'--b',LineWidth=3)
legend('Alice', 'Eve', 'Fontsize', 22)
xlabel({'Subcarrier Index','(d) Amplitude-normalized CSI w/o \itZ'}, 'Fontsize', 22)
ylabel('CSI Amplitude', 'Fontsize', 22)
set(gca,'FontSize',22)
xlim([x_label_low_bound, x_label_high_bound])
ylim([-0.05,1.05])

load('key_Alice.mat')
key_A = key_arr_fdd_a_norm_sort_ds_bit_disc_final;
load('key_Bob.mat')
key_B = key_arr_fdd_a_norm_sort_ds_bit_disc_final;
load('key_Eve.mat')
key_E = key_arr_fdd_e_norm_sort_ds_bit_disc_final;
load('key_Alice_wo_Z.mat')
key_A_wo_Z = key_arr_fdd_a_norm_sort_ds_bit_disc_final;

key_A_temp(key_A == 1) = 0.7;
key_A_temp(key_A == 0) = 0.3;

key_B_temp(key_B == 1) = 0.7;
key_B_temp(key_B == 0) = 0.3;

key_E_temp(key_E == 1) = 0.7;
key_E_temp(key_E == 0) = 0.3;

key_A_wo_Z_temp(key_A_wo_Z == 1) = 0.7;
key_A_wo_Z_temp(key_A_wo_Z == 0) = 0.3;

figure
plot(1:128, key_A,'or',LineWidth=1.5);
hold on
plot(1:128, key_B,'*K',LineWidth=1);
hold on
h1 = plot(18:3:106, key_A_temp(50:79),'or','MarkerSize',16,'MarkerFaceColor','none','MarkerEdgeColor','r','LineWidth',1.5);
hold on 
h2 = plot(18:3:106, key_B_temp(50:79),'*K','MarkerSize',16,'MarkerFaceColor','none','MarkerEdgeColor','k','LineWidth',1.5);
hold on
legend([h1, h2], {'Alice', 'Bob'}, 'FontSize', 35)
xlabel({'Key Bits Index','(e) Key bits w/ segmentation'}, 'Fontsize', 35)
ylabel('Key Value', 'Fontsize', 35)
xlim([1,128])
ylim([-0.05, 1.05])
set(gca,'FontSize',35)

figure
plot(key_A,'Or',LineWidth=1.5);
hold on
plot(key_E,'*b',LineWidth=1);
h1 = plot(18:3:106, key_A_temp(50:79),'or','MarkerSize',16,'MarkerFaceColor','none','MarkerEdgeColor','r','LineWidth',1.5);
hold on 
h2 = plot(18:3:106, key_E_temp(50:79),'*b','MarkerSize',16,'MarkerFaceColor','none','MarkerEdgeColor','b','LineWidth',1.5);
hold on
legend([h1, h2], {'Alice', 'Eve'}, 'FontSize', 35)
xlabel({'Key Bits Index','(f)  Key bits w/ segmentation'}, 'Fontsize', 35)
ylabel('Key Value', 'Fontsize', 35)
xlim([1,128])
ylim([-0.05, 1.05])
set(gca,'FontSize', 35)

figure
plot(key_A_wo_Z,'Or',LineWidth=1.5);
hold on
plot(key_E,'*b',LineWidth=1);
h1 = plot(18:3:106, key_A_wo_Z_temp(50:79),'or','MarkerSize',16,'MarkerFaceColor','none','MarkerEdgeColor','r','LineWidth',1.5);
hold on 
h2 = plot(18:3:106, key_E_temp(50:79),'*b','MarkerSize',16,'MarkerFaceColor','none','MarkerEdgeColor','b','LineWidth',1.5);
hold on
lgd = legend([h1, h2], {'Alice', 'Eve'}, 'FontSize', 35);
xlabel({'Key Bits Index','(g)  Key bits w/o segmentation'}, 'Fontsize', 35)
ylabel('Key Value', 'Fontsize', 35)
xlim([1,128])
ylim([-0.05, 1.05])
set(gca,'FontSize', 35)

