clear all; close all; clc

load("MIs_arr_simple.mat")
load("Gaussian_estimate_simple.mat")

gca_size = 22;

MIs_arr_simple_combine_mean = mean(MIs_arr_simple,2); 

MIs_arr_tifs_bit = log2(exp(1)) * MIs_arr_simple_combine_mean; % ln -> log2
Gaussian_estimate_tifs_bit = log2(exp(1)) * Gaussian_estimate_simple; % ln -> log2

LineWidth_val = 2; MarkerSize_val = 8;
figure
plot((0:30),MIs_arr_tifs_bit,'-^',LineWidth=LineWidth_val, MarkerSize=MarkerSize_val);
hold on 
plot((0:30),Gaussian_estimate_tifs_bit,'-s',LineWidth=LineWidth_val, MarkerSize=MarkerSize_val);
yticks([0:2:10])
xlabel({"SNR"; '(a)'})
ylabel({'{\it{I}} (bits/symbol)'})
legend("MINE", "Theory")
% xlim([0 26]);
% ylim([-0.01 0.5])
set(gca,'FontSize',gca_size);