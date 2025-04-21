% Plots from TIFS'10 paper.
close all
clear all
clc
tic
do_show_progress = 0;
do_flush = 0; % Set=0 for MATLAB, =1  for Octave
 
% Directory for caching results.
data_dir = 'data';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Figure 3: Single Antenna 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_equal = [];
data_adv = [];

% Equal SNR case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SNR = 15;
N = 1;
% Npath_range = [5];
Npath_range = [1 2 10];
d_max = 10;
d_min = 0.1;
d_N = 100;
d_range = 10.^((log10(d_max)-log10(d_min))/(d_N-1)*(0:d_N-1) + log10(d_min));

% Number of realizations (to get smooth curves) for different number of paths.
M_range = [10 5000 1000];

fname = strcat(data_dir, '/sing_ant_vs_d.mat');
if ~exist(fname)

mi_h_f = zeros(d_N, length(Npath_range));
mi_sa_f = mi_h_f;
mi_st_f = mi_h_f;
mi_h_t = zeros(d_N, length(Npath_range));
mi_sa_t = mi_h_t;
mi_st_t = mi_h_t;
 
for ii=1:length(Npath_range)
  for jj=1:length(d_range)
   
    if (do_show_progress)
      fprintf('[%d,%d] ', ii, jj);
      if (do_flush), fflush(1); end
    end
    
    clear ops
    ops.SNR1 = SNR;
    ops.SNR2 = SNR;
    ops.SNR3 = SNR;
    ops.N1 = 1;
    ops.N2 = 1;
    ops.N3 = 1;
    ops.d = d_range(jj);
    ops.Nray = Npath_range(ii);
    ops.M = M_range(ii);
    ops.diff = 0.2;
    [mi_h_f(jj,ii), mi_sa_f(jj,ii), mi_st_f(jj,ii), mi_h_t(jj,ii), mi_sa_t(jj,ii), mi_st_t(jj,ii)] = do_mi1_avg(ops);
  end
end

% save(fname, 'mi_h', 'mi_sa', 'mi_st');

else
  
% load(fname);

end

data_equal.mi_h_f = mi_h_f;
data_equal.mi_sa_f = mi_sa_f;
data_equal.mi_st_f = mi_st_f;

data_equal.mi_h_t = mi_h_t;
data_equal.mi_sa_t = mi_sa_t;
data_equal.mi_st_t = mi_st_t;

figure
for ii=1:length(Npath_range)
  semilogx(d_range, data_equal.mi_st_f(:, ii)./data_equal.mi_h_f(:, ii), 'r', LineWidth=1);
  hold on;
end
for ii=1:length(Npath_range)
  semilogx(d_range, data_equal.mi_st_t(:, ii)./data_equal.mi_h_t(:, ii), 'k', LineWidth=1);
  hold on;
end
hold off;
xlabel('Eavesdropper Separation w (wavelengths)');
ylabel('I_{VK}/I_{K}');
ylim([0,1]);

figure
h1 = semilogx(d_range, data_equal.mi_st_f(:, 1)./data_equal.mi_h_f(:, 1), '--r', LineWidth=1.5);
hold on;
h2 = semilogx(d_range, data_equal.mi_st_f(:, 2)./data_equal.mi_h_f(:, 2), '-.r', LineWidth=1.5);
hold on;
h3 = semilogx(d_range, data_equal.mi_st_f(:, 3)./data_equal.mi_h_f(:, 3), 'r', LineWidth=1.5);
hold on;
h4 = semilogx(d_range, data_equal.mi_st_t(:, 1)./data_equal.mi_h_t(:, 1), '--k', LineWidth=1.5);
hold on;
h5 = semilogx(d_range, data_equal.mi_st_t(:, 2)./data_equal.mi_h_t(:, 2), '-.k', LineWidth=1.5);
hold on;
h6 = semilogx(d_range, data_equal.mi_st_t(:, 3)./data_equal.mi_h_t(:, 3), 'k', LineWidth=1.5);
hold off;
% legend('LBKG, {\it{N}}_{path}=1','LBKG, {\it{N}}_{path}=2', 'LBKG, {\it{N}}_{path}=10', 'LB-TDD, {\it{N}}_{path}=1','LB-TDD, {\it{N}}_{path}=2', 'LB-TDD, {\it{N}}_{path}=10')
% xlabel('Eavesdropper Separation {\it{d}} (wavelengths)');
% ylabel('I_{VK}/I_{K}');
% ylim([0,0.9]);
xlabel('Eavesdropper Separation {\it{w}} (wavelengths)');
ylabel('I_{VK}/I_{K}');
ylim([0,0.9]);
set(gca,'FontSize',18);
lgd1 = legend([h1, h2, h3],'LBKG, {\it{N}}_{path}=1','LBKG, {\it{N}}_{path}=2', 'LBKG, {\it{N}}_{path}=10');
set(lgd1,'FontSize',14);
hold on
ah=axes('position',get(gca,'position'),'visible','off');
lgd2 = legend(ah,[h4, h5, h6],'LB-TDD, {\it{N}}_{path}=1','LB-TDD, {\it{N}}_{path}=2', 'LB-TDD, {\it{N}}_{path}=10');
set(lgd2,'FontSize',14);

error('OK')

figure
for ii=1:length(Npath_range)
  semilogx(d_range, data_equal.mi_h_f(:, ii), 'r');
  hold on;
end
for ii=1:length(Npath_range)
  semilogx(d_range, data_equal.mi_h_t(:, ii), 'k');
  hold on;
end
hold off;
xlabel('Eavesdropper Separation w (wavelengths)');
ylabel('I_{K}');

fname = strcat(data_dir, '/sing_ant_vs_d.mat');

snr_arr = 0:30;
M_range = [3000 3000 3000];

mi_h_f = zeros(length(snr_arr), length(Npath_range));
mi_sa_f = mi_h_f;
mi_st_f = mi_h_f;
mi_h_t = zeros(length(snr_arr), length(Npath_range));
mi_sa_t = mi_h_t;
mi_st_t = mi_h_t;

rff_diff_range = [0,0.1,0.2];
for ii=1:length(rff_diff_range)
  for jj=1:length(snr_arr)
   
    if (do_show_progress)
      fprintf('[%d,%d] ', ii, jj);
      if (do_flush), fflush(1); end
    end
    
    clear ops
    ops.SNR1 = snr_arr(jj);
    ops.SNR2 = snr_arr(jj);
    ops.SNR3 = snr_arr(jj);
    ops.N1 = 1;
    ops.N2 = 1;
    ops.N3 = 1;
%     ops.d = d_range(jj);
    ops.d  = 1;
%     ops.Nray = Npath_range(ii);
    ops.Nray = 2;
    ops.diff = rff_diff_range(ii);
    ops.M = M_range(ii);
    [mi_h_f(jj,ii), mi_sa_f(jj,ii), mi_st_f(jj,ii), mi_h_t(jj,ii), mi_sa_t(jj,ii), mi_st_t(jj,ii)] = do_mi1_avg(ops);
  end
end

data_equal.mi_h_f = mi_h_f;
data_equal.mi_sa_f = mi_sa_f;
data_equal.mi_st_f = mi_st_f;

data_equal.mi_h_t = mi_h_t;
data_equal.mi_sa_t = mi_sa_t;
data_equal.mi_st_t = mi_st_t;

figure
plot(snr_arr, data_equal.mi_st_f(:, 1)./data_equal.mi_h_f(:, 1), '--r', LineWidth=1);
hold on;
plot(snr_arr, data_equal.mi_st_f(:, 2)./data_equal.mi_h_f(:, 2), '-.r', LineWidth=1);
hold on;
plot(snr_arr, data_equal.mi_st_f(:, 3)./data_equal.mi_h_f(:, 3), 'r', LineWidth=1);
hold on;
plot(snr_arr, data_equal.mi_st_t(:, 1)./data_equal.mi_h_t(:, 1), '--k', LineWidth=1);
hold on;
plot(snr_arr, data_equal.mi_st_t(:, 2)./data_equal.mi_h_t(:, 2), '-.k', LineWidth=1);
hold on;
plot(snr_arr, data_equal.mi_st_t(:, 3)./data_equal.mi_h_t(:, 3), 'k', LineWidth=1);
hold off;
xlabel('snr');
ylabel('I_{VK}/I_{K}');
set(gca,'FontSize',24);
legend('FDD, {\it{d}}=0','FDD, {\it{d}}=0.1', 'FDD, {\it{d}}=0.2', 'TDD, {\it{d}}=0','TDD, {\it{d}}=0.1', 'TDD, {\it{d}}=0.2', 'Fontsize', 22);


figure
plot(snr_arr, data_equal.mi_h_f(:, 1), '--r', LineWidth=1);
hold on;
plot(snr_arr, data_equal.mi_h_f(:, 2), '-.r', LineWidth=1);
hold on;
plot(snr_arr, data_equal.mi_h_f(:, 3), 'r', LineWidth=1);
hold on;
plot(snr_arr, data_equal.mi_h_t(:, 1), '--k', LineWidth=1);
hold on;
plot(snr_arr, data_equal.mi_h_t(:, 2), '-.k', LineWidth=1);
hold on;
plot(snr_arr, data_equal.mi_h_t(:, 3), 'k', LineWidth=1);
hold off;
legend('FDD, {\it{d}}=0','FDD, {\it{d}}=0.1', 'FDD, {\it{d}}=0.2', 'TDD, {\it{d}}=0','TDD, {\it{d}}=0.1', 'TDD, {\it{d}}=0.2')
xlabel('SNR');
ylabel('I_{K}');

figure
plot(snr_arr, data_equal.mi_sa_f(:, 1), '--r', LineWidth=1.5);
hold on;
plot(snr_arr, data_equal.mi_sa_f(:, 2), '-.r', LineWidth=1.5);
hold on;
plot(snr_arr, data_equal.mi_sa_f(:, 3), 'r', LineWidth=1.5);
hold on;
plot(snr_arr, data_equal.mi_sa_t(:, 1), '--k', LineWidth=1.5);
hold on;
plot(snr_arr, data_equal.mi_sa_t(:, 2), '-.k', LineWidth=1.5);
hold on;
plot(snr_arr, data_equal.mi_sa_t(:, 3), 'k', LineWidth=1.5);
hold off;
xlabel('SNR');
ylabel('I_{sK}');
set(gca,'FontSize',24);

% figure
% for ii=1:length(Npath_range)
%   semilogx(snr_arr, data_equal.mi_h_f(:, ii), 'r');
%   hold on;
% end
% for ii=1:length(Npath_range)
%   semilogx(snr_arr, data_equal.mi_h_t(:, ii), 'k');
%   hold on;
% end
% hold off;
% xlabel('Eavesdropper Separation d (wavelengths)');
% ylabel('I_{K}');

toc
error('hh')

% SNR advantage case %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = data_equal.mi_st./data_equal.mi_h;
d_range = d_range.';

SNR = 15;
SNRe = 35;
N = 1;
Npath_range = [1 2 10];
d_max = 10;
d_min = 0.1;
d_N = 100;
d_range = 10.^((log10(d_max)-log10(d_min))/(d_N-1)*(0:d_N-1) + log10(d_min));
M_range = [10 5000 1000];

fname = strcat(data_dir, '/sing_ant_vs_d_snra.mat');
if ~exist(fname)

mi_h_f = zeros(d_N, length(Npath_range));
mi_sa_f = mi_h_f;
mi_st_f = mi_h_f;

for ii=1:length(Npath_range)
  for jj=1:length(d_range)
   
    if (do_show_progress)
      fprintf('[%d,%d] ', ii, jj);
      if (do_flush), fflush(1); end
    end
    
    clear ops
    ops.SNR1 = SNR;
    ops.SNR2 = SNR;
    ops.SNR3 = SNRe;
    ops.N1 = 1;
    ops.N2 = 1;
    ops.N3 = 1;
    ops.d = d_range(jj);
    ops.Nray = Npath_range(ii);
    ops.M = M_range(ii);
    [mi_h_f(jj,ii), mi_sa_f(jj,ii), mi_st_f(jj,ii)] = do_mi1_avg(ops);
  end
end

% save(fname, 'mi_h', 'mi_sa', 'mi_st');

else
  
% load(fname);

end

data_adv.mi_h = mi_h_f;
data_adv.mi_sa = mi_sa_f;
data_adv.mi_st = mi_st_f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(3);

for ii=1:length(Npath_range)
  semilogx(d_range, data_equal.mi_st(:, ii)./data_equal.mi_h(:, ii), 'k');
  hold on;
  semilogx(d_range, data_adv.mi_st(:, ii)./data_adv.mi_h(:, ii), 'b--');
end
hold off;
xlabel('Eavesdropper Separation d (wavelengths)');
ylabel('I_{VK}/I_{K}');
ylim([0 1.1]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 4: Four antenna (Eve: four antennas (equal) or 10 antennas)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data_equal = [];
data_adv = [];

%% Four antenna, no advantage %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SNR = 15;
N = 4;
Npath_range = [1 10 100];
d_max = 100;
d_min = 0.1;
d_N = 100;
d_range = 10.^((log10(d_max)-log10(d_min))/(d_N-1)*(0:d_N-1) + log10(d_min));
M_range = [10 5000 1000];

fname = strcat(data_dir, '/four_ant_vs_d.dat');
if ~exist(fname)

mi_h_f = zeros(d_N, length(Npath_range));
mi_sa_f = mi_h_f;
mi_st_f = mi_h_f;

for ii=1:length(Npath_range)
  for jj=1:length(d_range)

    if (do_show_progress)
      fprintf('[%d,%d] ', ii, jj);
      if (do_flush), fflush(1); end
    end
    
    clear ops
    ops.SNR1 = SNR;
    ops.SNR2 = SNR;
    ops.SNR3 = SNR;
    ops.N1 = N;
    ops.N2 = N;
    ops.N3 = N;
    ops.d = d_range(jj);
    ops.Nray = Npath_range(ii);
    ops.M = M_range(ii);
    
    [mi_h_f(jj,ii), mi_sa_f(jj,ii), mi_st_f(jj,ii)] = do_mi1_avg(ops);
  end
end

save(fname, 'mi_h_f', 'mi_sa_f', 'mi_st_f');

else
  
load(fname);

end

data_equal.mi_h = mi_h_f;
data_equal.mi_sa = mi_sa_f;
data_equal.mi_st = mi_st_f;

%% Four antenna, Antenna advantage at Eve %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SNR = 15;
N = 4;
Ne = 10;
Npath_range = [1 10 100];
d_max = 100;
d_min = 0.1;
d_N = 100;
d_range = 10.^((log10(d_max)-log10(d_min))/(d_N-1)*(0:d_N-1) + log10(d_min));
M_range = [10 5000 2000];

fname = strcat(data_dir, '/four_ant_vs_d_edant.dat');
if ~exist(fname)

mi_h_f = zeros(d_N, length(Npath_range));
mi_sa_f = mi_h_f;
mi_st_f = mi_h_f;

for ii=1:length(Npath_range)
  for jj=1:length(d_range)

    if (do_show_progress)
      fprintf('[%d,%d] ', ii, jj);
      if (do_flush), fflush(1); end
    end

    clear ops
    ops.SNR1 = SNR;
    ops.SNR2 = SNR;
    ops.SNR3 = SNR;
    ops.N1 = N;
    ops.N2 = N;
    ops.N3 = Ne;
    ops.d = d_range(jj);
    ops.Nray = Npath_range(ii);
    ops.M = M_range(ii);
    [mi_h_f(jj,ii), mi_sa_f(jj,ii), mi_st_f(jj,ii)] = do_mi1_avg(ops);
  end
end

save(fname, 'mi_h_f', 'mi_sa_f', 'mi_st_f');

else
  
load(fname);

end

data_adv.mi_h = mi_h_f;
data_adv.mi_sa = mi_sa_f;
data_adv.mi_st = mi_st_f;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(4);

for ii=1:length(Npath_range)
  semilogx(d_range, data_equal.mi_st(:, ii)./data_equal.mi_h(:, ii), 'k');
  hold on;
  semilogx(d_range, data_adv.mi_st(:, ii)./data_adv.mi_h(:, ii), 'b--');
end
hold off;
xlabel('Eavesdropper Separation d (wavelengths)');
ylabel('I_{VK}/I_{K}');
ylim([0 1.1]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 6: CQA/CQG symbol error rate performance.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



