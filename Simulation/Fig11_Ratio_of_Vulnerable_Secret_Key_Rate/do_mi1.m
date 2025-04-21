function [mi_h_f, mi_sa_f, mi_st_f, mi_h_t, mi_sa_t, mi_st_t] = do_mi1(ops)

% [mi_h mi_sa mi_st] = do_mi1(ops);
%
% Computes the mutual information for the varying channel and
% how much of this information is safe.  
%
% ops:
%	Nray		Number of rays in model
%	SNR1..3		SNRs of nodes in dB
%	N1..3		Number of antennas of users
%	d		Distance from eavesdropper to legit receiver (wl)
%	info		Print/plot info if =1

%% Model type
%mod_type = 0;	% SVA model
mod_type = 1;	% Fixed number of rays (iid arrivals)
%Nray = 10;	% Total number of rays, or rays per cluster
Nray = ops.Nray;

%% Estimation Error (prop to SNR)
SNR1 = ops.SNR1;
SNR2 = ops.SNR2;
SNR3 = ops.SNR3;

if isfield(ops, 'info')
  info = ops.info;
else
  info = 0;
end

% Distance from eavesdropper to legit user
d13 = ops.d;
%d13 = 0.5;
%d13 = 100;

%% SVA Model parameters

% Cluster decay rate
gamma = 1;

% Cluster angular spread
sigma = 25*pi/180;

% Basis type (0=gaussian, 1=laplacian)
basis = 1;

% Threshold to stop looking for clusters
beta_sq_th = 0.01;

%% Antenna Arrays

% Alice
N1 = ops.N1;
d1 = 0.5;
x1 = (0:N1-1)*d1 - (N1-1)*d1/2;
y1 = zeros(size(x1));

% Bob
N2 = ops.N2;
%d2 = 0.5;
d2 = 0.5;
x2 = (0:N2-1)*d2 - (N2-1)*d2/2;
y2 = zeros(size(x2));
x2 = 100; y2 = 100;

% Eve
N3 = ops.N3;
d3 = 0.5;
x3 = (0:N3-1)*d3 - (N3-1)*d3/2;
y3 = -ones(size(x3))*d13;  % Put below Alice

% Get combined Alice/Eve channel
x13 = [x1(:); x3(:)];
y13 = [y1(:); y3(:)];

% Get linear power
r = 1;
snr1 = 10^(SNR1/10); sigma1 = r/snr1;
snr2 = 10^(SNR2/10); sigma2 = r/snr2;
snr3 = 10^(SNR3/10); sigma3 = r/snr3;

if (mod_type == 0)
  % SVA model

  % Get clusters
  [phi_t0, phi_r0, beta0] = sva_clus(gamma, beta_sq_th);

  % Get rays in clusters
  [beta, phi_t, phi_r] = sva_ray(phi_t0, phi_r0, beta0, Nray, sigma, basis);
  beta = abs(beta).^2;
else
  % Simple ray model
  beta = ones(Nray, 1);
  phi_t = 2*pi*rand(Nray, 1); phi_t0 = phi_t;
  phi_r = 2*pi*rand(Nray, 1); phi_r0 = phi_r;
end
  
% Normalize angles to be between [-pi,pi]
phi_t = angle(exp(1i*phi_t));
phi_r = angle(exp(1i*phi_r));
phi_t0 = angle(exp(1i*phi_t0));
phi_r0 = angle(exp(1i*phi_r0));

% Make the channel unit power (on average)
beta = beta/sum(beta);

% Plot as sanity check
if (info)
figure(1);
beta_db = 10*log10(beta(:));
data1 = [phi_t(:) phi_r(:) beta_db(:)];
data2 = [phi_t0(:) phi_r0(:)];
eval(sprintf('gset zrange [%f:%f]', max(beta_db-40), max(beta_db)));
gset pointsize 3;
%gset autoscale x
%gset autoscale y
gset autoscale z;
gset nokey;
gset xrange [-pi:pi];
gset yrange [-pi:pi];
%gplot data1 with points3d 1, data2 with points 2;
end

% Compute the covariance for legitimate channel
%r12 = R_ray_dd(x2, y2, x1, y1, abs(beta(:)).^2, phi_t(:), phi_r(:));

% Compute the covariance for far eavesdropper channel
%r32 = R_ray_dd(x3, y3, x1, y1, abs(beta(:)).^2, phi_t(:), phi_r(:));

% Get the combined covariance of 1-2 and 2-3
%R_AC = R_ray_dd(x2, y2, x13, y13, abs(beta(:)).^2, phi_t(:), phi_r(:));
[R_AC_f, R_AC_t]= R_ray_dd(x13, y13, x2, y2, beta(:), phi_r(:), phi_t(:), ops.diff);

if (info)
figure(2);
gset zrange [0:1];
plotimgg(1:columns(R_AC_f), 1:rows(R_AC_f), abs(R_AC_f));
end

% % % Cut out individual parts
% % R_aa = R_AC(7);
% % R_cc = R_AC(10);
% % R_ac = R_AC(9);
% % 
% % % Compute mutual information of legit channel
% % Ia = eye(N1*N2);
% % Ic = eye(N3*N2);
% % Rh_aa = (R_aa + Ia*sigma2);
% % Rh_AA = (R_aa + Ia*sigma1);
% % Rh_aA = [R_aa + Ia*sigma2,	R_aa;
% % 	 R_aa			R_aa+Ia*sigma1];
% % I12 = log2(det(Rh_aa)*det(Rh_AA)/det(Rh_aA));
% % 
% % % Compute mutual information that is secret
% % Rh_ac = [R_aa + Ia*sigma2,	R_ac;
% % 	 R_ac',			R_cc + Ic*sigma3];
% % Rh_Ac = [R_aa + Ia*sigma1,	R_ac;
% % 	 R_ac',			R_cc + Ic*sigma3];
% % Rh_cc = (R_cc + Ic*sigma3);
% % Rh_aAc = [R_aa + Ia*sigma2,	R_aa,			R_ac;
% % 	  R_aa,			R_aa + Ia*sigma1,	R_ac;
% % 	  R_ac',		R_ac',			R_cc + Ic*sigma3];
% % I12s = log2(det(Rh_ac)*det(Rh_Ac)/det(Rh_cc)/det(Rh_aAc));
% % 
% % % Get stealible mutal information (vulnerable key bits)
% % I3 = I12-I12s;

% FDD
sigma1 = sigma1*R_AC_f(7);
sigma2 = sigma2*R_AC_f(8);
sigma3 = sigma3*R_AC_f(9);

% Cut out individual parts
R_aa = R_AC_f(1);
R_aa1= R_AC_f(2);
R_ac = R_AC_f(3);
R_a1a1 = R_AC_f(4);
R_a1c  = R_AC_f(5);
R_cc = R_AC_f(6);

% Compute mutual information of legit channel
Ia = eye(N1*N2);
Ic = eye(N3*N2);
% Rh_aa = (R_aa + Ia*sigma2);
% Rh_AA = (R_a1a1 + Ia*sigma1);
% Rh_aA = [R_aa + Ia*sigma2,	R_aa1;
% 	 R_aa1'			R_a1a1+Ia*sigma1];
% I12 = log2(det(Rh_aa)*det(Rh_AA)/det(Rh_aA));
Rh_aa = (R_aa + Ia*sigma1);
Rh_AA = (R_a1a1 + Ia*sigma2);
Rh_aA = [R_aa + Ia*sigma1,	R_aa1;
	 R_aa1'			R_a1a1+Ia*sigma2];
I12 = log2(det(Rh_aa)*det(Rh_AA)/det(Rh_aA));

% Compute mutual information that is secret
% Rh_ac = [R_aa + Ia*sigma2,	R_ac;
% 	 R_ac',			R_cc + Ic*sigma3];
% Rh_Ac = [R_a1a1 + Ia*sigma1,	R_a1c;
% 	 R_a1c',			R_cc + Ic*sigma3];
% Rh_cc = (R_cc + Ic*sigma3);
% Rh_aAc = [R_aa + Ia*sigma2,	R_aa1,			R_ac;
% 	  R_aa1',			R_a1a1 + Ia*sigma1,	R_a1c;
% 	  R_ac',		R_a1c',			R_cc + Ic*sigma3];
% I12s = log2(det(Rh_ac)*det(Rh_Ac)/det(Rh_cc)/det(Rh_aAc));
Rh_ac = [R_aa + Ia*sigma1,	R_ac;
	 R_ac',			R_cc + Ic*sigma3];
Rh_Ac = [R_a1a1 + Ia*sigma2,	R_a1c;
	 R_a1c',			R_cc + Ic*sigma3];
Rh_cc = (R_cc + Ic*sigma3);
Rh_aAc = [R_aa + Ia*sigma1,	R_aa1,			R_ac;
	  R_aa1',			R_a1a1 + Ia*sigma2,	R_a1c;
	  R_ac',		R_a1c',			R_cc + Ic*sigma3];
I12s = log2(det(Rh_ac)*det(Rh_Ac)/det(Rh_cc)/det(Rh_aAc));

% Get stealible mutal information (vulnerable key bits)
I3 = I12-I12s;
 
if (info)
fprintf('Shared Information:	%.2f bits\n', real(I12));
fprintf('Safe Information:	%.2f bits\n', real(I12s));
fprintf('Unsafe Information:	%.2f bits\n', real(I3));
end

mi_h_f = real(I12);
mi_sa_f = real(I12s);
mi_st_f = real(I3);

% TDD
sigma1 = sigma1*R_AC_t(7);
sigma2 = sigma2*R_AC_t(8);
sigma3 = sigma3*R_AC_t(9);

% Cut out individual parts
R_aa = R_AC_t(1);
R_aa1= R_AC_t(2);
R_ac = R_AC_t(3);
R_a1a1 = R_AC_t(4);
R_a1c  = R_AC_t(5);
R_cc = R_AC_t(6);

% Compute mutual information of legit channel
Ia = eye(N1*N2);
Ic = eye(N3*N2);
% Rh_aa = (R_aa + Ia*sigma2);
% Rh_AA = (R_a1a1 + Ia*sigma1);
% Rh_aA = [R_aa + Ia*sigma2,	R_aa1;
% 	 R_aa1'			R_a1a1+Ia*sigma1];
% I12 = log2(det(Rh_aa)*det(Rh_AA)/det(Rh_aA));
Rh_aa = (R_aa + Ia*sigma1);
Rh_AA = (R_a1a1 + Ia*sigma2);
Rh_aA = [R_aa + Ia*sigma1,	R_aa1;
	 R_aa1'			R_a1a1+Ia*sigma2];
I12 = log2(det(Rh_aa)*det(Rh_AA)/det(Rh_aA));

% Compute mutual information that is secret
% Rh_ac = [R_aa + Ia*sigma2,	R_ac;
% 	 R_ac',			R_cc + Ic*sigma3];
% Rh_Ac = [R_a1a1 + Ia*sigma1,	R_a1c;
% 	 R_a1c',			R_cc + Ic*sigma3];
% Rh_cc = (R_cc + Ic*sigma3);
% Rh_aAc = [R_aa + Ia*sigma2,	R_aa1,			R_ac;
% 	  R_aa1',			R_a1a1 + Ia*sigma1,	R_a1c;
% 	  R_ac',		R_a1c',			R_cc + Ic*sigma3];
% I12s = log2(det(Rh_ac)*det(Rh_Ac)/det(Rh_cc)/det(Rh_aAc));
Rh_ac = [R_aa + Ia*sigma1,	R_ac;
	 R_ac',			R_cc + Ic*sigma3];
Rh_Ac = [R_a1a1 + Ia*sigma2,	R_a1c;
	 R_a1c',			R_cc + Ic*sigma3];
Rh_cc = (R_cc + Ic*sigma3);
Rh_aAc = [R_aa + Ia*sigma1,	R_aa1,			R_ac;
	  R_aa1',			R_a1a1 + Ia*sigma2,	R_a1c;
	  R_ac',		R_a1c',			R_cc + Ic*sigma3];
I12s = log2(det(Rh_ac)*det(Rh_Ac)/det(Rh_cc)/det(Rh_aAc));

% Get stealible mutal information (vulnerable key bits)
I3 = I12-I12s;
 
if (info)
fprintf('Shared Information:	%.2f bits\n', real(I12));
fprintf('Safe Information:	%.2f bits\n', real(I12s));
fprintf('Unsafe Information:	%.2f bits\n', real(I3));
end

mi_h_t = real(I12);
mi_sa_t = real(I12s);
mi_st_t = real(I3);




