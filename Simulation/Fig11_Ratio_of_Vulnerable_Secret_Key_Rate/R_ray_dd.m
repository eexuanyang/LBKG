function [R_f, R_t] = R_ray_dd(xt, yt, xr, yr, beta, th_t, th_r, rff_diff)

% [R] = R_ray_dd(xt, yt, xr, yr, beta, th_t, th_r);
%
% Computes the double-directional covariance matrix for
% a set of rays.
%
% Inputs:
%	xt, yt		Transmit antenna positions
%	xr, yr		Receive antenna positions
%	beta		Ray powers
%	th_t, th_r	Ray departure and arrival angles
%
% Outputs:
%	R(i,j;k,l)	Covariance matrix

%%L = length(beta);

% Get stacked antenna positions
Nr = length(xr);
Nt = length(xt);

xt = xt(:).';
yt = yt(:).';
ot = ones(size(xt));

xr = xr(:).';
yr = yr(:).';
or = ones(size(xr));

xr1 = kron(kron(ot,or), kron(ot,xr));
xr2 = kron(kron(ot,xr), kron(ot,or));
yr1 = kron(kron(ot,or), kron(ot,yr));
yr2 = kron(kron(ot,yr), kron(ot,or));

xt1 = kron(kron(ot,or), kron(xt,or));
xt2 = kron(kron(xt,or), kron(ot,or));
yt1 = kron(kron(ot,or), kron(yt,or));
yt2 = kron(kron(yt,or), kron(ot,or));

B = beta(:)*ones(1, Nr*Nt*Nt*Nr);

th_t = th_t(:);
th_r = th_r(:);

if (length(beta)==1)
  R_f = reshape((B.*exp(1i*2*pi*(cos(th_r)*(xr1-xr2) + sin(th_r)*(yr1-yr2) + cos(th_t)*(xt1-xt2) + sin(th_t)*(yt1-yt2)))), Nr*Nt, Nr*Nt);
  R_t = reshape((B.*exp(1i*2*pi*(cos(th_r)*(xr1-xr2) + sin(th_r)*(yr1-yr2) + cos(th_t)*(xt1-xt2) + sin(th_t)*(yt1-yt2)))), Nr*Nt, Nr*Nt);
else
  R_f = reshape(sum(B.*exp(1i*2*pi*(cos(th_r)*(xr1-xr2) + sin(th_r)*(yr1-yr2) + cos(th_t)*(xt1-xt2) + sin(th_t)*(yt1-yt2)))), Nr*Nt, Nr*Nt);
  R_t = reshape(sum(B.*exp(1i*2*pi*(cos(th_r)*(xr1-xr2) + sin(th_r)*(yr1-yr2) + cos(th_t)*(xt1-xt2) + sin(th_t)*(yt1-yt2)))), Nr*Nt, Nr*Nt);
end

% th_r1 = th_r;
% th_t1 = th_t;
th_r1 = 2*pi*rand(length(beta), 1);
th_t1 = 2*pi*rand(length(beta), 1);
th_r2 = 2*pi*rand(length(beta), 1);
th_t2 = 2*pi*rand(length(beta), 1);

Ga1 = 1;
Gb1 = 1 - rff_diff;
Ga2 = 1 - rff_diff;
Gb2 = 1;
c = 1;

% Ga1 = 1 - rff_diff;
% Gb1 = 1;
% Ga2 = 1;
% Gb2 = 1 - rff_diff;
% c = 1;

% Ga1 = 1; Gb1 = 1; Ga2 = 1; Gb2 = 1;
% % Gb1 = normrnd(0.9,0.1,1,1);
% % Ga2 = normrnd(0.9,0.1,1,1);
r=1;
h1_temp = zeros(length(beta),1);
for idx_k = 1:length(beta)
    h1_temp(idx_k) = sqrt(r)*exp(1i*2*pi*(cos(th_r1(idx_k))*xr1(1)+sin(th_r1(idx_k))*yr1(1)+cos(th_t1(idx_k))*xt1(1)+sin(th_t1(idx_k))*yt1(1)));
end
% hab1 = Ga1*h1_temp;
% hba1 = Gb1*h1_temp;

% % he_temp = zeros(length(beta),1);
% % for idx_k = 1:length(beta)
% %     he_temp(idx_k) = sqrt(r)*exp(1i*2*pi*(cos(th_r1(idx_k))*xr2(1)+sin(th_r1(idx_k))*yr2(1)+cos(th_t1(idx_k))*xt2(1)+sin(th_t1(idx_k))*yt2(4)));
% % end
% hc = he_temp; 
% R = [sum(B(1)*hab1.*conj(hab1)) sum(B(1)*hab1.*conj(hba1)) sum(B(1)*hab1.*conj(hc)) sum(B(1)*hba1.*conj(hba1)) sum(B(1)*hba1.*conj(hc)) sum(B(1)*hc.*conj(hc)) R(:).']

h2_temp = zeros(length(beta),1);
for idx_k = 1:length(beta)
    h2_temp(idx_k) = sqrt(r)*exp(1i*2*pi*(cos(th_r2(idx_k))*xr1(1)+sin(th_r2(idx_k))*yr1(1)+cos(th_t2(idx_k))*xt1(1)+sin(th_t2(idx_k))*yt1(1)));
end
hab2 = Ga2*Gb1*h1_temp.*h2_temp;
hba2 = Gb2*Ga1*h1_temp.*h2_temp;

he_temp = zeros(length(beta),1);
for idx_k = 1:length(beta)
    he_temp(idx_k) = sqrt(r)*exp(1i*2*pi*(cos(th_r2(idx_k))*xr2(1)+sin(th_r2(idx_k))*yr2(1)+cos(th_t2(idx_k))*xt2(1)+sin(th_t2(idx_k))*yt2(4)));
end
hc = c*Ga1*h1_temp.*he_temp;
% R_t = [sum(B(1)*hab2.*conj(hab2)) sum(B(1)*hab2.*conj(hba2)) sum(B(1)*hab2.*conj(hc)) sum(B(1)*hba2.*conj(hba2)) sum(B(1)*hba2.*conj(hc)) sum(B(1)*hc.*conj(hc)) 1+Gb2^2 1+Ga2^2 1+Gb2^2];
% R_t = [sum(B(1)*hab2.*conj(hab2)) sum(B(1)*hab2.*conj(hba2)) sum(B(1)*hab2.*conj(hc)) sum(B(1)*hba2.*conj(hba2)) sum(B(1)*hba2.*conj(hc)) sum(B(1)*hc.*conj(hc)) 1 1 1];
R_t = [sum(B(1)*hab2.*conj(hab2)) sum(B(1)*hab2.*conj(hba2)) sum(B(1)*hab2.*conj(hc)) sum(B(1)*hba2.*conj(hba2)) sum(B(1)*hba2.*conj(hc)) sum(B(1)*hc.*conj(hc)) 1+Gb2^2 1+Ga2^2 1+c^2];

hab4 = Ga1*Ga2*Gb1*Gb2*h1_temp.*h2_temp.*h1_temp.*h2_temp;
hba4 = Ga1*Ga2*Gb1*Gb2*h1_temp.*h2_temp.*h1_temp.*h2_temp;
hc   = Ga1*Ga2*Gb1* c *h1_temp.*h2_temp.*h1_temp.*he_temp;
% R_f = [sum(B(1)*hab4.*conj(hab4)) sum(B(1)*hab4.*conj(hba4)) sum(B(1)*hab4.*conj(hc)) sum(B(1)*hba4.*conj(hba4)) sum(B(1)*hba4.*conj(hc)) sum(B(1)*hc.*conj(hc))...
%     (Gb1*Ga2*Gb2)^2+(Ga2*Gb2)^2+Gb2^2+1     (Ga1*Gb1*Ga2)^2+(Ga1*Gb1)^2+Ga1^2+1     (Gb1*Ga2*Gb2)^2+(Ga2*Gb2)^2+Gb2^2+1];
% R_f = [sum(B(1)*hab4.*conj(hab4)) sum(B(1)*hab4.*conj(hba4)) sum(B(1)*hab4.*conj(hc)) sum(B(1)*hba4.*conj(hba4)) sum(B(1)*hba4.*conj(hc)) sum(B(1)*hc.*conj(hc)) 1 1 1];
R_f = [sum(B(1)*hab4.*conj(hab4)) sum(B(1)*hab4.*conj(hba4)) sum(B(1)*hab4.*conj(hc)) sum(B(1)*hba4.*conj(hba4)) sum(B(1)*hba4.*conj(hc)) sum(B(1)*hc.*conj(hc)) ...
     2*(Gb1*Ga2*Gb2)^2+2*(Ga2*Gb2)^2+Gb2^2+1     2*(Ga1*Gb1*Ga2)^2+2*(Ga1*Gb1)^2+Ga1^2+1     (Gb1*Ga2*c)^2+(Ga2*c)^2+c^2+1];

% % % a=1;
% % 
% % h2_temp = 0;
% % for idx_k = 1:length(beta)
% %     h2_temp = h2_temp + B(idx_k)*exp(1i*2*pi*(cos(th_r2(idx_k))*xr1(1)+sin(th_r2(idx_k))*yr1(1)+cos(th_t2(idx_k))*xt1(1)+sin(th_t2(idx_k))*yt1(1)));
% % end
% % hba2 = Gb2*hab1*h2_temp;
% % hab2 = Ga2*hba1*h2_temp;
% % % % % hab1 = Ga1*B(1)*exp(1i*2*pi*(cos(th_r)*xr1(1)+sin(th_r)*yr1(1)+cos(th_t)*xt1(1)+sin(th_t)*yt1(1)));
% % % % % hba1 = Gb1*B(1)*exp(1i*2*pi*(cos(th_r)*xr1(1)+sin(th_r)*yr1(1)+cos(th_t)*xt1(1)+sin(th_t)*yt1(1)));
% % % % % hba2 = Gb2*hab1*B(1)*exp(1i*2*pi*(cos(th_r2)*xr1(1)+sin(th_r2)*yr1(1)+cos(th_t2)*xt1(1)+sin(th_t2)*yt1(1)));
% % % % % hab2 = Ga2*hba1*B(1)*exp(1i*2*pi*(cos(th_r2)*xr1(1)+sin(th_r2)*yr1(1)+cos(th_t2)*xt1(1)+sin(th_t2)*yt1(1)));
% % he_temp = 0;
% % for idx_k = 1:length(beta)
% %     he_temp = he_temp + B(idx_k)*exp(1i*2*pi*(cos(th_r2(idx_k))*xr2(1)+sin(th_r2(idx_k))*yr2(1)+cos(th_t2(idx_k))*xt2(1)+sin(th_t2(idx_k))*yt2(4)));
% % end
% % % hc = hba2*B(1)*exp(1i*2*pi*(cos(th_r2(idx_k))*xr2(1)+sin(th_r2(idx_k))*yr2(1)+cos(th_t2(idx_k))*xt2(1)+sin(th_t2(idx_k))*yt2(4)));
% % hc = Gb2*hab1*he_temp;
% % 
% % R = [hab2*hab2' hab2*hba2' hab2*hc' hba2*hba2' hba2*hc' hc*hc'];
