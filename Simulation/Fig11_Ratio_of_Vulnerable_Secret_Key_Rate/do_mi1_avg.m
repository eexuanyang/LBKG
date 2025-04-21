function [mi_h_f, mi_sa_f, mi_st_f, mi_h_t, mi_sa_t, mi_st_t] = do_mi1_avg(ops)

% [mi_h mi_sa mi_st] = do_mi1_avg(ops);
%
% Computes the mutual information for the varying channel and
% how much of this information is safe.  
%
% ops:
%	M		Number of realizations
%	Nray		Number of rays in model
%	SNR1..3		SNRs of nodes in dB
%	N1..3		Number of antennas of users
%	d		Distance from eavesdropper to legit receiver (wl)
%	info		Print/plot info if =1

mi_h_f = 0;
mi_sa_f = 0;
mi_st_f = 0;

mi_h_t = 0;
mi_sa_t = 0;
mi_st_t = 0;

for m=1:ops.M
  [mi_h1_f,mi_sa1_f,mi_st1_f,mi_h1_t,mi_sa1_t,mi_st1_t] = do_mi1(ops);
  mi_h_f = mi_h_f + mi_h1_f;
  mi_sa_f = mi_sa_f + mi_sa1_f;
  mi_st_f = mi_st_f + mi_st1_f;
  mi_h_t = mi_h_t + mi_h1_t;
  mi_sa_t = mi_sa_t + mi_sa1_t;
  mi_st_t = mi_st_t + mi_st1_t;
end

mi_h_f = mi_h_f/ops.M;
mi_sa_f = mi_sa_f/ops.M;
mi_st_f = mi_st_f/ops.M;

mi_h_t = mi_h_t/ops.M;
mi_sa_t = mi_sa_t/ops.M;
mi_st_t = mi_st_t/ops.M;