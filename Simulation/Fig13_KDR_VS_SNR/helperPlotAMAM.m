function hFig = helperPlotAMAM(mnl)
%HELPERPLOTAMAM Plot the AM/AM characteristic of a table-based amplifier
%  This function is only in support of
%  SimulateAndVerifyBackoffForAPowerAmplifierExample. It may change in a
%  future release. 
%
%  HFIG = HELPERPLOTAMAM(MNL) plots the AM/AM characteristic of a
%  table-based amplifier. MNL is a comm.MemorylessNonlinearity System
%  object. HFIG is a handle to the figure.
%

% Copyright 2020 The MathWorks, Inc.

% Calculate power in and power out
powXmin = (floor(mnl.Table(1,1)/10))*10;
powXmax = (ceil(mnl.Table(end,1)/10))*10;
pIndBm = [powXmin; mnl.Table(:,1); powXmax];
vin = 10.^((pIndBm-30)/20);
vout = mnl(vin);
pOutdBm = 20*log10(abs(vout))+30;

% AM/AM response
hLine = plot(pIndBm, pOutdBm, "-b", ...
  pIndBm(2:end-1), pOutdBm(2:end-1), "*b", ...
  'LineWidth', 1.6);
xlabel("Instantaneous P_i_n (dBm)"); ylabel("Instantaneous P_o_u_t (dBm)");
xlim([-30 10]);
title("Specified and Actual AM/AM");
grid on;
hFig = get(get(hLine(1), "Parent"), "Parent");

% EOF
