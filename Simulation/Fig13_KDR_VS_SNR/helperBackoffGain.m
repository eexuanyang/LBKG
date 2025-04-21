function gain = helperBackoffGain(ipPwrAtPkOut, PxdBm, IBO)
% HELPERBACKOFFGAIN Set the gain to back off a signal input to an amplifier
%  This function is only in support of
%  SimulateAndVerifyBackoffForAPowerAmplifierExample. It may change in a
%  future release. 
%
% GAIN = HELPERBACKOFFGAIN(IPPWRATPKOUT, PXDBM, IBO) calculates the gain
% needed to scale a signal to achieve a desired input backoff prior to
% inputting the signal to a nonlinear amplifier. IPPWRATPKOUT is the input
% power, in dBm, that results in the peak output power of the amplifier.
% PXDBM is the power, in dBm, of the signal.  IBO, the input backoff in dB,
% is the power reduction relative to IPPWRATPKOUT required to scale the
% signal.
%

gaindB = ipPwrAtPkOut - PxdBm - IBO;  % dB power gain
gain = 10^(gaindB/20);  % linear voltage gain

% EOF
