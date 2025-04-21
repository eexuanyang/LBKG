% close all;
clear all;
clc
tic

% Parameters Setting
N = 8; % Number of filter taps
num_freq_points = 512; % 

codeword_len = 255; % actual number of bits per BCH codeword (code length)
infor_len = 131;
key_len_final = 128;

v = 5; % v = 18       % UT velocity in km/h
seg_num = 100;
mv_len = 200;
exp_time_max = 15000;
% exp_time_max = 100;

ds_ratio = 16;
snr_arr = 0:2:30;
% snr_arr = 20;

rff_h_a1b2 = ones(2048,1);
rff_h_a1b2 = ifft(rff_h_a1b2);
rff_h_a2b1 = 0.8*ones(2048,1);
rff_h_a2b1 = ifft(rff_h_a2b1);

%% Simulation Configuration
NFrames = 1;                % Number of frames to simulate at each SNR
% SNRIn = [-4.1, -2.0, 0.1];  % SNR points to simulate
SNRIn = 10;  % SNR points to simulate

%% UE Configuration
% User Equipment (UE) settings are specified in a structure form.
ue.TotSubframes = 1; % Total number of subframes to generate a waveform for
ue.NCellID = 10;     % Cell identity
ue.RC = 'A3-7';      % FRC number

mul_para = 16;
CP_array = [10 9 9 9 9 9 9 10 9 9 9 9 9 9] * mul_para;
fft_len = 128 * mul_para;
subframe_len = 1920*mul_para;
if strcmp(ue.RC, 'A3-5')
    start_subcar = 726;
    end_subcar   = 1325;
    subcar_len   = 600;
else
    start_subcar = 726 - 300;
    end_subcar   = 1325+ 300;
    subcar_len   = 1200;
end

%% Channel Estimator Configuration
cec.PilotAverage = 'UserDefined'; % Type of pilot averaging 
cec.FreqWindow = 13;              % Frequency averaging windows in REs
cec.TimeWindow = 1;               % Time averaging windows in REs
cec.InterpType = 'cubic';         % Interpolation type
cec.Reference = 'Antennas';       % Reference for channel estimation

%% Uplink RMC Configuration
% Generate FRC configuration structure for A3-2
frc = lteRMCUL(ue);

% Transport block sizes for each subframe within a frame
trBlkSizes = frc.PUSCH.TrBlkSizes;
codedTrBlkSizes = frc.PUSCH.CodedTrBlkSizes;
% 
% Redundancy version sequence
rvSequence = frc.PUSCH.RVSeq;  

%% Generate transmit signal
% Initialize state of all HARQ processes
harqProcesses = hNewHARQProcess(frc); 

harqProcessSequence = 1;

offsetused = 0;
subframeNo = 0;

% Update subframe number
frc.NSubframe = subframeNo;

% Get HARQ process ID for the subframe from HARQ process sequence
harqID = harqProcessSequence(mod(subframeNo, length(harqProcessSequence))+1);

% Update current HARQ process
harqProcesses(harqID) = hHARQScheduling(...
    harqProcesses(harqID), subframeNo, rvSequence);

% Update the PUSCH transmission config with HARQ process state
frc.PUSCH = harqProcesses(harqID).txConfig;
data = harqProcesses(harqID).data;

% Create transmit waveform and get the HARQ scheduling ID sequence
% from 'frcOut' structure output which also contains the waveform
% configuration and OFDM modulation parameters
[txWaveform,~,frcOut] = lteRMCULTool(frc, data);
% txWaveform = 10*txWaveform;
% figure
% plot(abs(txWaveform))

if strcmp(ue.RC,'A3-2')
    fs_real = 1.92e6;
    power_mul = 1;
elseif strcmp(ue.RC,'A3-3')
    fs_real = 3.84e6;
    power_mul = 2;
elseif strcmp(ue.RC,'A3-4')
    fs_real = 7.68e6;
    power_mul = 4;
elseif strcmp(ue.RC,'A3-5')
    fs_real = 15.36e6;
    power_mul = 8;
elseif strcmp(ue.RC, 'A3-7')
    fs_real = 30.72e6;
    power_mul = 16;
end
fs_aim = (1.92e6)*16; 

txWaveform = resample(txWaveform,fs_aim,fs_real);

%% Initiate FDD channel 
% v = 5;                    % UT velocity in km/h
DelaySpread = 35e-9;
fc1 = 1.6e9;                    % carrier frequency in Hz
c = physconst('lightspeed'); % speed of light in m/s
fd1 = (v*1000/3600)/c*fc1;     % UT max Doppler frequency in Hz

tdl1 = nrTDLChannel;
tdl1.SampleRate = fs_real;
tdl1.MaximumDopplerShift = fd1;
tdl1.DelaySpread = DelaySpread;
tdl1.DelayProfile = 'TDL-D';

tdl1.NumTransmitAntennas = 1;
tdl1.NumReceiveAntennas = 1;
tdl1.Seed = 1;

fc2 = 1.9e9;                    % carrier frequency in Hz
c = physconst('lightspeed'); % speed of light in m/s
fd2 = (v*1000/3600)/c*fc2;     % UT max Doppler frequency in Hz

tdl2 = nrTDLChannel;
tdl2.SampleRate = fs_real;
tdl2.MaximumDopplerShift = fd2;
tdl2.DelaySpread = DelaySpread;
tdl2.DelayProfile = 'TDL-D';

tdl2.NumTransmitAntennas = 1;
tdl2.NumReceiveAntennas = 1;
tdl2.Seed = 2;

fczz = 1.9e9;                    % carrier frequency in Hz
c = physconst('lightspeed'); % speed of light in m/s
fdzz = (0*1000/3600)/c*fczz;     % UT max Doppler frequency in Hz

tdlzz = nrTDLChannel;
tdlzz.SampleRate = fs_real;
tdlzz.MaximumDopplerShift = fdzz;
tdlzz.DelaySpread = DelaySpread;
tdlzz.DelayProfile = 'TDL-D';

tdlzz.NumTransmitAntennas = 1;
tdlzz.NumReceiveAntennas = 1;
tdlzz.Seed = 30;

tdlae1 = nrTDLChannel;
tdlae1.SampleRate = fs_real;
tdlae1.MaximumDopplerShift = 0;
tdlae1.DelaySpread = DelaySpread;
tdlae1.DelayProfile = 'TDL-D';

tdlae1.NumTransmitAntennas = 1;
tdlae1.NumReceiveAntennas = 1;
tdlae1.Seed = 3;

tdlae2 = nrTDLChannel;
tdlae2.SampleRate = fs_real;
tdlae2.MaximumDopplerShift = 0;
tdlae2.DelaySpread = DelaySpread;
tdlae2.DelayProfile = 'TDL-D';

tdlae2.NumTransmitAntennas = 1;
tdlae2.NumReceiveAntennas = 1;
tdlae2.Seed = 4;

c = physconst('lightspeed'); % speed of light in m/s
fd1 = (v*1000/3600)/c*fc1;     % UT max Doppler frequency in Hz

tdlbe1 = nrTDLChannel;
tdlbe1.SampleRate = fs_real;
tdlbe1.MaximumDopplerShift = fd1;
tdlbe1.DelaySpread = DelaySpread;
tdlbe1.DelayProfile = 'TDL-D';

tdlbe1.NumTransmitAntennas = 1;
tdlbe1.NumReceiveAntennas = 1;
tdlbe1.Seed = 5;

c = physconst('lightspeed'); % speed of light in m/s
fd2 = (v*1000/3600)/c*fc2;     % UT max Doppler frequency in Hz

tdlbe2 = nrTDLChannel;
tdlbe2.SampleRate = fs_real;
tdlbe2.MaximumDopplerShift = fd2;
tdlbe2.DelaySpread = DelaySpread;
tdlbe2.DelayProfile = 'TDL-D';

tdlbe2.NumTransmitAntennas = 1;
tdlbe2.NumReceiveAntennas = 1;
tdlbe2.Seed = 6;

%% Declare array size
MSE_TDD = zeros(length(snr_arr), exp_time_max);
MSE_FDD = zeros(length(snr_arr), exp_time_max);
MSE_FDD_ZZ = zeros(length(snr_arr), exp_time_max);
frame_len_rsp = 30720;

key_arr_fdd_amp_a = zeros(subcar_len,exp_time_max);
key_arr_fdd_amp_b = zeros(subcar_len,exp_time_max);
key_arr_tdd_amp_a = zeros(subcar_len,exp_time_max);
key_arr_tdd_amp_b = zeros(subcar_len,exp_time_max);
key_arr_tdd_amp_b2e = zeros(subcar_len,exp_time_max);
key_arr_fdd_amp_b2e = zeros(subcar_len,exp_time_max);

key_arr_jrn_amp_a = zeros(subcar_len,exp_time_max);
key_arr_jrn_amp_b = zeros(subcar_len,exp_time_max);
key_arr_jrn_amp_b2e = zeros(subcar_len,exp_time_max);

key_arr_fdd_a_norm_sort = zeros(subcar_len,exp_time_max);
key_arr_fdd_b_norm_sort = zeros(subcar_len,exp_time_max);
key_arr_fdd_a_norm_sort_zz = zeros(subcar_len,exp_time_max);
key_arr_fdd_b_norm_sort_zz = zeros(subcar_len,exp_time_max);
key_arr_fdd_b2e_norm_sort = zeros(subcar_len,exp_time_max);

key_arr_fdd_a_norm_sort4 = zeros(subcar_len,exp_time_max);
key_arr_fdd_b_norm_sort4 = zeros(subcar_len,exp_time_max);
key_arr_fdd_b2e_norm_sort4 = zeros(subcar_len,exp_time_max);

key_arr_tdd_a_norm_sort = zeros(subcar_len,exp_time_max);
key_arr_tdd_b_norm_sort = zeros(subcar_len,exp_time_max);
key_arr_tdd_b2e_norm_sort = zeros(subcar_len,exp_time_max);

key_arr_fdd_a = zeros(subcar_len,exp_time_max);
key_arr_fdd_b = zeros(subcar_len,exp_time_max);
key_arr_fdd_b2e = zeros(subcar_len,exp_time_max);
key_arr_tdd_a = zeros(subcar_len,exp_time_max);
key_arr_tdd_b = zeros(subcar_len,exp_time_max);
key_arr_tdd_b2e = zeros(subcar_len,exp_time_max);
key_arr_jrn_b = zeros(subcar_len,exp_time_max);
key_arr_jrn_a = zeros(subcar_len,exp_time_max);
key_arr_jrn_b2e = zeros(subcar_len,exp_time_max);
key_arr_jrn_e_divide = zeros(subcar_len,exp_time_max);

key_arr_jrn_a_norm_sort = zeros(subcar_len,exp_time_max);
key_arr_jrn_b_norm_sort = zeros(subcar_len,exp_time_max);
key_arr_jrn_e_divide_norm_sort = zeros(subcar_len,exp_time_max);

sig_b2a4_csi_arr = zeros(exp_time_max,subcar_len);
sig_a2b4_csi_arr = zeros(exp_time_max,subcar_len);
sig_b2a4_csi_arr_zz = zeros(exp_time_max,subcar_len);
sig_a2b4_csi_arr_zz = zeros(exp_time_max,subcar_len);
sig_b2a2_csi_arr = zeros(exp_time_max,subcar_len);
sig_a2b2_csi_arr = zeros(exp_time_max,subcar_len);
sig_b2a2_JRN_csi_arr = zeros(exp_time_max,subcar_len);
sig_a2b2_JRN_csi_arr = zeros(exp_time_max,subcar_len);

sig_b2e2_TDD_csi_arr = zeros(exp_time_max,subcar_len);
sig_a2e2_TDD_csi_arr = zeros(exp_time_max,subcar_len);
sig_b2e4_FDD_csi_arr = zeros(exp_time_max,subcar_len);
sig_a2e4_FDD_csi_arr = zeros(exp_time_max,subcar_len);
sig_a2e1_JRN_csi_arr = zeros(exp_time_max,subcar_len);
sig_b2e1_JRN_csi_arr = zeros(exp_time_max,subcar_len);
sig_a2e2_JRN_csi_arr = zeros(exp_time_max,subcar_len);
sig_b2e2_JRN_csi_arr = zeros(exp_time_max,subcar_len);
sig_b2e2_JRN_divid_csi_arr = zeros(exp_time_max,subcar_len);

key_tdd_amp_kdr = zeros(1,length(snr_arr));
key_fdd_amp_kdr = zeros(1,length(snr_arr));
key_fdd_amp_kgr = zeros(1,length(snr_arr));
key_tdd_amp_kgr = zeros(1,length(snr_arr));
key_fdd_amp_kdr_ae = zeros(1,length(snr_arr));
key_tdd_amp_kdr_ae = zeros(1,length(snr_arr));
key_fdd_amp_kgr_ae = zeros(1,length(snr_arr));
key_tdd_amp_kgr_ae = zeros(1,length(snr_arr));

key_fdd_kdr = zeros(1,length(snr_arr));
key_tdd_kdr = zeros(1,length(snr_arr));
key_fdd_kgr = zeros(1,length(snr_arr));
key_tdd_kgr = zeros(1,length(snr_arr));
key_jrn_kdr = zeros(1,length(snr_arr));
key_jrn_kgr = zeros(1,length(snr_arr));

key_fdd_norm_sort_kdr = zeros(1,length(snr_arr));
key_fdd_norm_sort_kgr = zeros(1,length(snr_arr));
key_fdd_norm_sort_kdr_ae = zeros(1,length(snr_arr));
key_fdd_norm_sort_kgr_ae = zeros(1,length(snr_arr));
key_fdd_norm_sort_kdr_zz = zeros(1,length(snr_arr));
key_fdd_norm_sort_kgr_zz = zeros(1,length(snr_arr));
key_fdd_norm_sort_kdr_ae_zz = zeros(1,length(snr_arr));
key_fdd_norm_sort_kgr_ae_zz = zeros(1,length(snr_arr));

key_fdd_norm_sort4_kdr = zeros(1,length(snr_arr));
key_fdd_norm_sort4_kgr = zeros(1,length(snr_arr));
key_fdd_norm_sort4_kdr_ae = zeros(1,length(snr_arr));
key_fdd_norm_sort4_kgr_ae = zeros(1,length(snr_arr));

key_tdd_norm_sort_kdr = zeros(1,length(snr_arr));
key_tdd_norm_sort_kgr = zeros(1,length(snr_arr));
key_tdd_norm_sort_kdr_ae = zeros(1,length(snr_arr));
key_tdd_norm_sort_kgr_ae = zeros(1,length(snr_arr));

key_jrn_norm_sort_kdr = zeros(1,length(snr_arr));
key_jrn_norm_sort_kgr = zeros(1,length(snr_arr));
key_jrn_norm_sort_kdr_ae_divid = zeros(1,length(snr_arr));
key_jrn_norm_sort_kgr_ae_divid = zeros(1,length(snr_arr));

key_jrn_amp_kdr = zeros(1,length(snr_arr));
key_jrn_amp_kgr = zeros(1,length(snr_arr));
key_jrn_amp_kdr_ae_divid = zeros(1,length(snr_arr));
key_jrn_amp_kgr_ae_divid = zeros(1,length(snr_arr));

key_fdd_kdr_ae = zeros(1,length(snr_arr));
key_tdd_kdr_ae = zeros(1,length(snr_arr));
key_jrn_kdr_ae = zeros(1,length(snr_arr));
key_jrn_kdr_ae_divid = zeros(1,length(snr_arr));

key_fdd_kgr_ae = zeros(1,length(snr_arr));
key_tdd_kgr_ae = zeros(1,length(snr_arr));
key_jrn_kgr_ae = zeros(1,length(snr_arr));
key_jrn_kgr_ae_divid = zeros(1,length(snr_arr));

%% Hardware Impairments
% Frequency response of hardware impairments--A1
rng(1); 
random_magnitude = rand(1, num_freq_points/2 + 1); 
random_magnitude = [random_magnitude, fliplr(random_magnitude(2:end-1))]; 
random_response = random_magnitude; 
h = real(ifft(random_response)); 
h = h(1:N);
h_a1 = h / max(abs(h)); % normalization
h_a1 = roundn(h_a1, -2);

% Frequency response of hardware impairments--A2
rng(2);
random_magnitude = rand(1, num_freq_points/2 + 1); 
random_magnitude = [random_magnitude, fliplr(random_magnitude(2:end-1))];
random_response = random_magnitude; 
h = real(ifft(random_response)); 
h = h(1:N); 
h_a2 = h / max(abs(h));
h_a2 = roundn(h_a2, -2);

% Frequency response of hardware impairments--B1
rng(3); 
random_magnitude = rand(1, num_freq_points/2 + 1);
random_magnitude = [random_magnitude, fliplr(random_magnitude(2:end-1))]; 
random_response = random_magnitude;
h = real(ifft(random_response)); 
h = h(1:N); 
h_b1 = h / max(abs(h));
h_b1 = roundn(h_b1, -2);

% Frequency response of hardware impairments--B2
rng(6); 
random_magnitude = rand(1, num_freq_points/2 + 1); 
random_magnitude = [random_magnitude, fliplr(random_magnitude(2:end-1))];
random_response = random_magnitude; 
h = real(ifft(random_response)); 
h = h(1:N); 
h_b2 = h / max(abs(h));
h_b2 = roundn(h_b2, -2);

rff_h_a1 = h_a1;
rff_h_a2 = h_a2;
rff_h_b1 = h_b1;
rff_h_b2 = h_b2;

% Generate Term Z
sig_z_csi_norm_arr = zeros(exp_time_max, 1200);
for seed_z = 1:exp_time_max
    rng((seed_z-1)*5+1)
    z_temp1 = randn(1,1200);

    rng((seed_z-1)*5+2)
    z_temp2 = randn(1,1200);

    rng((seed_z-1)*5+3)
    z_temp3 = randn(1,1200);

    rng((seed_z-1)*5+4)
    z_temp4 = randn(1,1200);

    rng((seed_z-1)*5+5)
    z_temp5 = randn(1,1200);

    mid_temp = z_temp1 + z_temp2 + z_temp3 + z_temp4 + z_temp5;

    mid_temp_min = min(mid_temp);
    mid_temp_max = max(mid_temp);
    k_val = (1.05-0.95)/(mid_temp_max - mid_temp_min);
    b_val = 0.95-k_val*mid_temp_min;

    sig_z_csi_norm_arr(seed_z,:) = (k_val*mid_temp + b_val);
end
sig_z_csi_norm_arr = sig_z_csi_norm_arr.';
sig_z_csi_norm_arr_5km = sig_z_csi_norm_arr;
save('sig_z_csi_norm_arr_5km.mat','sig_z_csi_norm_arr_5km')
toc
			
for snr_idx = 1:length(snr_arr)
    snr = snr_arr(snr_idx);
    disp(snr)
    noise_p = mean(abs(txWaveform).^2)/(10^(snr/10));

    parfor exp_time = 1:exp_time_max
        rng(exp_time)

        %% FDD loop-back 1
        sig_rff_a1 = conv(txWaveform,rff_h_a1);
        sig_rff_b1 = conv(txWaveform,rff_h_b2);

        sig_rff_pa_a1 = sig_rff_a1;
        sig_rff_pa_b1 = sig_rff_b1;

        sig_a2b1 = tdl1(sig_rff_pa_a1);

        ap_mul = sqrt(sum(abs(txWaveform(1:subframe_len)).^2) / sum(abs(sig_a2b1(1:subframe_len)).^2));
        sig_a2b1 = sig_a2b1*ap_mul;
        sig_a2b1 = sig_a2b1 + normrnd(0,sqrt(noise_p/2),[length(sig_a2b1),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_a2b1),1]);

        sig_b2a1   = tdl2(sig_rff_pa_b1);
        ap_mul = sqrt(sum(abs(txWaveform(1:subframe_len)).^2) / sum(abs(sig_b2a1(1:subframe_len)).^2));
        sig_b2a1 = sig_b2a1*ap_mul;
        sig_b2a1 = sig_b2a1 + normrnd(0,sqrt(noise_p/2),[length(sig_b2a1),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_b2a1),1]);

        % loop-back 2
        sig_rff_a2 = conv(sig_b2a1,rff_h_a2);
        sig_rff_b2 = conv(sig_a2b1,rff_h_b1);

        sig_rff_pa_a2 = sig_rff_a2;
        sig_rff_pa_b2 = sig_rff_b2;

        sig_a2b2   = tdl2(sig_rff_pa_a2);
        ap_mul = sqrt(sum(abs(sig_b2a1(1:subframe_len)).^2) / sum(abs(sig_a2b2(1:subframe_len)).^2));
        sig_a2b2 = sig_a2b2*ap_mul;
        sig_a2b2 = sig_a2b2 + normrnd(0,sqrt(noise_p/2),[length(sig_a2b2),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_a2b2),1]);

        sig_b2a2   = tdl1(sig_rff_pa_b2);
        ap_mul = sqrt(sum(abs(sig_a2b1(1:subframe_len)).^2) / sum(abs(sig_b2a2(1:subframe_len)).^2));
        sig_b2a2 = sig_b2a2*ap_mul;

        sig_b2a2 = sig_b2a2 + normrnd(0,sqrt(noise_p/2),[length(sig_b2a2),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_b2a2),1]);

        % loop-back 3
        sig_rff_a3 = conv(sig_b2a2,rff_h_a2);
        sig_rff_b3 = conv(sig_a2b2,rff_h_b1);

        sig_rff_pa_a3 = sig_rff_a3;
        sig_rff_pa_b3 = sig_rff_b3;

        sig_a2b3   = tdl2(sig_rff_pa_a3);
        ap_mul = sqrt(sum(abs(sig_b2a2(1:subframe_len)).^2) / sum(abs(sig_a2b3(1:subframe_len)).^2));
        sig_a2b3 = sig_a2b3*ap_mul;
        sig_a2b3 = sig_a2b3 + normrnd(0,sqrt(noise_p/2),[length(sig_a2b3),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_a2b3),1]);

        sig_b2a3   = tdl1(sig_rff_pa_b3);
        ap_mul = sqrt(sum(abs(sig_a2b2(1:subframe_len)).^2) / sum(abs(sig_b2a3(1:subframe_len)).^2));
        sig_b2a3 = sig_b2a3*ap_mul;
        sig_b2a3 = sig_b2a3 + normrnd(0,sqrt(noise_p/2),[length(sig_b2a3),1])  + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_b2a3),1]);

        % loop-back 4
        sig_rff_a4 = conv(sig_b2a3,rff_h_a1);
        sig_rff_b4 = conv(sig_a2b3,rff_h_b2);

        sig_rff_pa_a4 = sig_rff_a4;
        sig_rff_pa_b4 = sig_rff_b4;

        sig_a2b4 = tdl1(sig_rff_pa_a4);
        ap_mul = sqrt(sum(abs(sig_b2a3(1:subframe_len)).^2) / sum(abs(sig_a2b4(1:subframe_len)).^2));
        sig_a2b4 = sig_a2b4*ap_mul;
        sig_a2b4 = sig_a2b4 + normrnd(0,sqrt(noise_p/2),[length(sig_a2b4),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_a2b4),1]);
        sig_a2b4_csi = CSI_get_CFO_exist(sig_a2b4, subframe_len, CP_array, fft_len, fs_aim);
        sig_a2b4_csi_T = sig_a2b4_csi.';

        sig_b2a4 = tdl2(sig_rff_pa_b4);
        ap_mul = sqrt(sum(abs(sig_a2b3(1:subframe_len)).^2) / sum(abs(sig_b2a4(1:subframe_len)).^2));
        sig_b2a4 = sig_b2a4*ap_mul;
        sig_b2a4 = sig_b2a4 + normrnd(0,sqrt(noise_p/2),[length(sig_b2a4),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_b2a4),1]);
        sig_b2a4_csi = CSI_get_CFO_exist(sig_b2a4, subframe_len, CP_array, fft_len, fs_aim);
        sig_b2a4_csi_T = sig_b2a4_csi.';

        % Term Z ***********************************************************
        sig_a2b4_csi_zz = sig_a2b4_csi;
        sig_b2a4_csi_zz = sig_b2a4_csi;

        sig_a2e1 = tdlae1(sig_rff_pa_a1);
        ap_mul = sqrt(sum(abs(txWaveform(1:subframe_len)).^2) / sum(abs(sig_a2e1(1:subframe_len)).^2));
        sig_a2e1 = sig_a2e1*ap_mul;
        sig_a2e1 = sig_a2e1 + normrnd(0,sqrt(noise_p/2),[length(sig_a2e1),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_a2e1),1]);
        sig_a2e1_csi = CSI_get_CFO_exist(sig_a2e1, subframe_len, CP_array, fft_len, fs_aim);
    
        sig_b2e1 = tdlbe2(sig_rff_pa_b1);
        ap_mul = sqrt(sum(abs(txWaveform(1:subframe_len)).^2) / sum(abs(sig_b2e1(1:subframe_len)).^2));
        sig_b2e1 = sig_b2e1*ap_mul;
        sig_b2e1 = sig_b2e1 + normrnd(0,sqrt(noise_p/2),[length(sig_b2e1),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_b2e1),1]);
        sig_b2e1_csi = CSI_get_CFO_exist(sig_b2e1, subframe_len, CP_array, fft_len, fs_aim);

        sig_a2e2 = tdlae2(sig_rff_pa_a2);
        ap_mul = sqrt(sum(abs(sig_b2a1(1:subframe_len)).^2) / sum(abs(sig_a2e2(1:subframe_len)).^2));
        sig_a2e2 = sig_a2e2*ap_mul;
        sig_a2e2 = sig_a2e2 + normrnd(0,sqrt(noise_p/2),[length(sig_a2e2),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_a2e2),1]);
        sig_a2e2_csi = CSI_get_CFO_exist(sig_a2e2, subframe_len, CP_array, fft_len, fs_aim);
    
        sig_b2e2 = tdlbe1(sig_rff_pa_b2);
        ap_mul = sqrt(sum(abs(sig_a2b1(1:subframe_len)).^2) / sum(abs(sig_b2e2(1:subframe_len)).^2));
        sig_b2e2 = sig_b2e2*ap_mul;
        sig_b2e2 = sig_b2e2 + normrnd(0,sqrt(noise_p/2),[length(sig_b2e2),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_b2e2),1]);
        sig_b2e2_csi = CSI_get_CFO_exist(sig_b2e2, subframe_len, CP_array, fft_len, fs_aim);

        sig_a2e3 = tdlae2(sig_rff_pa_a3);
        ap_mul = sqrt(sum(abs(sig_b2a2(1:subframe_len)).^2) / sum(abs(sig_a2e3(1:subframe_len)).^2));
        sig_a2e3 = sig_a2e3*ap_mul;
        sig_a2e3 = sig_a2e3 + normrnd(0,sqrt(noise_p/2),[length(sig_a2e3),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_a2e3),1]);
        sig_a2e3_csi = CSI_get_CFO_exist(sig_a2e3, subframe_len, CP_array, fft_len, fs_aim);
    
        sig_b2e3 = tdlbe1(sig_rff_pa_b3);
        ap_mul = sqrt(sum(abs(sig_a2b2(1:subframe_len)).^2) / sum(abs(sig_b2e3(1:subframe_len)).^2));
        sig_b2e3 = sig_b2e3*ap_mul;
        sig_b2e3 = sig_b2e3 + normrnd(0,sqrt(noise_p/2),[length(sig_b2e3),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_b2e3),1]);
        sig_b2e3_csi = CSI_get_CFO_exist(sig_b2e3, subframe_len, CP_array, fft_len, fs_aim);

        sig_a2e4 = tdlae1(sig_rff_pa_a4);
        ap_mul = sqrt(sum(abs(sig_b2a3(1:subframe_len)).^2) / sum(abs(sig_a2e4(1:subframe_len)).^2));
        sig_a2e4 = sig_a2e4*ap_mul;
        sig_a2e4 = sig_a2e4 + normrnd(0,sqrt(noise_p/2),[length(sig_a2e4),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_a2e4),1]);
        sig_a2e4_csi = CSI_get_CFO_exist(sig_a2e4, subframe_len, CP_array, fft_len, fs_aim);
    
        sig_b2e4 = tdlbe2(sig_rff_pa_b4);
        ap_mul = sqrt(sum(abs(sig_a2b3(1:subframe_len)).^2) / sum(abs(sig_b2e4(1:subframe_len)).^2));
        sig_b2e4 = sig_b2e4*ap_mul;
        sig_b2e4 = sig_b2e4 + normrnd(0,sqrt(noise_p/2),[length(sig_b2e4),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_b2e4),1]);
        sig_b2e4_csi = CSI_get_CFO_exist(sig_b2e4, subframe_len, CP_array, fft_len, fs_aim);

        sig_b2a4_csi_arr(exp_time,:) = sig_b2a4_csi(start_subcar : end_subcar);
        sig_a2b4_csi_arr(exp_time,:) = sig_a2b4_csi(start_subcar : end_subcar);

        sig_b2a4_csi_arr_zz(exp_time,:) = sig_b2a4_csi_zz(start_subcar : end_subcar);
        sig_a2b4_csi_arr_zz(exp_time,:) = sig_a2b4_csi_zz(start_subcar : end_subcar);

        sig_b2e4_FDD_csi_arr(exp_time,:) = sig_b2e4_csi(start_subcar : end_subcar);
        sig_a2e4_FDD_csi_arr(exp_time,:) = sig_a2e4_csi(start_subcar : end_subcar);

        sig_b2e4_FDD_csi_arr(exp_time,:) = sig_b2e4_csi(start_subcar : end_subcar) .* sig_a2e4_csi(start_subcar : end_subcar) .* sig_b2e2_csi(start_subcar : end_subcar) .* sig_a2e2_csi(start_subcar : end_subcar) ./ ...
            (sig_b2e1_csi(start_subcar : end_subcar) .* sig_a2e1_csi(start_subcar : end_subcar) .* sig_b2e3_csi(start_subcar : end_subcar) .* sig_a2e3_csi(start_subcar : end_subcar));

        %% TDD loop-back 1
        sig_rff_a1 = conv(txWaveform,rff_h_a1);
        sig_rff_b1 = conv(txWaveform,rff_h_b1);

        sig_rff_pa_a1 = sig_rff_a1;
        sig_rff_pa_b1 = sig_rff_b1;

        sig_a2b1   = tdl1(sig_rff_pa_a1);
        ap_mul = sqrt(sum(abs(txWaveform(1:subframe_len)).^2) / sum(abs(sig_a2b1(1:subframe_len)).^2));
        sig_a2b1 = sig_a2b1*ap_mul;
        sig_a2b1 = sig_a2b1 + normrnd(0,sqrt(noise_p/2),[length(sig_a2b1),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_a2b1),1]);

        sig_b2a1   = tdl1(sig_rff_pa_b1);
        ap_mul = sqrt(sum(abs(txWaveform(1:subframe_len)).^2) / sum(abs(sig_b2a1(1:subframe_len)).^2));
        sig_b2a1 = sig_b2a1*ap_mul;
        sig_b2a1 = sig_b2a1 + normrnd(0,sqrt(noise_p/2),[length(sig_b2a1),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_b2a1),1]);

        % TDD loop-back 2
        sig_rff_a2 = conv(sig_b2a1,rff_h_a2);
        sig_rff_b2 = conv(sig_a2b1,rff_h_b2);

        sig_rff_pa_a2 = sig_rff_a2;
        sig_rff_pa_b2 = sig_rff_b2;

        sig_a2b2   = tdl2(sig_rff_pa_a2);
        ap_mul = sqrt(sum(abs(sig_b2a1(1:subframe_len)).^2) / sum(abs(sig_a2b2(1:subframe_len)).^2));
        sig_a2b2 = sig_a2b2*ap_mul;
        sig_a2b2 = sig_a2b2 + normrnd(0,sqrt(noise_p/2),[length(sig_a2b2),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_a2b2),1]);
        sig_a2b2_csi = CSI_get_CFO_exist(sig_a2b2, subframe_len, CP_array, fft_len, fs_aim);
        sig_a2b2_csi_T = sig_a2b2_csi.';

        sig_b2a2   = tdl2(sig_rff_pa_b2);
        ap_mul = sqrt(sum(abs(sig_a2b1(1:subframe_len)).^2) / sum(abs(sig_b2a2(1:subframe_len)).^2));
        sig_b2a2 = sig_b2a2*ap_mul;
        sig_b2a2 = sig_b2a2 + normrnd(0,sqrt(noise_p/2),[length(sig_b2a2),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_b2a2),1]);
        sig_b2a2_csi = CSI_get_CFO_exist(sig_b2a2, subframe_len, CP_array, fft_len, fs_aim);
        sig_b2a2_csi_T = sig_b2a2_csi.';

        sig_a2e2 = tdlae2(sig_rff_pa_a2);
        ap_mul = sqrt(sum(abs(sig_b2a1(1:subframe_len)).^2) / sum(abs(sig_a2e2(1:subframe_len)).^2));
        sig_a2e2 = sig_a2e2*ap_mul;
        sig_a2e2 = sig_a2e2 + normrnd(0,sqrt(noise_p/2),[length(sig_a2e2),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_a2e2),1]);
        sig_a2e2_csi = CSI_get_CFO_exist(sig_a2e2, subframe_len, CP_array, fft_len, fs_aim);

        sig_b2e2 = tdlbe2(sig_rff_pa_b2);
        ap_mul = sqrt(sum(abs(sig_a2b1(1:subframe_len)).^2) / sum(abs(sig_b2e2(1:subframe_len)).^2));
        sig_b2e2 = sig_b2e2*ap_mul;
        sig_b2e2 = sig_b2e2 + normrnd(0,sqrt(noise_p/2),[length(sig_a2e2),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_a2e2),1]);
        sig_b2e2_csi = CSI_get_CFO_exist(sig_b2e2, subframe_len, CP_array, fft_len, fs_aim);

        sig_b2a2_csi_arr(exp_time,:) = sig_b2a2_csi(start_subcar : end_subcar);
        sig_a2b2_csi_arr(exp_time,:) = sig_a2b2_csi(start_subcar : end_subcar);

        sig_b2e2_TDD_csi_arr(exp_time,:) = sig_b2e2_csi(start_subcar : end_subcar);
        sig_a2e2_TDD_csi_arr(exp_time,:) = sig_a2e2_csi(start_subcar : end_subcar);

        %% JRNSO loop-back 1
        sig_rff_a1 = conv(txWaveform,rff_h_a1);
        sig_rff_b1 = conv(txWaveform,rff_h_b2);

        sig_rff_pa_a1 = sig_rff_a1;
        sig_rff_pa_b1 = sig_rff_b1;

        sig_a2b1 = tdl1(sig_rff_pa_a1);
        ap_mul = sqrt(sum(abs(txWaveform(1:subframe_len)).^2) / sum(abs(sig_a2b1(1:subframe_len)).^2));
        sig_a2b1 = sig_a2b1*ap_mul;
        sig_a2b1 = sig_a2b1 + normrnd(0,sqrt(noise_p/2),[length(sig_a2b1),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_a2b1),1]);
        sig_a2b1_csi = CSI_get_CFO_exist(sig_a2b1, subframe_len, CP_array, fft_len, fs_aim);

        sig_b2a1 = tdl2(sig_rff_pa_b1);
        ap_mul = sqrt(sum(abs(txWaveform(1:subframe_len)).^2) / sum(abs(sig_b2a1(1:subframe_len)).^2));
        sig_b2a1 = sig_b2a1*ap_mul;
        sig_b2a1 = sig_b2a1 + normrnd(0,sqrt(noise_p/2),[length(sig_b2a1),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_a2b1),1]);
        sig_b2a1_csi = CSI_get_CFO_exist(sig_b2a1, subframe_len, CP_array, fft_len, fs_aim);

        sig_a2e1 = tdlae1(sig_rff_pa_a1);
        ap_mul = sqrt(sum(abs(txWaveform(1:subframe_len)).^2) / sum(abs(sig_a2e1(1:subframe_len)).^2));
        sig_a2e1 = sig_a2e1*ap_mul;
        sig_a2e1 = sig_a2e1 + normrnd(0,sqrt(noise_p/2),[length(sig_b2a1),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_a2b1),1]);
        sig_a2e1_csi = CSI_get_CFO_exist(sig_a2e1, subframe_len, CP_array, fft_len, fs_aim);

        sig_b2e1 = tdlbe2(sig_rff_pa_b1);
        ap_mul = sqrt(sum(abs(txWaveform(1:subframe_len)).^2) / sum(abs(sig_b2e1(1:subframe_len)).^2));
        sig_b2e1 = sig_b2e1*ap_mul;
        sig_b2e1 = sig_b2e1 + normrnd(0,sqrt(noise_p/2),[length(sig_b2a1),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_a2b1),1]);
        sig_b2e1_csi = CSI_get_CFO_exist(sig_b2e1, subframe_len, CP_array, fft_len, fs_aim);

        % JRNSO loop-back 2
        sig_rff_a2 = conv(sig_b2a1,rff_h_a1);
        sig_rff_b2 = conv(sig_a2b1,rff_h_b2);

        sig_rff_pa_a2 = sig_rff_a2;
        sig_rff_pa_b2 = sig_rff_b2;

        sig_a2b2   = tdl1(sig_rff_pa_a2);
        ap_mul = sqrt(sum(abs(sig_b2a1(1:subframe_len)).^2) / sum(abs(sig_a2b2(1:subframe_len)).^2));
        sig_a2b2 = sig_a2b2*ap_mul;
        sig_a2b2 = sig_a2b2 + normrnd(0,sqrt(noise_p/2),[length(sig_a2b2),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_a2b2),1]);
        sig_a2b2_csi = CSI_get_CFO_exist(sig_a2b2, subframe_len, CP_array, fft_len, fs_aim);

        sig_b2a2   = tdl2(sig_rff_pa_b2);
        ap_mul = sqrt(sum(abs(sig_a2b1(1:subframe_len)).^2) / sum(abs(sig_b2a2(1:subframe_len)).^2));
        sig_b2a2 = sig_b2a2*ap_mul;
        sig_b2a2 = sig_b2a2 + normrnd(0,sqrt(noise_p/2),[length(sig_b2a2),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_b2a2),1]);
        sig_b2a2_csi = CSI_get_CFO_exist(sig_b2a2, subframe_len, CP_array, fft_len, fs_aim);

        sig_a2e2   = tdlae1(sig_rff_pa_a2);
        ap_mul = sqrt(sum(abs(sig_b2a1(1:subframe_len)).^2) / sum(abs(sig_a2e2(1:subframe_len)).^2));
        sig_a2e2 = sig_a2e2*ap_mul;
        sig_a2e2 = sig_a2e2 + normrnd(0,sqrt(noise_p/2),[length(sig_b2a2),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_b2a2),1]);
        sig_a2e2_csi = CSI_get_CFO_exist(sig_a2e2, subframe_len, CP_array, fft_len, fs_aim);

        sig_b2e2   = tdlbe2(sig_rff_pa_b2);
        ap_mul = sqrt(sum(abs(sig_a2b1(1:subframe_len)).^2) / sum(abs(sig_b2e2(1:subframe_len)).^2));
        sig_b2e2 = sig_b2e2*ap_mul;
        sig_b2e2 = sig_b2e2 + normrnd(0,sqrt(noise_p/2),[length(sig_b2a2),1]) + 1j*normrnd(0,sqrt(noise_p/2),[length(sig_b2a2),1]);
        sig_b2e2_csi = CSI_get_CFO_exist(sig_b2e2, subframe_len, CP_array, fft_len, fs_aim);

        sig_b2a2_JRN_csi_arr(exp_time,:) = sig_b2a2_csi(start_subcar : end_subcar);
        sig_a2b2_JRN_csi_arr(exp_time,:) = sig_a2b2_csi(start_subcar : end_subcar);

        sig_a2e1_JRN_csi_arr(exp_time,:) = sig_a2e1_csi(start_subcar : end_subcar);
        sig_b2e1_JRN_csi_arr(exp_time,:) = sig_b2e1_csi(start_subcar : end_subcar);
        sig_a2e2_JRN_csi_arr(exp_time,:) = sig_a2e2_csi(start_subcar : end_subcar);
        sig_b2e2_JRN_csi_arr(exp_time,:) = sig_b2e2_csi(start_subcar : end_subcar);

        sig_b2e2_JRN_divid_csi_arr(exp_time,:) = sig_a2e2_JRN_csi_arr(exp_time,:).*sig_b2e2_JRN_csi_arr(exp_time,:)./(sig_a2e1_JRN_csi_arr(exp_time,:).*sig_b2e1_JRN_csi_arr(exp_time,:));

    end
    
    %% FDD Normalization
    ampli_tem = sqrt(1./mean((abs(sig_b2a4_csi_arr)).^2, 2));
    sig_b2a4_csi_arr = sig_b2a4_csi_arr .* repmat(ampli_tem,1,1200);
    ampli_tem = sqrt(1./mean((abs(sig_a2b4_csi_arr)).^2, 2));
    sig_a2b4_csi_arr = sig_a2b4_csi_arr .* repmat(ampli_tem,1,1200);
    ampli_tem = sqrt(1./mean((abs(sig_b2a4_csi_arr_zz)).^2, 2));
    sig_b2a4_csi_arr_zz = sig_b2a4_csi_arr_zz .* repmat(ampli_tem,1,1200);
    ampli_tem = sqrt(1./mean((abs(sig_a2b4_csi_arr_zz)).^2, 2));
    sig_a2b4_csi_arr_zz = sig_a2b4_csi_arr_zz .* repmat(ampli_tem,1,1200);
    ampli_tem = sqrt(1./mean((abs(sig_b2e4_FDD_csi_arr)).^2, 2));
    sig_b2e4_FDD_csi_arr = sig_b2e4_FDD_csi_arr .* repmat(ampli_tem,1,1200);
 
    %% FDD Smoothing
    sig_b2a4_csi_proc = abs(sig_b2a4_csi_arr).';
    sig_b2a4_csi_proc = smoothdata(sig_b2a4_csi_proc,1,'movmean',mv_len);
    sig_a2b4_csi_proc = abs(sig_a2b4_csi_arr).';
    sig_a2b4_csi_proc = smoothdata(sig_a2b4_csi_proc,1,'movmean',mv_len);

    sig_b2a4_csi_proc_zz = abs(sig_b2a4_csi_arr_zz).';
    sig_b2a4_csi_proc_zz = smoothdata(sig_b2a4_csi_proc_zz,1,'movmean',mv_len);
    sig_a2b4_csi_proc_zz = abs(sig_a2b4_csi_arr_zz).';
    sig_a2b4_csi_proc_zz = smoothdata(sig_a2b4_csi_proc_zz,1,'movmean',mv_len);

    sig_b2e4_FDD_csi_proc = abs(sig_b2e4_FDD_csi_arr).';
    sig_b2e4_FDD_csi_proc = smoothdata(sig_b2e4_FDD_csi_proc,1,'movmean',mv_len);

    %% Term Z
    sig_b2a4_csi_proc_zz = sig_b2a4_csi_proc_zz .* sig_z_csi_norm_arr;
    sig_a2b4_csi_proc_zz = sig_a2b4_csi_proc_zz .* sig_z_csi_norm_arr;

    %% FDD Segmentation and Sorting
    seg_len = size(sig_b2a4_csi_proc,1)/seg_num;
    thre_num = 4;

    parfor col_idx = 1:exp_time_max
        % sig_b2a4_csi_proc
        temp_arr_norm = zeros(size(sig_b2a4_csi_proc,1),1);
        for seg_idx = 1:seg_num
            a = sig_b2a4_csi_proc(1+(seg_idx-1)*seg_len:seg_idx*seg_len,col_idx);
            temp_arr_norm(1+(seg_idx-1)*seg_len:seg_idx*seg_len) = (a - min(a))/(max(a) - min(a));
        end
        temp_arr_norm_sort = sort(temp_arr_norm);

        thre_arr_b2a4_1 = mean([temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_b2a4_2 = mean([temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_b2a4_3 = mean([temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num + 1)]);

        temp_arr = zeros(length(temp_arr_norm),1);
        temp_arr(temp_arr_norm<(thre_arr_b2a4_1)) = 1;
        temp_arr((temp_arr_norm>(thre_arr_b2a4_1)) & (temp_arr_norm<(thre_arr_b2a4_2))) = 2;
        temp_arr((temp_arr_norm>(thre_arr_b2a4_2)) & (temp_arr_norm<(thre_arr_b2a4_3))) = 3;
        temp_arr(temp_arr_norm>(thre_arr_b2a4_3)) = 4;

        key_arr_fdd_a_norm_sort(:,col_idx) = temp_arr;

        % sig_a2b4_csi_proc
        temp_arr_norm = zeros(size(sig_a2b4_csi_proc,1),1);
        for seg_idx = 1:seg_num
            a = sig_a2b4_csi_proc(1+(seg_idx-1)*seg_len:seg_idx*seg_len,col_idx);
            temp_arr_norm(1+(seg_idx-1)*seg_len:seg_idx*seg_len) = (a - min(a))/(max(a) - min(a));
        end
        temp_arr_norm_sort = sort(temp_arr_norm);

        thre_arr_a2b4_1 = mean([temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_a2b4_2 = mean([temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_a2b4_3 = mean([temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num + 1)]);

        temp_arr = zeros(length(temp_arr_norm),1);
        temp_arr(temp_arr_norm<(thre_arr_b2a4_1)) = 1;
        temp_arr((temp_arr_norm>(thre_arr_b2a4_1)) & (temp_arr_norm<(thre_arr_b2a4_2))) = 2;
        temp_arr((temp_arr_norm>(thre_arr_b2a4_2)) & (temp_arr_norm<(thre_arr_b2a4_3))) = 3;
        temp_arr(temp_arr_norm>(thre_arr_b2a4_3)) = 4;

        key_arr_fdd_b_norm_sort(:,col_idx) = temp_arr;

        % sig_b2a4_csi_proc_zz
        temp_arr_norm = zeros(size(sig_b2a4_csi_proc_zz,1),1);
        for seg_idx = 1:seg_num
            a = sig_b2a4_csi_proc_zz(1+(seg_idx-1)*seg_len:seg_idx*seg_len,col_idx);
            temp_arr_norm(1+(seg_idx-1)*seg_len:seg_idx*seg_len) = (a - min(a))/(max(a) - min(a));
        end
        temp_arr_norm_sort = sort(temp_arr_norm);

        thre_arr_b2a4_zz_1 = mean([temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_b2a4_zz_2 = mean([temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_b2a4_zz_3 = mean([temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num + 1)]);

        temp_arr = zeros(length(temp_arr_norm),1);
        temp_arr(temp_arr_norm<(thre_arr_b2a4_zz_1)) = 1;
        temp_arr((temp_arr_norm>(thre_arr_b2a4_zz_1)) & (temp_arr_norm<(thre_arr_b2a4_zz_2))) = 2;
        temp_arr((temp_arr_norm>(thre_arr_b2a4_zz_2)) & (temp_arr_norm<(thre_arr_b2a4_zz_3))) = 3;
        temp_arr(temp_arr_norm>(thre_arr_b2a4_zz_3)) = 4;

        key_arr_fdd_a_norm_sort_zz(:,col_idx) = temp_arr;

        % sig_a2b4_csi_proc
        temp_arr_norm = zeros(size(sig_a2b4_csi_proc_zz,1),1);
        for seg_idx = 1:seg_num
            a = sig_a2b4_csi_proc_zz(1+(seg_idx-1)*seg_len:seg_idx*seg_len,col_idx);
            temp_arr_norm(1+(seg_idx-1)*seg_len:seg_idx*seg_len) = (a - min(a))/(max(a) - min(a));
        end
        temp_arr_norm_sort = sort(temp_arr_norm);

        thre_arr_a2b4_zz_1 = mean([temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_a2b4_zz_2 = mean([temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_a2b4_zz_3 = mean([temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num + 1)]);

        temp_arr = zeros(length(temp_arr_norm),1);
        temp_arr(temp_arr_norm<(thre_arr_a2b4_zz_1)) = 1;
        temp_arr((temp_arr_norm>(thre_arr_a2b4_zz_1)) & (temp_arr_norm<(thre_arr_a2b4_zz_2))) = 2;
        temp_arr((temp_arr_norm>(thre_arr_a2b4_zz_2)) & (temp_arr_norm<(thre_arr_a2b4_zz_3))) = 3;
        temp_arr(temp_arr_norm>(thre_arr_a2b4_zz_3)) = 4;

        key_arr_fdd_b_norm_sort_zz(:,col_idx) = temp_arr;

        % sig_b2e4_csi_proc
        temp_arr_norm = zeros(size(sig_b2e4_FDD_csi_proc,1),1);
        for seg_idx = 1:seg_num
            a = sig_b2e4_FDD_csi_proc(1+(seg_idx-1)*seg_len:seg_idx*seg_len,col_idx);
            temp_arr_norm(1+(seg_idx-1)*seg_len:seg_idx*seg_len) = (a - min(a))/(max(a) - min(a));
        end
        temp_arr_norm_sort = sort(temp_arr_norm);

        thre_arr_b2e4_1 = mean([temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_b2e4_2 = mean([temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_b2e4_3 = mean([temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num + 1)]);

        temp_arr = zeros(length(temp_arr_norm),1);
        temp_arr(temp_arr_norm<(thre_arr_b2e4_1)) = 1;
        temp_arr((temp_arr_norm>(thre_arr_b2e4_1)) & (temp_arr_norm<(thre_arr_b2e4_2))) = 2;
        temp_arr((temp_arr_norm>(thre_arr_b2e4_2)) & (temp_arr_norm<(thre_arr_b2e4_3))) = 3;
        temp_arr(temp_arr_norm>(thre_arr_b2e4_3)) = 4;

        key_arr_fdd_b2e_norm_sort(:,col_idx) = temp_arr;
    end

    key_arr_fdd_a_norm_sort_ds = key_arr_fdd_a_norm_sort(1:ds_ratio:end, :);
    key_arr_fdd_b_norm_sort_ds = key_arr_fdd_b_norm_sort(1:ds_ratio:end, :);
    key_arr_fdd_a_norm_sort_ds_zz = key_arr_fdd_a_norm_sort_zz(1:ds_ratio:end, :);
    key_arr_fdd_b_norm_sort_ds_zz = key_arr_fdd_b_norm_sort_zz(1:ds_ratio:end, :);
    key_arr_fdd_b2e_norm_sort_ds = key_arr_fdd_b2e_norm_sort(1:ds_ratio:end, :);

    key_arr_fdd_a_norm_sort_ds_bit1 = zeros(1, numel(key_arr_fdd_a_norm_sort_ds));
    key_arr_fdd_a_norm_sort_ds_bit2 = zeros(1, numel(key_arr_fdd_a_norm_sort_ds));
    parfor ele_idx = 1:numel(key_arr_fdd_a_norm_sort_ds)
        switch key_arr_fdd_a_norm_sort_ds(ele_idx)
            case 0
                disp('???')
                key_arr_fdd_a_norm_sort_ds_bit1(ele_idx) = randi([0 1],1,1);
                key_arr_fdd_a_norm_sort_ds_bit2(ele_idx) = randi([0 1],1,1);
            case 1
                key_arr_fdd_a_norm_sort_ds_bit1(ele_idx) = 0;
                key_arr_fdd_a_norm_sort_ds_bit2(ele_idx) = 0;
            case 2
                key_arr_fdd_a_norm_sort_ds_bit1(ele_idx) = 0;
                key_arr_fdd_a_norm_sort_ds_bit2(ele_idx) = 1;
            case 3
                key_arr_fdd_a_norm_sort_ds_bit1(ele_idx) = 1;
                key_arr_fdd_a_norm_sort_ds_bit2(ele_idx) = 1;
            case 4
                key_arr_fdd_a_norm_sort_ds_bit1(ele_idx) = 1;
                key_arr_fdd_a_norm_sort_ds_bit2(ele_idx) = 0;
        end
    end
    key_arr_fdd_a_norm_sort_ds_bit = [key_arr_fdd_a_norm_sort_ds_bit1; key_arr_fdd_a_norm_sort_ds_bit2];
    key_arr_fdd_a_norm_sort_ds_bit = key_arr_fdd_a_norm_sort_ds_bit(:);

    key_arr_fdd_b_norm_sort_ds_bit1 = zeros(1, numel(key_arr_fdd_b_norm_sort_ds));
    key_arr_fdd_b_norm_sort_ds_bit2 = zeros(1, numel(key_arr_fdd_b_norm_sort_ds));
    parfor ele_idx = 1:numel(key_arr_fdd_b_norm_sort_ds)
        switch key_arr_fdd_b_norm_sort_ds(ele_idx)
            case 0
                disp('???')
                key_arr_fdd_b_norm_sort_ds_bit1(ele_idx) = randi([0 1],1,1);
                key_arr_fdd_b_norm_sort_ds_bit2(ele_idx) = randi([0 1],1,1);
            case 1
                key_arr_fdd_b_norm_sort_ds_bit1(ele_idx) = 0;
                key_arr_fdd_b_norm_sort_ds_bit2(ele_idx) = 0;
            case 2
                key_arr_fdd_b_norm_sort_ds_bit1(ele_idx) = 0;
                key_arr_fdd_b_norm_sort_ds_bit2(ele_idx) = 1;
            case 3
                key_arr_fdd_b_norm_sort_ds_bit1(ele_idx) = 1;
                key_arr_fdd_b_norm_sort_ds_bit2(ele_idx) = 1;
            case 4
                key_arr_fdd_b_norm_sort_ds_bit1(ele_idx) = 1;
                key_arr_fdd_b_norm_sort_ds_bit2(ele_idx) = 0;
        end
    end
    key_arr_fdd_b_norm_sort_ds_bit = [key_arr_fdd_b_norm_sort_ds_bit1; key_arr_fdd_b_norm_sort_ds_bit2];
    key_arr_fdd_b_norm_sort_ds_bit = key_arr_fdd_b_norm_sort_ds_bit(:);

    key_arr_fdd_a_norm_sort_ds_bit_zz1 = zeros(1, numel(key_arr_fdd_a_norm_sort_ds_zz));
    key_arr_fdd_a_norm_sort_ds_bit_zz2 = zeros(1, numel(key_arr_fdd_a_norm_sort_ds_zz));
    parfor ele_idx = 1:numel(key_arr_fdd_a_norm_sort_ds_zz)
        switch key_arr_fdd_a_norm_sort_ds_zz(ele_idx)
            case 0
                disp('???')
                key_arr_fdd_a_norm_sort_ds_bit_zz1(ele_idx) = randi([0 1],1,1);
                key_arr_fdd_a_norm_sort_ds_bit_zz2(ele_idx) = randi([0 1],1,1);
            case 1
                key_arr_fdd_a_norm_sort_ds_bit_zz1(ele_idx) = 0;
                key_arr_fdd_a_norm_sort_ds_bit_zz2(ele_idx) = 0;
            case 2
                key_arr_fdd_a_norm_sort_ds_bit_zz1(ele_idx) = 0;
                key_arr_fdd_a_norm_sort_ds_bit_zz2(ele_idx) = 1;
            case 3
                key_arr_fdd_a_norm_sort_ds_bit_zz1(ele_idx) = 1;
                key_arr_fdd_a_norm_sort_ds_bit_zz2(ele_idx) = 1;
            case 4
                key_arr_fdd_a_norm_sort_ds_bit_zz1(ele_idx) = 1;
                key_arr_fdd_a_norm_sort_ds_bit_zz2(ele_idx) = 0;
        end
    end
    key_arr_fdd_a_norm_sort_ds_bit_zz = [key_arr_fdd_a_norm_sort_ds_bit_zz1; key_arr_fdd_a_norm_sort_ds_bit_zz2];
    key_arr_fdd_a_norm_sort_ds_bit_zz = key_arr_fdd_a_norm_sort_ds_bit_zz(:);

    key_arr_fdd_b_norm_sort_ds_bit_zz1 = zeros(1, numel(key_arr_fdd_b_norm_sort_ds_zz));
    key_arr_fdd_b_norm_sort_ds_bit_zz2 = zeros(1, numel(key_arr_fdd_b_norm_sort_ds_zz));
    parfor ele_idx = 1:numel(key_arr_fdd_b_norm_sort_ds_zz)
        switch key_arr_fdd_b_norm_sort_ds_zz(ele_idx)
            case 0
                disp('???')
                key_arr_fdd_b_norm_sort_ds_bit_zz1(ele_idx) = randi([0 1],1,1);
                key_arr_fdd_b_norm_sort_ds_bit_zz2(ele_idx) = randi([0 1],1,1);
            case 1
                key_arr_fdd_b_norm_sort_ds_bit_zz1(ele_idx) = 0;
                key_arr_fdd_b_norm_sort_ds_bit_zz2(ele_idx) = 0;
            case 2
                key_arr_fdd_b_norm_sort_ds_bit_zz1(ele_idx) = 0;
                key_arr_fdd_b_norm_sort_ds_bit_zz2(ele_idx) = 1;
            case 3
                key_arr_fdd_b_norm_sort_ds_bit_zz1(ele_idx) = 1;
                key_arr_fdd_b_norm_sort_ds_bit_zz2(ele_idx) = 1;
            case 4
                key_arr_fdd_b_norm_sort_ds_bit_zz1(ele_idx) = 1;
                key_arr_fdd_b_norm_sort_ds_bit_zz2(ele_idx) = 0;
        end
    end
    key_arr_fdd_b_norm_sort_ds_bit_zz = [key_arr_fdd_b_norm_sort_ds_bit_zz1; key_arr_fdd_b_norm_sort_ds_bit_zz2];
    key_arr_fdd_b_norm_sort_ds_bit_zz = key_arr_fdd_b_norm_sort_ds_bit_zz(:);

    key_arr_fdd_e_norm_sort_ds_bit1 = zeros(1, numel(key_arr_fdd_b2e_norm_sort_ds));
    key_arr_fdd_e_norm_sort_ds_bit2 = zeros(1, numel(key_arr_fdd_b2e_norm_sort_ds));

    real_idx = 1;
    parfor ele_idx = 1:numel(key_arr_fdd_b2e_norm_sort_ds)
        switch key_arr_fdd_b2e_norm_sort_ds(ele_idx)
            case 0
                disp('???')
                key_arr_fdd_e_norm_sort_ds_bit1(ele_idx) = randi([0 1],1,1);
                key_arr_fdd_e_norm_sort_ds_bit2(ele_idx) = randi([0 1],1,1);
            case 1
                key_arr_fdd_e_norm_sort_ds_bit1(ele_idx) = 0;
                key_arr_fdd_e_norm_sort_ds_bit2(ele_idx) = 0;
            case 2
                key_arr_fdd_e_norm_sort_ds_bit1(ele_idx) = 0;
                key_arr_fdd_e_norm_sort_ds_bit2(ele_idx) = 1;
            case 3
                key_arr_fdd_e_norm_sort_ds_bit1(ele_idx) = 1;
                key_arr_fdd_e_norm_sort_ds_bit2(ele_idx) = 1;
            case 4
                key_arr_fdd_e_norm_sort_ds_bit1(ele_idx) = 1;
                key_arr_fdd_e_norm_sort_ds_bit2(ele_idx) = 0;
        end
    end
    key_arr_fdd_e_norm_sort_ds_bit = [key_arr_fdd_e_norm_sort_ds_bit1; key_arr_fdd_e_norm_sort_ds_bit2];
    key_arr_fdd_e_norm_sort_ds_bit = key_arr_fdd_e_norm_sort_ds_bit(:);

    % reconciliation
    key_per_len = 1200*2/ds_ratio;
    key_num = length(key_arr_fdd_a_norm_sort_ds_bit)/key_per_len;

    key_arr_fdd_a_norm_sort_ds_bit_disc = zeros(1,key_num*infor_len);
    key_arr_fdd_b_norm_sort_ds_bit_disc = zeros(1,key_num*infor_len);
    key_arr_fdd_a_norm_sort_ds_bit_zz_disc = zeros(1,key_num*infor_len);
    key_arr_fdd_b_norm_sort_ds_bit_zz_disc = zeros(1,key_num*infor_len);
    key_arr_fdd_e_norm_sort_ds_bit_disc = zeros(1,key_num*infor_len);

    for key_idx = 1:key_num
        key_arr_fdd_a_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:(key_idx)*infor_len) = key_arr_fdd_a_norm_sort_ds_bit((key_idx-1)*key_per_len+1:(key_idx-1)*key_per_len+infor_len);
        key_arr_fdd_b_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:(key_idx)*infor_len) = key_arr_fdd_b_norm_sort_ds_bit((key_idx-1)*key_per_len+1:(key_idx-1)*key_per_len+infor_len);
        key_arr_fdd_a_norm_sort_ds_bit_zz_disc((key_idx-1)*infor_len+1:(key_idx)*infor_len) = key_arr_fdd_a_norm_sort_ds_bit_zz((key_idx-1)*key_per_len+1:(key_idx-1)*key_per_len+infor_len);
        key_arr_fdd_b_norm_sort_ds_bit_zz_disc((key_idx-1)*infor_len+1:(key_idx)*infor_len) = key_arr_fdd_b_norm_sort_ds_bit_zz((key_idx-1)*key_per_len+1:(key_idx-1)*key_per_len+infor_len);
        key_arr_fdd_e_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:(key_idx)*infor_len) = key_arr_fdd_e_norm_sort_ds_bit((key_idx-1)*key_per_len+1:(key_idx-1)*key_per_len+infor_len);
    end

    for key_idx = 1:key_num
        key_temp = key_arr_fdd_a_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:key_idx*infor_len);
        msgTx = gf(key_temp);
        key_coded = bchenc(msgTx, codeword_len, infor_len);
        key_coded(1:infor_len) = logical(key_arr_fdd_b_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:key_idx*infor_len))';
        [key_corr,~] = bchdec(key_coded, codeword_len, infor_len); 
        key_arr_fdd_b_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:key_idx*infor_len) = key_corr.x;
    end

    for key_idx = 1:key_num
        key_temp = key_arr_fdd_a_norm_sort_ds_bit_zz_disc((key_idx-1)*infor_len+1:key_idx*infor_len);
        msgTx = gf(key_temp);
        key_coded = bchenc(msgTx, codeword_len, infor_len);
        key_coded(1:infor_len) = logical(key_arr_fdd_b_norm_sort_ds_bit_zz_disc((key_idx-1)*infor_len+1:key_idx*infor_len))';
        [key_corr,~] = bchdec(key_coded, codeword_len, infor_len); 
        key_arr_fdd_b_norm_sort_ds_bit_zz_disc((key_idx-1)*infor_len+1:key_idx*infor_len) = key_corr.x;
    end

    key_arr_fdd_a_norm_sort_ds_bit_disc_final = zeros(1,key_num*key_len_final);
    key_arr_fdd_b_norm_sort_ds_bit_disc_final = zeros(1,key_num*key_len_final);
    key_arr_fdd_a_norm_sort_ds_bit_zz_disc_final = zeros(1,key_num*key_len_final);
    key_arr_fdd_b_norm_sort_ds_bit_zz_disc_final = zeros(1,key_num*key_len_final);
    key_arr_fdd_e_norm_sort_ds_bit_disc_final = zeros(1,key_num*key_len_final);

    for key_idx = 1:key_num
        key_arr_fdd_a_norm_sort_ds_bit_disc_final((key_idx-1)*key_len_final+1:(key_idx)*key_len_final) = key_arr_fdd_a_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:(key_idx-1)*infor_len+key_len_final);
        key_arr_fdd_b_norm_sort_ds_bit_disc_final((key_idx-1)*key_len_final+1:(key_idx)*key_len_final) = key_arr_fdd_b_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:(key_idx-1)*infor_len+key_len_final);
        key_arr_fdd_a_norm_sort_ds_bit_zz_disc_final((key_idx-1)*key_len_final+1:(key_idx)*key_len_final) = key_arr_fdd_a_norm_sort_ds_bit_zz_disc((key_idx-1)*infor_len+1:(key_idx-1)*infor_len+key_len_final);
        key_arr_fdd_b_norm_sort_ds_bit_zz_disc_final((key_idx-1)*key_len_final+1:(key_idx)*key_len_final) = key_arr_fdd_b_norm_sort_ds_bit_zz_disc((key_idx-1)*infor_len+1:(key_idx-1)*infor_len+key_len_final);
        key_arr_fdd_e_norm_sort_ds_bit_disc_final((key_idx-1)*key_len_final+1:(key_idx)*key_len_final) = key_arr_fdd_e_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:(key_idx-1)*infor_len+key_len_final);
    end

    key_fdd_norm_sort_kdr(snr_idx) = 1 - length(find((key_arr_fdd_a_norm_sort_ds_bit_disc_final ~= -1) & (key_arr_fdd_a_norm_sort_ds_bit_disc_final == key_arr_fdd_b_norm_sort_ds_bit_disc_final)))/length(find(key_arr_fdd_a_norm_sort_ds_bit_disc_final ~= -1 & key_arr_fdd_b_norm_sort_ds_bit_disc_final ~= -1));
    key_fdd_norm_sort_kdr_ae(snr_idx) = 1 - length(find((key_arr_fdd_a_norm_sort_ds_bit_disc_final ~= -1) & (key_arr_fdd_a_norm_sort_ds_bit_disc_final == key_arr_fdd_e_norm_sort_ds_bit_disc_final)))/length(find(key_arr_fdd_a_norm_sort_ds_bit_disc_final ~= -1 & key_arr_fdd_e_norm_sort_ds_bit_disc_final ~= -1));
    key_fdd_norm_sort_kdr_zz(snr_idx) = 1 - length(find((key_arr_fdd_a_norm_sort_ds_bit_zz_disc_final ~= -1) & (key_arr_fdd_a_norm_sort_ds_bit_zz_disc_final == key_arr_fdd_b_norm_sort_ds_bit_zz_disc_final)))/length(find(key_arr_fdd_a_norm_sort_ds_bit_zz_disc_final ~= -1 & key_arr_fdd_b_norm_sort_ds_bit_zz_disc_final ~= -1));
    key_fdd_norm_sort_kdr_ae_zz(snr_idx) = 1 - length(find((key_arr_fdd_a_norm_sort_ds_bit_zz_disc_final ~= -1) & (key_arr_fdd_a_norm_sort_ds_bit_zz_disc_final == key_arr_fdd_e_norm_sort_ds_bit_disc_final))) /length(find(key_arr_fdd_a_norm_sort_ds_bit_zz_disc_final ~= -1 & key_arr_fdd_e_norm_sort_ds_bit_disc_final ~= -1));

    %% TDD Normalization
    ampli_tem = sqrt(1./mean((abs(sig_b2a2_csi_arr)).^2, 2));
    sig_b2a2_csi_arr = sig_b2a2_csi_arr .* repmat(ampli_tem,1,1200);
    ampli_tem = sqrt(1./mean((abs(sig_a2b2_csi_arr)).^2, 2));
    sig_a2b2_csi_arr = sig_a2b2_csi_arr .* repmat(ampli_tem,1,1200);
    ampli_tem = sqrt(1./mean((abs(sig_b2e2_TDD_csi_arr)).^2, 2));
    sig_b2e2_TDD_csi_arr = sig_b2e2_TDD_csi_arr .* repmat(ampli_tem,1,1200);

    %% TDD Smoothing
    sig_b2a2_csi_proc = abs(sig_b2a2_csi_arr).';
    sig_b2a2_csi_proc = smoothdata(sig_b2a2_csi_proc,1,'movmean',mv_len);
    sig_a2b2_csi_proc = abs(sig_a2b2_csi_arr).';
    sig_a2b2_csi_proc = smoothdata(sig_a2b2_csi_proc,1,'movmean',mv_len);
    sig_b2e2_TDD_csi_proc = abs(sig_b2e2_TDD_csi_arr).';
    sig_b2e2_TDD_csi_proc = smoothdata(sig_b2e2_TDD_csi_proc,1,'movmean',mv_len);

    %% TDD Segmentation and Sorting
    parfor col_idx = 1:exp_time_max
        % sig_b2a2_csi_proc
        temp_arr_norm = zeros(size(sig_b2a2_csi_proc,1),1);
        for seg_idx = 1:seg_num
            a = sig_b2a2_csi_proc(1+(seg_idx-1)*seg_len:seg_idx*seg_len,col_idx);
            temp_arr_norm(1+(seg_idx-1)*seg_len:seg_idx*seg_len) = (a - min(a))/(max(a) - min(a));
        end
        temp_arr_norm_sort = sort(temp_arr_norm);

        thre_arr_b2a2_1 = mean([temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_b2a2_2 = mean([temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_b2a2_3 = mean([temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num + 1)]);

        temp_arr = zeros(length(temp_arr_norm),1);
        temp_arr(temp_arr_norm<(thre_arr_b2a2_1)) = 1;
        temp_arr((temp_arr_norm>(thre_arr_b2a2_1)) & (temp_arr_norm<(thre_arr_b2a2_2))) = 2;
        temp_arr((temp_arr_norm>(thre_arr_b2a2_2)) & (temp_arr_norm<(thre_arr_b2a2_3))) = 3;
        temp_arr(temp_arr_norm>(thre_arr_b2a2_3)) = 4;

        key_arr_tdd_a_norm_sort(:,col_idx) = temp_arr;

        % sig_a2b2_csi_proc
        temp_arr_norm = zeros(size(sig_a2b2_csi_proc,1),1);
        for seg_idx = 1:seg_num
            a = sig_a2b2_csi_proc(1+(seg_idx-1)*seg_len:seg_idx*seg_len,col_idx);
            temp_arr_norm(1+(seg_idx-1)*seg_len:seg_idx*seg_len) = (a - min(a))/(max(a) - min(a));
        end
        temp_arr_norm_sort = sort(temp_arr_norm);

        thre_arr_a2b2_1 = mean([temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_a2b2_2 = mean([temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_a2b2_3 = mean([temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num + 1)]);

        temp_arr = zeros(length(temp_arr_norm),1);
        temp_arr(temp_arr_norm<(thre_arr_a2b2_1)) = 1;
        temp_arr((temp_arr_norm>(thre_arr_a2b2_1)) & (temp_arr_norm<(thre_arr_a2b2_2))) = 2;
        temp_arr((temp_arr_norm>(thre_arr_a2b2_2)) & (temp_arr_norm<(thre_arr_a2b2_3))) = 3;
        temp_arr(temp_arr_norm>(thre_arr_a2b2_3)) = 4;

        key_arr_tdd_b_norm_sort(:,col_idx) = temp_arr;

        % sig_b2e2_csi_proc
        temp_arr_norm = zeros(size(sig_b2e2_TDD_csi_proc,1),1);
        for seg_idx = 1:seg_num
            a = sig_b2e2_TDD_csi_proc(1+(seg_idx-1)*seg_len:seg_idx*seg_len,col_idx);
            temp_arr_norm(1+(seg_idx-1)*seg_len:seg_idx*seg_len) = (a - min(a))/(max(a) - min(a));
        end
        temp_arr_norm_sort = sort(temp_arr_norm);

        thre_arr_b2e2_1 = mean([temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_b2e2_2 = mean([temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_b2e2_3 = mean([temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num + 1)]);

        temp_arr = zeros(length(temp_arr_norm),1);
        temp_arr(temp_arr_norm<(thre_arr_b2e2_1)) = 1;
        temp_arr(temp_arr_norm>(thre_arr_b2e2_1) & (temp_arr_norm<(thre_arr_b2e2_2))) = 2;
        temp_arr((temp_arr_norm>(thre_arr_b2e2_2)) & (temp_arr_norm<(thre_arr_b2e2_3))) = 3;
        temp_arr(temp_arr_norm>(thre_arr_b2e2_3)) = 4;

        key_arr_tdd_b2e_norm_sort(:,col_idx) = temp_arr;
    end

    key_arr_tdd_a_norm_sort_ds = key_arr_tdd_a_norm_sort(1:ds_ratio:end, :);
    key_arr_tdd_b_norm_sort_ds = key_arr_tdd_b_norm_sort(1:ds_ratio:end, :);
    key_arr_tdd_b2e_norm_sort_ds = key_arr_tdd_b2e_norm_sort(1:ds_ratio:end, :);

    key_arr_tdd_a_norm_sort_ds_bit1 = zeros(1, numel(key_arr_tdd_a_norm_sort_ds));
    key_arr_tdd_a_norm_sort_ds_bit2 = zeros(1, numel(key_arr_tdd_a_norm_sort_ds));
    parfor ele_idx = 1:numel(key_arr_tdd_a_norm_sort_ds)
        switch key_arr_tdd_a_norm_sort_ds(ele_idx)
            case 0
                disp('???')
                key_arr_tdd_a_norm_sort_ds_bit1(ele_idx) = randi([0 1],1,1);
                key_arr_tdd_a_norm_sort_ds_bit2(ele_idx) = randi([0 1],1,1);
            case 1
                key_arr_tdd_a_norm_sort_ds_bit1(ele_idx) = 0;
                key_arr_tdd_a_norm_sort_ds_bit2(ele_idx) = 0;
            case 2
                key_arr_tdd_a_norm_sort_ds_bit1(ele_idx) = 0;
                key_arr_tdd_a_norm_sort_ds_bit2(ele_idx) = 1;
            case 3
                key_arr_tdd_a_norm_sort_ds_bit1(ele_idx) = 1;
                key_arr_tdd_a_norm_sort_ds_bit2(ele_idx) = 1;
            case 4
                key_arr_tdd_a_norm_sort_ds_bit1(ele_idx) = 1;
                key_arr_tdd_a_norm_sort_ds_bit2(ele_idx) = 0;
        end
    end
    key_arr_tdd_a_norm_sort_ds_bit = [key_arr_tdd_a_norm_sort_ds_bit1; key_arr_tdd_a_norm_sort_ds_bit2];
    key_arr_tdd_a_norm_sort_ds_bit = key_arr_tdd_a_norm_sort_ds_bit(:);

    key_arr_tdd_b_norm_sort_ds_bit1 = zeros(1, numel(key_arr_tdd_b_norm_sort_ds));
    key_arr_tdd_b_norm_sort_ds_bit2 = zeros(1, numel(key_arr_tdd_b_norm_sort_ds));
    for ele_idx = 1:numel(key_arr_tdd_b_norm_sort_ds)
        switch key_arr_tdd_b_norm_sort_ds(ele_idx)
            case 0
                disp('???')
                key_arr_tdd_b_norm_sort_ds_bit1(ele_idx) = randi([0 1],1,1);
                key_arr_tdd_b_norm_sort_ds_bit2(ele_idx) = randi([0 1],1,1);
            case 1
                key_arr_tdd_b_norm_sort_ds_bit1(ele_idx) = 0;
                key_arr_tdd_b_norm_sort_ds_bit2(ele_idx) = 0;
            case 2
                key_arr_tdd_b_norm_sort_ds_bit1(ele_idx) = 0;
                key_arr_tdd_b_norm_sort_ds_bit2(ele_idx) = 1;
            case 3
                key_arr_tdd_b_norm_sort_ds_bit1(ele_idx) = 1;
                key_arr_tdd_b_norm_sort_ds_bit2(ele_idx) = 1;
            case 4
                key_arr_tdd_b_norm_sort_ds_bit1(ele_idx) = 1;
                key_arr_tdd_b_norm_sort_ds_bit2(ele_idx) = 0;
        end
    end
    key_arr_tdd_b_norm_sort_ds_bit = [key_arr_tdd_b_norm_sort_ds_bit1; key_arr_tdd_b_norm_sort_ds_bit2];
    key_arr_tdd_b_norm_sort_ds_bit = key_arr_tdd_b_norm_sort_ds_bit(:);

    key_arr_tdd_e_norm_sort_ds_bit1 = zeros(1, numel(key_arr_tdd_b2e_norm_sort_ds));
    key_arr_tdd_e_norm_sort_ds_bit2 = zeros(1, numel(key_arr_tdd_b2e_norm_sort_ds));
    for ele_idx = 1:numel(key_arr_tdd_b2e_norm_sort_ds)
        switch key_arr_tdd_b2e_norm_sort_ds(ele_idx)
            case 0
                disp('???')
                key_arr_tdd_e_norm_sort_ds_bit1(ele_idx) = randi([0 1],1,1);
                key_arr_tdd_e_norm_sort_ds_bit2(ele_idx) = randi([0 1],1,1);
            case 1
                key_arr_tdd_e_norm_sort_ds_bit1(ele_idx) = 0;
                key_arr_tdd_e_norm_sort_ds_bit2(ele_idx) = 0;
            case 2
                key_arr_tdd_e_norm_sort_ds_bit1(ele_idx) = 0;
                key_arr_tdd_e_norm_sort_ds_bit2(ele_idx) = 1;
            case 3
                key_arr_tdd_e_norm_sort_ds_bit1(ele_idx) = 1;
                key_arr_tdd_e_norm_sort_ds_bit2(ele_idx) = 1;
            case 4
                key_arr_tdd_e_norm_sort_ds_bit1(ele_idx) = 1;
                key_arr_tdd_e_norm_sort_ds_bit2(ele_idx) = 0;
        end
    end
    key_arr_tdd_e_norm_sort_ds_bit = [key_arr_tdd_e_norm_sort_ds_bit1; key_arr_tdd_e_norm_sort_ds_bit2];
    key_arr_tdd_e_norm_sort_ds_bit = key_arr_tdd_e_norm_sort_ds_bit(:);

    key_arr_tdd_a_norm_sort_ds_bit_disc = zeros(1,key_num*infor_len);
    key_arr_tdd_b_norm_sort_ds_bit_disc = zeros(1,key_num*infor_len);
    key_arr_tdd_a_norm_sort_ds_bit_zz_disc = zeros(1,key_num*infor_len);
    key_arr_tdd_b_norm_sort_ds_bit_zz_disc = zeros(1,key_num*infor_len);
    key_arr_tdd_e_norm_sort_ds_bit_disc = zeros(1,key_num*infor_len);

    for key_idx = 1:key_num
        key_arr_tdd_a_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:(key_idx)*infor_len) = key_arr_tdd_a_norm_sort_ds_bit((key_idx-1)*key_per_len+1:(key_idx-1)*key_per_len+infor_len);
        key_arr_tdd_b_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:(key_idx)*infor_len) = key_arr_tdd_b_norm_sort_ds_bit((key_idx-1)*key_per_len+1:(key_idx-1)*key_per_len+infor_len);
        key_arr_tdd_e_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:(key_idx)*infor_len) = key_arr_tdd_e_norm_sort_ds_bit((key_idx-1)*key_per_len+1:(key_idx-1)*key_per_len+infor_len);
    end

    for key_idx = 1:key_num
        key_temp = key_arr_tdd_a_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:key_idx*infor_len);
        msgTx = gf(key_temp);
        key_coded = bchenc(msgTx, codeword_len, infor_len);
        key_coded(1:infor_len) = logical(key_arr_tdd_b_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:key_idx*infor_len))';
        [key_corr,~] = bchdec(key_coded, codeword_len, infor_len); 
        key_arr_tdd_b_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:key_idx*infor_len) = key_corr.x;
    end

    key_arr_tdd_a_norm_sort_ds_bit_disc_final = zeros(1,key_num*key_len_final);
    key_arr_tdd_b_norm_sort_ds_bit_disc_final = zeros(1,key_num*key_len_final);
    key_arr_tdd_e_norm_sort_ds_bit_disc_final = zeros(1,key_num*key_len_final);

    for key_idx = 1:key_num
        key_arr_tdd_a_norm_sort_ds_bit_disc_final((key_idx-1)*key_len_final+1:(key_idx)*key_len_final) = key_arr_tdd_a_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:(key_idx-1)*infor_len+key_len_final);
        key_arr_tdd_b_norm_sort_ds_bit_disc_final((key_idx-1)*key_len_final+1:(key_idx)*key_len_final) = key_arr_tdd_b_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:(key_idx-1)*infor_len+key_len_final);
        key_arr_tdd_e_norm_sort_ds_bit_disc_final((key_idx-1)*key_len_final+1:(key_idx)*key_len_final) = key_arr_tdd_e_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:(key_idx-1)*infor_len+key_len_final);
    end

    key_tdd_norm_sort_kdr(snr_idx) = 1 - length(find((key_arr_tdd_a_norm_sort_ds_bit_disc_final ~= -1) & (key_arr_tdd_a_norm_sort_ds_bit_disc_final == key_arr_tdd_b_norm_sort_ds_bit_disc_final)))/length(find(key_arr_tdd_a_norm_sort_ds_bit_disc_final ~= -1 & key_arr_tdd_b_norm_sort_ds_bit_disc_final ~= -1));
    key_tdd_norm_sort_kdr_ae(snr_idx) = 1 - length(find((key_arr_tdd_a_norm_sort_ds_bit_disc_final ~= -1) & (key_arr_tdd_a_norm_sort_ds_bit_disc_final == key_arr_tdd_e_norm_sort_ds_bit_disc_final)))/length(find(key_arr_tdd_a_norm_sort_ds_bit_disc_final ~= -1 & key_arr_tdd_e_norm_sort_ds_bit_disc_final ~= -1));
    
    %% JRNSO Normalization
    ampli_tem = sqrt(1./mean((abs(sig_b2a2_JRN_csi_arr)).^2, 2));
    sig_b2a2_JRN_csi_arr = sig_b2a2_JRN_csi_arr .* repmat(ampli_tem,1,1200);
    ampli_tem = sqrt(1./mean((abs(sig_a2b2_JRN_csi_arr)).^2, 2));
    sig_a2b2_JRN_csi_arr = sig_a2b2_JRN_csi_arr .* repmat(ampli_tem,1,1200);
    ampli_tem = sqrt(1./mean((abs(sig_b2e2_JRN_csi_arr)).^2, 2));
    sig_b2e2_JRN_csi_arr = sig_b2e2_JRN_csi_arr .* repmat(ampli_tem,1,1200);
    ampli_tem = sqrt(1./mean((abs(sig_b2e2_JRN_divid_csi_arr)).^2, 2));
    sig_b2e2_JRN_divid_csi_arr = sig_b2e2_JRN_divid_csi_arr .* repmat(ampli_tem,1,1200);

    %% JRNSO Smoothing
    sig_b2a2_JRN_csi_proc = abs(sig_b2a2_JRN_csi_arr).';
    sig_b2a2_JRN_csi_proc = smoothdata(sig_b2a2_JRN_csi_proc,1,'movmean',mv_len);
    sig_a2b2_JRN_csi_proc = abs(sig_a2b2_JRN_csi_arr).';
    sig_a2b2_JRN_csi_proc = smoothdata(sig_a2b2_JRN_csi_proc,1,'movmean',mv_len);
    sig_b2e2_JRN_csi_proc = abs(sig_b2e2_JRN_csi_arr).';
    sig_b2e2_JRN_csi_proc = smoothdata(sig_b2e2_JRN_csi_proc,1,'movmean',mv_len);
    sig_b2e2_JRN_divid_csi_proc = abs(sig_b2e2_JRN_divid_csi_arr).';
    sig_b2e2_JRN_divid_csi_proc = smoothdata(sig_b2e2_JRN_divid_csi_proc,1,'movmean',mv_len);

    %% JRNSO Segmentation and Sorting
    parfor col_idx = 1:exp_time_max
        % sig_b2a2_csi_proc
        temp_arr_norm = zeros(size(sig_b2a2_JRN_csi_proc,1),1);
        for seg_idx = 1:seg_num
            a = sig_b2a2_JRN_csi_proc(1+(seg_idx-1)*seg_len:seg_idx*seg_len,col_idx);
            temp_arr_norm(1+(seg_idx-1)*seg_len:seg_idx*seg_len) = (a - min(a))/(max(a) - min(a));
        end
        temp_arr_norm_sort = sort(temp_arr_norm);

        thre_arr_jrn_b2a2_1 = mean([temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_jrn_b2a2_2 = mean([temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_jrn_b2a2_3 = mean([temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num + 1)]);

        temp_arr = zeros(length(temp_arr_norm),1);
        temp_arr(temp_arr_norm<(thre_arr_jrn_b2a2_1)) = 1;
        temp_arr((temp_arr_norm>(thre_arr_jrn_b2a2_1)) & (temp_arr_norm<(thre_arr_jrn_b2a2_2))) = 2;
        temp_arr((temp_arr_norm>(thre_arr_jrn_b2a2_2)) & (temp_arr_norm<(thre_arr_jrn_b2a2_3))) = 3;
        temp_arr(temp_arr_norm>(thre_arr_jrn_b2a2_3)) = 4;

        key_arr_jrn_a_norm_sort(:,col_idx) = temp_arr;

        % sig_a2b2_csi_proc
        temp_arr_norm = zeros(size(sig_a2b2_JRN_csi_proc,1),1);
        for seg_idx = 1:seg_num
            a = sig_a2b2_JRN_csi_proc(1+(seg_idx-1)*seg_len:seg_idx*seg_len,col_idx);
            temp_arr_norm(1+(seg_idx-1)*seg_len:seg_idx*seg_len) = (a - min(a))/(max(a) - min(a));
        end
        temp_arr_norm_sort = sort(temp_arr_norm);

        thre_arr_jrn_a2b2_1 = mean([temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_jrn_a2b2_2 = mean([temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_jrn_a2b2_3 = mean([temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num + 1)]);
    
        temp_arr = zeros(length(temp_arr_norm),1);
        temp_arr(temp_arr_norm<(thre_arr_jrn_a2b2_1)) = 1;
        temp_arr((temp_arr_norm>(thre_arr_jrn_a2b2_1)) & (temp_arr_norm<(thre_arr_jrn_a2b2_2))) = 2;
        temp_arr((temp_arr_norm>(thre_arr_jrn_a2b2_2)) & (temp_arr_norm<(thre_arr_jrn_a2b2_3))) = 3;
        temp_arr(temp_arr_norm>(thre_arr_jrn_a2b2_3)) = 4;

        key_arr_jrn_b_norm_sort(:,col_idx) = temp_arr;

        % sig_b2e2_csi_proc
        temp_arr_norm = zeros(size(sig_b2e2_JRN_divid_csi_proc,1),1);
        for seg_idx = 1:seg_num
            a = sig_b2e2_JRN_divid_csi_proc(1+(seg_idx-1)*seg_len:seg_idx*seg_len,col_idx);
            temp_arr_norm(1+(seg_idx-1)*seg_len:seg_idx*seg_len) = (a - min(a))/(max(a) - min(a));
        end
        temp_arr_norm_sort = sort(temp_arr_norm);

        thre_arr_jrn_b2e2_1 = mean([temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_jrn_b2e2_2 = mean([temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num + 1)]);
        thre_arr_jrn_b2e2_3 = mean([temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num + 1)]);

        temp_arr = zeros(length(temp_arr_norm),1);
        temp_arr(temp_arr_norm<(thre_arr_jrn_b2e2_1)) = 1;
        temp_arr((temp_arr_norm>(thre_arr_jrn_b2e2_1)) & (temp_arr_norm<(thre_arr_jrn_b2e2_2))) = 2;
        temp_arr((temp_arr_norm>(thre_arr_jrn_b2e2_2)) & (temp_arr_norm<(thre_arr_jrn_b2e2_3))) = 3;
        temp_arr(temp_arr_norm>(thre_arr_jrn_b2e2_3)) = 4;

        key_arr_jrn_e_divide_norm_sort(:,col_idx) = temp_arr;
    end

    key_arr_jrn_a_norm_sort_ds = key_arr_jrn_a_norm_sort(1:ds_ratio:end, :);
    key_arr_jrn_b_norm_sort_ds = key_arr_jrn_b_norm_sort(1:ds_ratio:end, :);
    key_arr_jrn_b2e_norm_sort_ds = key_arr_jrn_e_divide_norm_sort(1:ds_ratio:end, :);

    key_arr_jrn_a_norm_sort_ds_bit1 = zeros(1, numel(key_arr_jrn_a_norm_sort_ds));
    key_arr_jrn_a_norm_sort_ds_bit2 = zeros(1, numel(key_arr_jrn_a_norm_sort_ds));
    for ele_idx = 1:numel(key_arr_jrn_a_norm_sort_ds)
        switch key_arr_jrn_a_norm_sort_ds(ele_idx)
            case 0
                disp('???')
                key_arr_jrn_a_norm_sort_ds_bit1(ele_idx) = randi([0 1],1,1);
                key_arr_jrn_a_norm_sort_ds_bit2(ele_idx) = randi([0 1],1,1);
            case 1
                key_arr_jrn_a_norm_sort_ds_bit1(ele_idx) = 0;
                key_arr_jrn_a_norm_sort_ds_bit2(ele_idx) = 0;
            case 2
                key_arr_jrn_a_norm_sort_ds_bit1(ele_idx) = 0;
                key_arr_jrn_a_norm_sort_ds_bit2(ele_idx) = 1;
            case 3
                key_arr_jrn_a_norm_sort_ds_bit1(ele_idx) = 1;
                key_arr_jrn_a_norm_sort_ds_bit2(ele_idx) = 1;
            case 4
                key_arr_jrn_a_norm_sort_ds_bit1(ele_idx) = 1;
                key_arr_jrn_a_norm_sort_ds_bit2(ele_idx) = 0;
        end
    end
    key_arr_jrn_a_norm_sort_ds_bit = [key_arr_jrn_a_norm_sort_ds_bit1; key_arr_jrn_a_norm_sort_ds_bit2];
    key_arr_jrn_a_norm_sort_ds_bit = key_arr_jrn_a_norm_sort_ds_bit(:);

    key_arr_jrn_b_norm_sort_ds_bit1 = zeros(1, numel(key_arr_jrn_b_norm_sort_ds));
    key_arr_jrn_b_norm_sort_ds_bit2 = zeros(1, numel(key_arr_jrn_b_norm_sort_ds));
    for ele_idx = 1:numel(key_arr_jrn_b_norm_sort_ds)
        switch key_arr_jrn_b_norm_sort_ds(ele_idx)
            case 0
                disp('???')
                key_arr_jrn_b_norm_sort_ds_bit1(ele_idx) = randi([0 1],1,1);
                key_arr_jrn_b_norm_sort_ds_bit2(ele_idx) = randi([0 1],1,1);
            case 1
                key_arr_jrn_b_norm_sort_ds_bit1(ele_idx) = 0;
                key_arr_jrn_b_norm_sort_ds_bit2(ele_idx) = 0;
            case 2
                key_arr_jrn_b_norm_sort_ds_bit1(ele_idx) = 0;
                key_arr_jrn_b_norm_sort_ds_bit2(ele_idx) = 1;
            case 3
                key_arr_jrn_b_norm_sort_ds_bit1(ele_idx) = 1;
                key_arr_jrn_b_norm_sort_ds_bit2(ele_idx) = 1;
            case 4
                key_arr_jrn_b_norm_sort_ds_bit1(ele_idx) = 1;
                key_arr_jrn_b_norm_sort_ds_bit2(ele_idx) = 0;
        end
    end
    key_arr_jrn_b_norm_sort_ds_bit = [key_arr_jrn_b_norm_sort_ds_bit1; key_arr_jrn_b_norm_sort_ds_bit2];
    key_arr_jrn_b_norm_sort_ds_bit = key_arr_jrn_b_norm_sort_ds_bit(:);

    key_arr_jrn_e_norm_sort_ds_bit1 = zeros(1, numel(key_arr_jrn_b2e_norm_sort_ds));
    key_arr_jrn_e_norm_sort_ds_bit2 = zeros(1, numel(key_arr_jrn_b2e_norm_sort_ds));
    for ele_idx = 1:numel(key_arr_jrn_b2e_norm_sort_ds)
        switch key_arr_jrn_b2e_norm_sort_ds(ele_idx)
            case 0
                disp('???')
                key_arr_jrn_e_norm_sort_ds_bit1(ele_idx) = randi([0 1],1,1);
                key_arr_jrn_e_norm_sort_ds_bit2(ele_idx) = randi([0 1],1,1);
            case 1
                key_arr_jrn_e_norm_sort_ds_bit1(ele_idx) = 0;
                key_arr_jrn_e_norm_sort_ds_bit2(ele_idx) = 0;
            case 2
                key_arr_jrn_e_norm_sort_ds_bit1(ele_idx) = 0;
                key_arr_jrn_e_norm_sort_ds_bit2(ele_idx) = 1;
            case 3
                key_arr_jrn_e_norm_sort_ds_bit1(ele_idx) = 1;
                key_arr_jrn_e_norm_sort_ds_bit2(ele_idx) = 1;
            case 4
                key_arr_jrn_e_norm_sort_ds_bit1(ele_idx) = 1;
                key_arr_jrn_e_norm_sort_ds_bit2(ele_idx) = 0;
        end
    end
    key_arr_jrn_e_norm_sort_ds_bit = [key_arr_jrn_e_norm_sort_ds_bit1; key_arr_jrn_e_norm_sort_ds_bit2];
    key_arr_jrn_e_norm_sort_ds_bit = key_arr_jrn_e_norm_sort_ds_bit(:);

    key_arr_jrn_a_norm_sort_ds_bit_disc = zeros(1,key_num*infor_len);
    key_arr_jrn_b_norm_sort_ds_bit_disc = zeros(1,key_num*infor_len);
    key_arr_jrn_a_norm_sort_ds_bit_zz_disc = zeros(1,key_num*infor_len);
    key_arr_jrn_b_norm_sort_ds_bit_zz_disc = zeros(1,key_num*infor_len);
    key_arr_jrn_e_norm_sort_ds_bit_disc = zeros(1,key_num*infor_len);

    for key_idx = 1:key_num
        key_arr_jrn_a_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:(key_idx)*infor_len) = key_arr_jrn_a_norm_sort_ds_bit((key_idx-1)*key_per_len+1:(key_idx-1)*key_per_len+infor_len);
        key_arr_jrn_b_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:(key_idx)*infor_len) = key_arr_jrn_b_norm_sort_ds_bit((key_idx-1)*key_per_len+1:(key_idx-1)*key_per_len+infor_len);
        key_arr_jrn_e_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:(key_idx)*infor_len) = key_arr_jrn_e_norm_sort_ds_bit((key_idx-1)*key_per_len+1:(key_idx-1)*key_per_len+infor_len);
    end

    for key_idx = 1:key_num
        key_temp = key_arr_jrn_a_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:key_idx*infor_len);
        msgTx = gf(key_temp);
        key_coded = bchenc(msgTx, codeword_len, infor_len);
        key_coded(1:infor_len) = logical(key_arr_jrn_b_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:key_idx*infor_len))';
        [key_corr,~] = bchdec(key_coded, codeword_len, infor_len); 
        key_arr_jrn_b_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:key_idx*infor_len) = key_corr.x;
    end

    key_arr_jrn_a_norm_sort_ds_bit_disc_final = zeros(1,key_num*key_len_final);
    key_arr_jrn_b_norm_sort_ds_bit_disc_final = zeros(1,key_num*key_len_final);
    key_arr_jrn_a_norm_sort_ds_bit_zz_disc_final = zeros(1,key_num*key_len_final);
    key_arr_jrn_b_norm_sort_ds_bit_zz_disc_final = zeros(1,key_num*key_len_final);
    key_arr_jrn_e_norm_sort_ds_bit_disc_final = zeros(1,key_num*key_len_final);

    for key_idx = 1:key_num
        key_arr_jrn_a_norm_sort_ds_bit_disc_final((key_idx-1)*key_len_final+1:(key_idx)*key_len_final) = key_arr_jrn_a_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:(key_idx-1)*infor_len+key_len_final);
        key_arr_jrn_b_norm_sort_ds_bit_disc_final((key_idx-1)*key_len_final+1:(key_idx)*key_len_final) = key_arr_jrn_b_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:(key_idx-1)*infor_len+key_len_final);
        key_arr_jrn_e_norm_sort_ds_bit_disc_final((key_idx-1)*key_len_final+1:(key_idx)*key_len_final) = key_arr_jrn_e_norm_sort_ds_bit_disc((key_idx-1)*infor_len+1:(key_idx-1)*infor_len+key_len_final);
    end

    key_jrn_norm_sort_kdr(snr_idx) = 1 - length(find((key_arr_jrn_a_norm_sort_ds_bit_disc_final ~= -1) & (key_arr_jrn_a_norm_sort_ds_bit_disc_final == key_arr_jrn_b_norm_sort_ds_bit_disc_final)))/length(find(key_arr_jrn_a_norm_sort_ds_bit_disc_final ~= -1 & key_arr_jrn_b_norm_sort_ds_bit_disc_final ~= -1));
    key_jrn_norm_sort_kdr_ae_divid(snr_idx) = 1 - length(find((key_arr_jrn_a_norm_sort_ds_bit_disc_final ~= -1) & (key_arr_jrn_a_norm_sort_ds_bit_disc_final == key_arr_jrn_e_norm_sort_ds_bit_disc_final)))/length(find(key_arr_jrn_a_norm_sort_ds_bit_disc_final ~= -1 & key_arr_jrn_e_norm_sort_ds_bit_disc_final ~= -1));

    toc
end

MSE_DIFF = MSE_FDD - MSE_TDD;
MSE_DIFF = mean(MSE_DIFF,2);
kdr_diff = key_fdd_kdr - key_tdd_kdr;

MSE_FDD_MEAN = mean(MSE_FDD,2);
MSE_TDD_MEAN = mean(MSE_TDD,2);

kdr_term_z_recon_no_gap_5km_15000_new_Z_new_HI_new = [key_fdd_norm_sort_kdr; key_fdd_norm_sort_kdr_zz; key_tdd_norm_sort_kdr; key_jrn_norm_sort_kdr;...
                  key_fdd_norm_sort_kdr_ae; key_fdd_norm_sort_kdr_ae_zz; key_tdd_norm_sort_kdr_ae; key_jrn_norm_sort_kdr_ae_divid];

save('kdr_term_z_recon_no_gap_5km_15000_new_Z_new_HI_new.mat','kdr_term_z_recon_no_gap_5km_15000_new_Z_new_HI_new');

figure
semilogy(snr_arr, (key_fdd_norm_sort_kdr),'-r',LineWidth=1)
hold on
semilogy(snr_arr, (key_fdd_norm_sort_kdr_zz),'-g',LineWidth=1)
hold on
semilogy(snr_arr, (key_tdd_norm_sort_kdr),'-b',LineWidth=1)
hold on
semilogy(snr_arr, (key_jrn_norm_sort_kdr),'-k',LineWidth=1)
hold on
semilogy(snr_arr, (key_fdd_norm_sort_kdr_ae),'--r',LineWidth=1)
hold on
semilogy(snr_arr, (key_fdd_norm_sort_kdr_ae_zz),'--g',LineWidth=1)
hold on
semilogy(snr_arr, (key_tdd_norm_sort_kdr_ae),'--b',LineWidth=1)
hold on
semilogy(snr_arr, (key_jrn_norm_sort_kdr_ae_divid),'--k',LineWidth=1)
legend('FDD', 'FDD_ZZ','TDD','JRNSO','FDD-Eve', 'FDD-Eve_ZZ', 'TDD-Eve','JRNSO-Eve', 'Fontsize', 11)
xlabel('SNR', 'Fontsize', 11)
ylabel('KDR', 'Fontsize', 11)
title([num2str(v) ' km'])
set(gca,'FontSize',11);
ylim([10^(-4),1])
