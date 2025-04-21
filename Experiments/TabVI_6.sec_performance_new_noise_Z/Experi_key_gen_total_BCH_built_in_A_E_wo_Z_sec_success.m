close all; 
clc;
clear all;

N = 8; % 
num_freq_points = 512; % 
fs = 1; % 

codeword_len = 255; % actual number of bits per BCH codeword (code length)
infor_len = 131;
key_len_final = 128;

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

%% Set Propagation Channel Model Sampling Rate
% info = lteSCFDMAInfo(frc);
% chcfg.SamplingRate = info.SamplingRate;     

%% Processing Loop
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

v = 0;
DelaySpread = 35e-9;
fczz = 1.9e9;                    % carrier frequency in Hz
c = physconst('lightspeed'); % speed of light in m/s
fdzz = (v*1000/3600)/c*fczz;     % UT max Doppler frequency in Hz

tdlzz = nrTDLChannel;
tdlzz.SampleRate = fs_real;
tdlzz.MaximumDopplerShift = fdzz;
tdlzz.DelaySpread = DelaySpread;
% tdl2.DelaySpread = 3000e-9;
% tdl2.DelaySpread = 30000e-9;
tdlzz.DelayProfile = 'TDL-D';

tdlzz.NumTransmitAntennas = 1;
tdlzz.NumReceiveAntennas = 1;

data_set_arr = "230421";
% data_set_arr = ["230617"];
dev_arr = [1 3];
mat_arr = 4;
mv_len = 200;

seg_num = 100;
ds_ratio = 16;
thre_num = 4;
sample_len = 1;
final_key_len = 128;

key_arr_a = [];
key_arr_e = [];

seed_z_arr = 1:10;


for data_set = data_set_arr
    if data_set == "230421"
        loc_arr = 1:17;
    elseif data_set == "230428"
        loc_arr = 1:94;
    elseif data_set == "230602"
        loc_arr = 1:15;
    elseif data_set == "230617"
        loc_arr = 1:16;
    end
    for dev = dev_arr
        for mat_index = mat_arr
            key_arr = zeros(1200*length(loc_arr),sample_len);
            sample_idx = 1;
            for loc = loc_arr
                if dev == 1
                    load(['File Path:' num2str(str2num(data_set)) '\proceed\user' num2str(dev) '\data' num2str(loc) '\proceed_'  num2str(mat_index)  '.mat']);
                    DMRS_arr = DMRS_arr_par(1:sample_len,:);
                    
                    % Power normalization
                    ampli_tem = sqrt(1./mean((abs(DMRS_arr)).^2));
                    DMRS_arr = DMRS_arr * ampli_tem;
                    DMRS_arr_smooth = smoothdata(abs(DMRS_arr.'),1,'movmean',mv_len);
                    
%                         DMRS_arr_proc = DMRS_arr_smooth.*sig_z_csi_norm(:, seed_z);
                    DMRS_arr_proc = DMRS_arr_smooth;
                else
                    if loc == 4 || loc == 7 || loc == 12 || loc == 17
                        dev_exact = 3;
                    else
                        dev_exact = 31;
                    end
        
                    load(['I:\Project\Key\' num2str(str2num(data_set)) '\proceed\user' num2str(dev_exact) '\data' num2str(loc) '\proceed_1.mat']);
                    DMRS_arr1 = DMRS_arr_par(1:sample_len,:);
                    DMRS_arr_proc1 = smoothdata(abs(DMRS_arr1.'),1,'movmean',mv_len);
                    load(['I:\Project\Key\' num2str(str2num(data_set)) '\proceed\user' num2str(dev_exact) '\data' num2str(loc) '\proceed_2.mat']);
                    DMRS_arr2 = DMRS_arr_par(1:sample_len,:);
                    DMRS_arr_proc2 = smoothdata(abs(DMRS_arr2.'),1,'movmean',mv_len);
                    load(['I:\Project\Key\' num2str(str2num(data_set)) '\proceed\user' num2str(dev_exact) '\data' num2str(loc) '\proceed_3.mat']);
                    DMRS_arr3 = DMRS_arr_par(1:sample_len,:);
                    DMRS_arr_proc3 = smoothdata(abs(DMRS_arr3.'),1,'movmean',mv_len);
                    load(['I:\Project\Key\' num2str(str2num(data_set)) '\proceed\user' num2str(dev_exact) '\data' num2str(loc) '\proceed_4.mat']);
                    DMRS_arr4 = DMRS_arr_par(1:sample_len,:);
                    DMRS_arr_proc4 = smoothdata(abs(DMRS_arr4.'),1,'movmean',mv_len);
                    load(['I:\Project\Key\' num2str(str2num(data_set)) '\proceed\user' num2str(dev_exact) '\data' num2str(loc) '\proceed_5.mat']);
                    DMRS_arr5 = DMRS_arr_par(1:sample_len,:);
                    DMRS_arr_proc5 = smoothdata(abs(DMRS_arr5.'),1,'movmean',mv_len);
                    load(['I:\Project\Key\' num2str(str2num(data_set)) '\proceed\user' num2str(dev_exact) '\data' num2str(loc) '\proceed_6.mat']);
                    DMRS_arr6 = DMRS_arr_par(1:sample_len,:);
                    DMRS_arr_proc6 = smoothdata(abs(DMRS_arr6.'),1,'movmean',mv_len);
                    load(['I:\Project\Key\' num2str(str2num(data_set)) '\proceed\user' num2str(dev_exact) '\data' num2str(loc) '\proceed_7.mat']);
                    DMRS_arr7 = DMRS_arr_par(1:sample_len,:);
                    DMRS_arr_proc7 = smoothdata(abs(DMRS_arr7.'),1,'movmean',mv_len);
                    load(['I:\Project\Key\' num2str(str2num(data_set)) '\proceed\user' num2str(dev_exact) '\data' num2str(loc) '\proceed_8.mat']);
                    DMRS_arr8 = DMRS_arr_par(1:sample_len,:);
                    DMRS_arr_proc8 = smoothdata(abs(DMRS_arr8.'),1,'movmean',mv_len);
                    
                    DMRS_arr = DMRS_arr_proc7 .* DMRS_arr_proc8 .* DMRS_arr_proc3 .* DMRS_arr_proc4 ./ (DMRS_arr_proc1 .* DMRS_arr_proc2 .* DMRS_arr_proc5 .* DMRS_arr_proc6);
                
                    % Power normalization
                    ampli_tem = sqrt(1./mean((abs(DMRS_arr)).^2));
                    DMRS_arr = DMRS_arr * ampli_tem;
                    DMRS_arr_smooth = smoothdata(abs(DMRS_arr.'),1,'movmean',mv_len);
                    
                    DMRS_arr_proc = DMRS_arr_smooth.';
                end

                for col_idx = 1:size(DMRS_arr_proc,2)
                    temp_arr_norm = zeros(size(DMRS_arr_proc,1),1);
                    
                    seg_len = round(length(temp_arr_norm)/seg_num);
                    for seg_idx = 1:seg_num
                        a = DMRS_arr_proc(1+(seg_idx-1)*seg_len:seg_idx*seg_len,col_idx);
                        temp_arr_norm(1+(seg_idx-1)*seg_len:seg_idx*seg_len) = (a - min(a))/(max(a) - min(a));
                    end
                    temp_arr_norm_sort = sort(temp_arr_norm);
                    
                    thre_arr(1) = mean([temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(length(temp_arr_norm_sort)/thre_num + 1)]);
                    thre_arr(2) = mean([temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(2*length(temp_arr_norm_sort)/thre_num + 1)]);
                    thre_arr(3) = mean([temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num), temp_arr_norm_sort(3*length(temp_arr_norm_sort)/thre_num + 1)]);
                    
                    temp_arr = zeros(length(temp_arr_norm),1);
                    temp_arr(temp_arr_norm<(thre_arr(1))) = 1;
                    temp_arr((temp_arr_norm>(thre_arr(1))) & (temp_arr_norm<(thre_arr(2)))) = 2;
                    temp_arr((temp_arr_norm>(thre_arr(2))) & (temp_arr_norm<(thre_arr(3)))) = 3;
                    temp_arr(temp_arr_norm>(thre_arr(3))) = 4;
            
                    key_arr(1+(sample_idx - 1)*length(temp_arr_norm):sample_idx*length(temp_arr_norm),col_idx) = temp_arr;
                end
                sample_idx = sample_idx + 1;
            end
        end
        if dev == 1
            key_arr_a = [key_arr_a; key_arr];
        else
            key_arr_e = [key_arr_e; key_arr];
        end
    end
end

key_arr_a_ds = key_arr_a(1:ds_ratio:end, :);
key_arr_e_ds = key_arr_e(1:ds_ratio:end, :);

key_arr_a_bit = zeros(1, 2*numel(key_arr_a_ds));
real_idx = 1;
for ele_idx = 1:numel(key_arr_a_ds)
    switch key_arr_a_ds(ele_idx)
        case 0
            disp('???')
            key_arr_a_bit(1+(real_idx-1)*2) = randi([0 1],1,1);
            key_arr_a_bit(2+(real_idx-1)*2) = randi([0 1],1,1);
            real_idx = real_idx + 1;
        case 1
            key_arr_a_bit(1+(real_idx-1)*2) = 0;
            key_arr_a_bit(2+(real_idx-1)*2) = 0;
            real_idx = real_idx + 1;
        case 2
            key_arr_a_bit(1+(real_idx-1)*2) = 0;
            key_arr_a_bit(2+(real_idx-1)*2) = 1;
            real_idx = real_idx + 1;
        case 3
            key_arr_a_bit(1+(real_idx-1)*2) = 1;
            key_arr_a_bit(2+(real_idx-1)*2) = 1;
            real_idx = real_idx + 1;
        case 4
            key_arr_a_bit(1+(real_idx-1)*2) = 1;
            key_arr_a_bit(2+(real_idx-1)*2) = 0;
            real_idx = real_idx + 1;
    end
end

key_arr_e_bit = zeros(1, 2*numel(key_arr_e_ds));
real_idx = 1;
for ele_idx = 1:numel(key_arr_e_ds)
    switch key_arr_e_ds(ele_idx)
        case 0
            disp('???')
            key_arr_e_bit(1+(real_idx-1)*2) = randi([0 1],1,1);
            key_arr_e_bit(2+(real_idx-1)*2) = randi([0 1],1,1);
            real_idx = real_idx + 1;
        case 1
            key_arr_e_bit(1+(real_idx-1)*2) = 0;
            key_arr_e_bit(2+(real_idx-1)*2) = 0;
            real_idx = real_idx + 1;
        case 2
            key_arr_e_bit(1+(real_idx-1)*2) = 0;
            key_arr_e_bit(2+(real_idx-1)*2) = 1;
            real_idx = real_idx + 1;
        case 3
            key_arr_e_bit(1+(real_idx-1)*2) = 1;
            key_arr_e_bit(2+(real_idx-1)*2) = 1;
            real_idx = real_idx + 1;
        case 4
            key_arr_e_bit(1+(real_idx-1)*2) = 1;
            key_arr_e_bit(2+(real_idx-1)*2) = 0;
            real_idx = real_idx + 1;
    end
end

key_arr_a_bit_vec = (key_arr_a_bit(:)).';
key_arr_e_bit_vec = (key_arr_e_bit(:)).';

key_per_len = 1200*2/ds_ratio;
key_num = length(key_arr_a_bit_vec)/key_per_len;

key_arr_a_bit_vec_discard = zeros(1,infor_len*key_num);
key_arr_e_bit_vec_discard = zeros(1,infor_len*key_num);

for key_idx = 1:key_num
    key_arr_a_bit_vec_discard((key_idx-1)*infor_len+1:infor_len*key_idx) = key_arr_a_bit_vec((key_idx-1)*key_per_len+1:(key_idx-1)*key_per_len+infor_len);
    key_arr_e_bit_vec_discard((key_idx-1)*infor_len+1:infor_len*key_idx) = key_arr_e_bit_vec((key_idx-1)*key_per_len+1:(key_idx-1)*key_per_len+infor_len);
end

key_arr_fdd_a_norm_sort_ds_bit_disc_final = zeros(1,key_num*key_len_final);
key_arr_fdd_e_norm_sort_ds_bit_disc_final = zeros(1,key_num*key_len_final);

for key_idx = 1:key_num
    key_arr_fdd_a_norm_sort_ds_bit_disc_final((key_idx-1)*key_len_final+1:(key_idx)*key_len_final) = key_arr_a_bit_vec_discard((key_idx-1)*infor_len+1:(key_idx-1)*infor_len+key_len_final);
    key_arr_fdd_e_norm_sort_ds_bit_disc_final((key_idx-1)*key_len_final+1:(key_idx)*key_len_final) = key_arr_e_bit_vec_discard((key_idx-1)*infor_len+1:(key_idx-1)*infor_len+key_len_final);
end

key_fdd_norm_sort_A_E_wo_Z_kdr = 1 - length(find((key_arr_fdd_a_norm_sort_ds_bit_disc_final == key_arr_fdd_e_norm_sort_ds_bit_disc_final)))/length(key_arr_fdd_a_norm_sort_ds_bit_disc_final);
key_fdd_norm_sort_A_E_wo_Z_kgr = 2*length(find((key_arr_fdd_a_norm_sort_ds_bit_disc_final == key_arr_fdd_e_norm_sort_ds_bit_disc_final)))/length(key_arr_fdd_a_norm_sort_ds_bit_disc_final);
 
save('key_fdd_norm_sort_A_E_wo_Z_kdr.mat','key_fdd_norm_sort_A_E_wo_Z_kdr')
save('key_fdd_norm_sort_A_E_wo_Z_kgr.mat','key_fdd_norm_sort_A_E_wo_Z_kgr')

figure
plot(key_arr_a_bit_vec_discard(1:1000))