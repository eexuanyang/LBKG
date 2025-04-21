% close all; 
clc;
clear all;
tic

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

% data_set_arr = ["230421" "230428" "230602" "230617"];
data_set_arr = "230421";
dev_arr = 1:2;
mat_arr = 4;
mv_len = 200;
% seg_num = 50;
ds_ratio = 16;
thre_num = 4;
sample_len = 1;
final_key_len = 128;

seed_z_arr = 1:10;

load('sig_z_csi_norm.mat')

seg_arr = 100;
% seg_arr = 5;
for seg_idx = 1:length(seg_arr)
    seg_num = seg_arr(seg_idx);
    key_fdd_norm_sort_kdr_zz = zeros(1,length(seed_z_arr));
    key_fdd_norm_sort_kgr_zz = zeros(1,length(seed_z_arr));
    for seed_z = seed_z_arr
        key_arr_a = [];
        key_arr_b = [];
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
                        load(['File Path:' num2str(str2num(data_set)) '\proceed\user' num2str(dev) '\data' num2str(loc) '\proceed_'  num2str(mat_index)  '.mat']);
        
                        %% Power normalization
                        DMRS_arr = DMRS_arr_par(1:sample_len,:);
                        ampli_tem = sqrt(1./mean((abs(DMRS_arr)).^2));
                        DMRS_arr = DMRS_arr * ampli_tem;
                        DMRS_arr_smooth = smoothdata(abs(DMRS_arr.'),1,'movmean',mv_len);
                        
                        DMRS_arr_proc = DMRS_arr_smooth.*sig_z_csi_norm(:, seed_z);

                        for col_idx = 1:size(DMRS_arr_proc,2)
                            temp_arr_norm = zeros(size(DMRS_arr_proc,1),1);
                            
                            seg_len = round(length(temp_arr_norm)/seg_num);
                            for seg_idx_temp = 1:seg_num
                                a = DMRS_arr_proc(1+(seg_idx_temp-1)*seg_len:seg_idx_temp*seg_len,col_idx);
                                temp_arr_norm(1+(seg_idx_temp-1)*seg_len:seg_idx_temp*seg_len) = (a - min(a))/(max(a) - min(a));
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
                    key_arr_b = [key_arr_b; key_arr];
                end
            end
        end
        
        key_arr_a_ds = key_arr_a(1:ds_ratio:end, :);
        key_arr_b_ds = key_arr_b(1:ds_ratio:end, :);
        
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
        
        key_arr_b_bit = zeros(1, 2*numel(key_arr_b_ds));
        real_idx = 1;
        for ele_idx = 1:numel(key_arr_b_ds)
            switch key_arr_b_ds(ele_idx)
                case 0
                    disp('???')
                    key_arr_b_bit(1+(real_idx-1)*2) = randi([0 1],1,1);
                    key_arr_b_bit(2+(real_idx-1)*2) = randi([0 1],1,1);
                    real_idx = real_idx + 1;
                case 1
                    key_arr_b_bit(1+(real_idx-1)*2) = 0;
                    key_arr_b_bit(2+(real_idx-1)*2) = 0;
                    real_idx = real_idx + 1;
                case 2
                    key_arr_b_bit(1+(real_idx-1)*2) = 0;
                    key_arr_b_bit(2+(real_idx-1)*2) = 1;
                    real_idx = real_idx + 1;
                case 3
                    key_arr_b_bit(1+(real_idx-1)*2) = 1;
                    key_arr_b_bit(2+(real_idx-1)*2) = 1;
                    real_idx = real_idx + 1;
                case 4
                    key_arr_b_bit(1+(real_idx-1)*2) = 1;
                    key_arr_b_bit(2+(real_idx-1)*2) = 0;
                    real_idx = real_idx + 1;
            end
        end
        
        key_arr_a_bit_vec = (key_arr_a_bit(:)).';
        key_arr_b_bit_vec = (key_arr_b_bit(:)).';
        
        key_per_len = 1200*2/ds_ratio;
        key_num = length(key_arr_a_bit_vec)/key_per_len;
        
        key_arr_a_bit_vec_discard = zeros(1,infor_len*key_num);
        key_arr_b_bit_vec_discard = zeros(1,infor_len*key_num);
        
        for key_idx = 1:key_num
            key_arr_a_bit_vec_discard((key_idx-1)*infor_len+1:infor_len*key_idx) = key_arr_a_bit_vec((key_idx-1)*key_per_len+1:(key_idx-1)*key_per_len+infor_len);
            key_arr_b_bit_vec_discard((key_idx-1)*infor_len+1:infor_len*key_idx) = key_arr_b_bit_vec((key_idx-1)*key_per_len+1:(key_idx-1)*key_per_len+infor_len);
        end
            
        for key_idx = 1:key_num
            key_temp = key_arr_a_bit_vec_discard((key_idx-1)*infor_len+1:key_idx*infor_len);
            msgTx = gf(key_temp);
            key_coded = bchenc(msgTx, codeword_len, infor_len);
            key_coded(1:infor_len) = logical(key_arr_b_bit_vec_discard((key_idx-1)*infor_len+1:key_idx*infor_len))';
            [key_corr,~] = bchdec(key_coded, codeword_len, infor_len); 
            key_arr_b_bit_vec_discard((key_idx-1)*infor_len+1:key_idx*infor_len) = key_corr.x;
        end
    
        key_arr_fdd_a_norm_sort_ds_bit_disc_final = zeros(1,key_num*key_len_final);
        key_arr_fdd_b_norm_sort_ds_bit_disc_final = zeros(1,key_num*key_len_final);
    
        for key_idx = 1:key_num
            key_arr_fdd_a_norm_sort_ds_bit_disc_final((key_idx-1)*key_len_final+1:(key_idx)*key_len_final) = key_arr_a_bit_vec_discard((key_idx-1)*infor_len+1:(key_idx-1)*infor_len+key_len_final);
            key_arr_fdd_b_norm_sort_ds_bit_disc_final((key_idx-1)*key_len_final+1:(key_idx)*key_len_final) = key_arr_b_bit_vec_discard((key_idx-1)*infor_len+1:(key_idx-1)*infor_len+key_len_final);
        end
        
        key_fdd_norm_sort_kdr_zz(seed_z) = 1 - length(find((key_arr_fdd_a_norm_sort_ds_bit_disc_final == key_arr_fdd_b_norm_sort_ds_bit_disc_final)))/length(key_arr_fdd_a_norm_sort_ds_bit_disc_final);
        key_fdd_norm_sort_kgr_zz(seed_z) = 2*length(find((key_arr_fdd_a_norm_sort_ds_bit_disc_final == key_arr_fdd_b_norm_sort_ds_bit_disc_final)))/length(key_arr_fdd_a_norm_sort_ds_bit_disc_final);
    
        dlmwrite(['key_arr_fdd_seg' num2str(seg_num) '_seed' num2str(seed_z) '.out'], key_arr_fdd_a_norm_sort_ds_bit_disc_final, ' ');
    
    end

    save(['key_fdd_norm_sort_kdr_zz_seg' num2str(seg_num) '.mat'],'key_fdd_norm_sort_kdr_zz')
    save(['key_fdd_norm_sort_kgr_zz_seg' num2str(seg_num) '.mat'],'key_fdd_norm_sort_kgr_zz')

    toc
end

figure
plot(key_arr_a_bit_vec_discard(1:1000))

toc
