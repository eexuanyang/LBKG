function sig_csi = CSI_get_CFO_exist(data_rsp, subframe_len, CP_array, fft_len, fs_aim);
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

start_p_time = 1;
% [rxWaveform_rm_CFO, CFO_Val] = CFO_Offset(data_rsp(start_p_time : start_p_time + subframe_len - 1), CP_array, fft_len, fs_aim);
% [rxWaveform_rm_CFO, CFO_Val] = CFO_Offset(rxWaveform_rm_CFO, CP_array, fft_len, fs_aim);
rxWaveform_rm_CFO = (data_rsp(start_p_time : start_p_time + subframe_len - 1)).';
% figure
% plot(abs(rxWaveform_rm_CFO(1:CP_array(1))),'b');
% hold on
% plot(abs(rxWaveform_rm_CFO(1+fft_len:CP_array(1)+fft_len)),'r');

nSymbolsPerSlot = 7;
cpFraction = 0.55;
idx = 0:fft_len-1;
offset= fft_len*3+sum(CP_array(1:3));
for symbol=4:7:11

% Get cyclic prefix length in samples for the current symbol.
cpLength=CP_array(symbol);

% Position the FFT part of the way through the cyclic
% prefix; the value of cpFraction should ensure that the 
% reception is unaffected by the windowing applied in the 
% lteOFDMModulate function. Default is 55%.                
fftStart=fix(cpLength*cpFraction);

% Create vector of phase corrections, one per FFT sample,
% to compensate for FFT performed away from zero phase
% point on the original subcarriers.
phaseCorrection=repmat(exp(-1i*2*pi*(cpLength-fftStart)/fft_len*idx),size(rxWaveform_rm_CFO,1),1);                   

% Extract the appropriate section of WAVEFORM, perform the
% FFT and half-subcarrier shifting and apply the phase correction.

halfsc=repmat(exp(1i*pi/fft_len*(idx+fftStart-cpLength)),size(rxWaveform_rm_CFO,1),1);                
fftOutput = fftshift(fft(rxWaveform_rm_CFO(:,offset+fftStart+(1:fft_len)).*halfsc,[],2).*phaseCorrection,1);
% fftOutput = fftshift(fft(rxWaveform_rm_CFO(:,offset+fftStart+(1:fft_len)).*halfsc,[],2),1);

fftOutput = [fftOutput(fft_len/2+1:end), fftOutput(1:fft_len/2)];
% temp_arr  = find(abs(fftOutput)>fftOutput_thre);
% sample_arr_width_cent_temp(k,1) = temp_arr(end) - temp_arr(1);
% sample_arr_width_cent_temp(k,2) = (temp_arr(end) + temp_arr(1))/2;
% if symbol == 4
%     sample_arr(1+(k-1)*2,:) = fftOutput;
% else
%     sample_arr(2+(k-1)*2,:) = fftOutput;
% end
if symbol == 4
%                     samples1 = fftOutput(temp_arr(1):temp_arr(end));
% %     figure
% %     plot(abs(fftOutput).','k','LineWidth',2)
% %     xlabel(['Subcarrier' sprintf('\n') '()']);
% %     ylabel('Amplitude');

%     DMRS = ((abs(fftOutput)).');
%     positve_loc_arr = find(DMRS>0.5);
    sig_csi = fftOutput;

% %     [c,l] = wavedec(DMRS(positve_loc_arr(1):positve_loc_arr(end)),6,'db6');       %6层db6小波分解
% %     a6=appcoef(c,l,'db6',6);
% %     d6=detcoef(c,l,6);
% %     d5=detcoef(c,l,5);
% %     d4=detcoef(c,l,4);
% %     d3=detcoef(c,l,3);
% %     d2=detcoef(c,l,2);
% %     d1=detcoef(c,l,1);
% %     c4=[zeros(size(a6));d6;d5;d4;d3;d2;d1];
% % 
% %     wavelet_type = 'db6';
% %     a6=wrcoef('a',c,l,wavelet_type,6);         %重构第3层细节系数  
% %     d6=wrcoef('d',c,l,wavelet_type,6);         %重构第3层细节系数  
% %     d5=wrcoef('d',c,l,wavelet_type,5);         %重构第2层细节系数  
% %     d4=wrcoef('d',c,l,wavelet_type,4);         %重构第1层细节系数
% %     d3=wrcoef('d',c,l,wavelet_type,3);         %重构第3层细节系数  
% %     d2=wrcoef('d',c,l,wavelet_type,2);         %重构第2层细节系数  
% %     d1=wrcoef('d',c,l,wavelet_type,1);         %重构第1层细节系数
% % 
% %     start_p_new =   positve_loc_arr(1);
% %     end_p_new   =   positve_loc_arr(end);
% % 
% %     signal = d1+d2+d3+d4+d5+d6;
% %     figure
% %     plot(positve_loc_arr(1):positve_loc_arr(end),signal,'k','LineWidth',1);
% %     xlabel('Subcarrier');
% %     ylabel('Amplitude');

end
    offset = offset + fft_len * nSymbolsPerSlot + sum(CP_array(1:nSymbolsPerSlot));

end

end