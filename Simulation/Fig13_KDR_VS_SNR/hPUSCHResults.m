%hPUSCHResults Display PUSCH conformance test results
%   hPUSCHResults(SNRIN,NFRAMES,TRBLKSIZES,THROUGHPUT,BITTHROUGHPUT) plots
%   the conformance test results given a range of SNR points SNRIN, a
%   number of frames simulated for each point NFRAMES, the transport block
%   sizes used for each subframe TRBLKSIZES and the results THROUGHPUT and
%   BITTHROUGHPUT.

%   Copyright 2010-2014 The MathWorks, Inc.

function hPUSCHResults(SNRIn, NFrames, trBlkSizes, throughput, bitThroughput)

    figure;
    plot(SNRIn, mean(bitThroughput, 2)/1000,'-*');
    title(['Throughput for ', num2str(NFrames) ' Frame(s)']);
    xlabel('SNR (dB)'); ylabel('Throughput (Mbps)');
    grid on;
    hold on;
    plot(SNRIn, mean(trBlkSizes)*0.7*ones(1, numel(SNRIn))/1000, '--rs');
    set(gca, 'YLim', [0 mean(trBlkSizes)*1.05/1000]);
    legend('Simulation Result', '70 Percent Throughput', ...
        'Location', 'SouthEast')
    
    figure;
    plot(SNRIn, throughput,'-*');
    title(['Throughput for ', num2str(NFrames) ' Frame(s)']);
    xlabel('SNR (dB)'); ylabel('Throughput (%)');
    grid on;
    hold on;
    plot(SNRIn, 70*ones(1, numel(SNRIn)), '--rs');
    set(gca, 'YLim', [0 105]);
    legend('Simulation Result', '70 Percent Throughput', ...
        'Location', 'SouthEast');

end