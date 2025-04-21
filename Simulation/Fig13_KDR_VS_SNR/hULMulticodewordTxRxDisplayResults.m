%hULMulticodewordTxRxDisplayResults Display results
%   hULMulticodewordTxRxDisplayResults(CRC,TXCQI,TXRI,TXACK,RXCQI,RXRI,RXACK)
%   This functions displays the crc results in vector CRC and displays the
%   contents of the transmitted and received CQI, RI and ACK vectors
%   TXCQI, TXRI, TXACK, RXCQI, RXRI and RXACK

%   Copyright 2011-2014 The MathWorks, Inc.

function hULMulticodewordTxRxDisplayResults(crc,txCQI,txRI,txACK,rxCQI,rxRI,rxACK)

disp('CRCs:');
disp(['Codeword 1: ' num2str(crc(1))]);
fprintf(['Codeword 2: ' num2str(crc(2)) '\n\n'])

disp('CQI:')
disp(['transmitted: ' num2str(txCQI.')])
fprintf(['received   : ' num2str(rxCQI.') '\n\n'])

disp('RI:')
disp(['transmitted: ' num2str(txRI.')])
fprintf(['received   : ' num2str(rxRI.') '\n\n'])

disp('ACK:')
disp(['transmitted: ' num2str(txACK.')])
fprintf(['received   : ' num2str(rxACK.') '\n\n'])

end
