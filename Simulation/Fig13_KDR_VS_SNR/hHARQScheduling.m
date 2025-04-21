%hHARQScheduling Update HARQ process
%   HARQPROCESS = hHARQScheduling(HARQPROCESS,SF,RVSEQUENCE) updates the
%   state of a hybrid automatic repeat request processes given the subframe
%   number SF and the redundancy version sequence RVSEQUENCE.

%   Copyright 2009-2014 The MathWorks, Inc.

function harqProcess = hHARQScheduling(harqProcess,sf,rvSequence)
    
	blkpass = ~(harqProcess.blkerr);
    ncw = harqProcess.ncw;
    
    % Allow for one RV sequence to be shared by two codewords
    if(ncw==2 && size(rvSequence, 1)==1)
        rvSequence = repmat(rvSequence, 2, 1);
    end
    
    % Transport blk size for given subframe
    trblksize = harqProcess.txConfig.TrBlkSizes(:, mod(sf, 10)+1).';
    
    % Advance to next RV and generate new data
    harqProcess.txConfig.RVIdx = mod(harqProcess.txConfig.RVIdx + 1, ...
        size(rvSequence, 2)+1);
    blkpass(harqProcess.txConfig.RVIdx == 0) = true;
    if(ncw==2)
       if(length(trblksize)==1)
           trblksize = [trblksize trblksize];
       end
       tempData = [{randi([0 1], trblksize(1), 1)} ...
           {randi([0 1], trblksize(2), 1)}];
    else
       tempData = {randi([0 1], trblksize(blkpass(1)), 1)};
    end        
    
    % Transmit new data if required
    harqProcess.data(blkpass) = tempData(blkpass);
    
    % Clear soft buffers for new transmission
    if(all(blkpass))
       harqProcess.decState = struct();
       [harqProcess.decState(1:ncw)]= deal(struct());
       [harqProcess.decState(1:ncw).CBSBuffers] = deal({});
       [harqProcess.decState(1:ncw).CBSCRC] = deal([]);
       [harqProcess.decState(1:ncw).BLKCRC] = deal([]);
    elseif(any(blkpass))
        harqProcess.decState(blkpass).CBSBuffers = {};
    end
    
    % Update the RV for new transmission
    harqProcess.txConfig.RVIdx(blkpass) = 1;
    
    % Update RVSeq for transmit side
    harqProcess.txConfig.RVSeq = rvSequence(1, ...
        harqProcess.txConfig.RVIdx(1)).';
    if(ncw==2)
        harqProcess.txConfig.RVSeq(2) = rvSequence(2, ...
            harqProcess.txConfig.RVIdx(2));
        harqProcess.txConfig.RVSeq = harqProcess.txConfig.RVSeq.';
    end
    
    % Updating RV for local decoder
    harqProcess.txConfig.RV = harqProcess.txConfig.RVSeq;
    
end
