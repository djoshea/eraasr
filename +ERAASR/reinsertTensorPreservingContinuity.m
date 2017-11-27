function into = reinsertTensorPreservingContinuity(into, insert, idxStart)
% into is either:
%   - nTrials x 1 cell array of time x channels matrices
%   - nTrials x time x channels tensor
% insert is nTrials x time x channels tensor
% idxStart is nTrials x 1 array of indices to start extracting data from
% each trial

    nSamplesInsert = size(insert, 2);

    if iscell(into)
        nTrials = numel(into);
        
        for iR = 1:nTrials
            % preserve existing delta to the first sample of the inserted
            % region and to the first sample after the inserted region
            correctionLeft = shiftdim(into{iR}(idxStart(iR), :), -1) - insert(iR, 1, :);
            lastInd = idxStart(iR) + nSamplesInsert - 1;
            if lastInd < size(into{iR}, 1)
                correctionRight = insert(iR, end, :) + correctionLeft - shiftdim(into{iR}(lastInd, :), -1);
                into{iR}(lastInd+1:end, :) = bsxfun(@plus, into{iR}(lastInd+1:end, :), shiftdim(correctionRight, 1));
            end
            into{iR}( idxStart(iR) + (1:nSamplesInsert), :) = bsxfun(@plus, insert(iR, :, :), correctionLeft);

        end
        
    else
        nTrials = size(into, 1);
        for iR = 1:nTrials
            correctionLeft = into(iR, idxStart(iR), :) - insert(iR, 1, :);
            lastInd = idxStart(iR) + nSamplesInsert - 1;
            if lastInd < size(into, 2)
                correctionRight = insert(iR, end, :) + correctionLeft - into(iR, lastInd, :);
                into(iR, lastInd+1:end, :) = bsxfun(@plus, into(iR, lastInd+1:end, :), correctionRight);
            end
            into(iR, idxStart(iR) + (1:nSamplesInsert), :) = bsxfun(@plus, insert(iR, :, :), correctionLeft);
        end
    end
    
end

% was 1 to 3 at the right edge, now insert ends with
% 5+correctionLeft...need to add 5+correctionLeft - 1 to right edge