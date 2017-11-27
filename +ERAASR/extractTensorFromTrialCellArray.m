function tensor = extractTensorFromTrialCellArray(data, idxStart, nSamples, channelIdx)
% data is either:
%  - nTrials x 1 cell array of time x channels matrices
%  - nTrials x time x channels tensor
%
% idxStart is nTrials x 1 array of indices to start extracting data from each trial
% nSamples is a scalar that indicates how many samples to take
% channelIdx is an optional mask over which channels to keep
    
    if iscell(data)
        nTrials = numel(data);
        if nargin < 4
            channelIdx = 1:size(data{1}, 2);
        end
        if isscalar(idxStart)
            idxStart = repmat(idxStart, nTrials, 1);
        end
        idxStop = idxStart + nSamples - 1;
        
        tensor = nan(nTrials, nSamples, nnz(channelIdx));
        for iTrial = 1:nTrials
            tensor(iTrial, :, :) = data{iTrial}(idxStart(iTrial):idxStop(iTrial), channelIdx);
        end
    else
        nTrials = size(data, 1);
        if nargin < 4
            channelIdx = 1:size(data, 3);
        end
        if isscalar(idxStart)
            idxStart = repmat(idxStart, nTrials, 1);
        end
        idxStop = idxStart + nSamples - 1;
        
        tensor = nan(nTrials, nSamples, nnz(channelIdx));
        for iTrial = 1:nTrials
            tensor(iTrial, :, :) = data(iTrial, idxStart(iTrial):idxStop(iTrial), channelIdx);
        end
    end

end