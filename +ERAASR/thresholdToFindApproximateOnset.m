function startIndex = thresholdToFindApproximateOnset(data, Fs, thresholdValue, varargin)
% Threshold a single channel to get a slightly better estimate of start time
% data is either:
%   - nTrials x 1 cell of nTime x nChannels matrices
%   - nTrials x nTime x nChannels tensor
%
% startIndex is nTrials x 1 vector of threshold crossing indices, NaN if
% threshold not crossed

    p = inputParser();
    p.addRequired('data', @(x) iscell(x) || isnumeric(x));
    p.addRequired('Fs', @isscalar);
    p.addRequired('thresholdValue', @isscalar);
    p.addParameter('thresholdChannel', 1, @isscalar);
    p.addParameter('hpCornerHz', 250, @isscalar);
    p.addParameter('hpOrder', 4, @isscalar);
    p.parse(data, Fs, thresholdValue, varargin{:});

    if iscell(data)
        dataThreshold = cellfun(@(d) d(:, p.Results.thresholdChannel), data, 'UniformOutput', false);
    else
        dataThreshold = data(:, :, p.Results.thresholdChannel);
    end
 
    if p.Results.hpCornerHz > 0
        dataThreshold = ERAASR.highPassFilter(dataThreshold, Fs, 'cornerHz', p.Results.hpCornerHz, ...
            'order', p.Results.hpOrder, 'filtfilt', true, 'subtractFirstSample', true);
    end

    thresh = p.Results.thresholdValue;

    % Find threshold crossings
    if iscell(dataThreshold)
        startIndex = nan(numel(data), 1);

        for iT = 1:numel(data)
            if thresh > 0
                ind = find(dataThreshold{iT} > thresh, 1, 'first');
            else
                ind = find(dataThreshold{iT} < thresh, 1, 'first');
            end
            if ~isempty(ind)
                startIndex(iT) = ind;
            end
        end
        
    else
        startIndex = nan(size(data, 1), 1);

        for iT = 1:size(dataThreshold, 1)
            if thresh > 0
                ind = find(dataThreshold(iT, :) > thresh, 1, 'first');
            else
                ind = find(dataThreshold(iT, :) < thresh, 1, 'first');
            end
            if ~isempty(ind)
                startIndex(iT) = ind;
            end
        end
    end

end