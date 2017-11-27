 function [tensorAligned, extractIndStart, extractIndStop] = preciseAlignExtractArtifacts(data, startIndexGuess, varargin)
% Extract and precisely align over trials a segment of data
% data is either:
%   - nTrials x 1 cell of nTime x nChannels matrices
%   - nTrials x nTime x nChannels tensor
% startIndexGuess is a guess at the start index for each trial, which will
%   then be precisely aligned
% tensorAligned will be nTrials x extractWindowDuration x nChannels tensor

    p = inputParser();
    p.addRequired('data', @(x) iscell(x) || isnumeric(x));
    p.addRequired('startIndexGuess', @(x) isnumeric(x) && isvector(x));
    p.addOptional('alignChannel', 1, @isscalar);
    p.addOptional('alignWindowPre', 15, @isscalar);
    p.addOptional('alignWindowDuration', 360, @isscalar);
    p.addOptional('alignUpsampleBy', 10, @isscalar);

    p.addOptional('extractWindowPre', 15, @isscalar);
    p.addOptional('extractWindowDuration', 360, @isscalar);
    
    p.addParameter('quiet', false, @islogical);
    p.parse(data, startIndexGuess, varargin{:});

    %% Extract a time window from each trial around the approximate start time that generously encompasses the full artifact

    % extract nTrials x alignWindowDuration x nChannels tensor around each
    % trial's startIndexGuess, taking alignWindowPre sampels before and a
    % total of alignWindowDuration samples
    tensorAlignWindow = ERAASR.extractTensorFromTrialCellArray(data, startIndexGuess - p.Results.alignWindowPre, p.Results.alignWindowDuration, p.Results.alignChannel);
    
    % Determine the precise sub-sample alignment for a subset of channels
    delays = ERAASR.Utils.findDelaysTensor(tensorAlignWindow, ...
        'maxLag', p.Results.alignWindowPre, 'timeDimension', 2, 'simultaneousDimensions', 3, ...
        'alignToEachOtherDimensions', 1, 'upsample', p.Results.alignUpsampleBy, ...
        'quiet', p.Results.quiet);

    %% Extract cleaning window based on trial-to-trial delays just computed

    indStartClean = -p.Results.extractWindowPre + startIndexGuess + floor(delays);
    indStopClean = indStartClean + p.Results.extractWindowDuration -1;
    remainingDelays = delays - floor(delays);

    % we'll generate time vectors in sample units for each trial so that we can
    % do the per-trial resampling
    
    if iscell(data)
        nTrials = numel(data);
        dataClean = deal(cell(nTrials, 1));

        for iTrial = 1:nTrials
            dataClean{iTrial} = data{iTrial}(indStartClean(iTrial):indStopClean(iTrial), :);
        end
        % nTrials x nTime x nChannels
        dataTensor = permute(cat(3, dataClean{:}), [3 1 2]);
    else
        nTrials = size(data, 1);
        nChannels = size(data, 3);
        dataTensor = nan(nTrials, p.Results.extractWindowDuration, nChannels);
        
        for iTrial = 1:nTrials
            dataTensor(iTrial, :, :) = data(iTrial, indStartClean(iTrial):indStopClean(iTrial), :);
        end
    end

    %% Remove remaining subsample delays by resampling
    tensorAligned = ERAASR.Utils.removeDelaysTensor(dataTensor, remainingDelays, 'fillMode', 'hold', ...
        'alignToEachOtherDimensions', 1, 'timeDimension', 2, 'quiet', p.Results.quiet, ...
        'simultaneousDimensions', 3, 'upsample', p.Results.alignUpsampleBy);

    extractIndStart = indStartClean;
    extractIndStop = indStopClean;
end