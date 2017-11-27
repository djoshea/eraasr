function [spikeTimes, waveforms] = extractSpikesCrossingThreshold(data, threshold, varargin)
% data is either:
%   nTrials x time x channels tensor
%   nTrials x 1 cell of nTime x nChannels

    p = inputParser();
    p.addParameter('mode', 'largestFirst', @ischar); % largestFirst or leftToRight
    p.addParameter('waveformSamplesPrePost', [10 38], @isvector); % should add to total number of waveform samples
    p.addParameter('lockoutPrePost', [], @(x) isempty(x) || isvector(x));
    p.addParameter('extremumWithin', 7, @islogical);
    p.addParameter('quiet', false, @islogical);
    p.parse(varargin{:});
    quiet = p.Results.quiet;
    snippetPre = p.Results.waveformSamplesPrePost(1);
    snippetPost = p.Results.waveformSamplesPrePost(2);
    snippetLength = snippetPre + snippetPost;

    if isempty(p.Results.lockoutPrePost)
        lockoutPre = samplesPre;
        lockoutPost = samplesPost;
    else
        lockoutPre = p.Results.lockoutPrePost(1);
        lockoutPost = p.Results.lockoutPrePost(2);
    end
    extremumWithin = p.Results.extremumWithin;
    
    if iscell(data)
        nTrials = numel(data);
        nChannels = size(data{1}, 2);
    else
        nTrials = size(data, 1);
        nChannels = size(data, 3);
    end
    
    threshold = ERAASR.TensorUtils.singletonExpandToSize(threshold, [nTrials nChannels]);

    % utility for finding all crossing times of thresh
    function idx = findCrossings(dataVec, thresh)
        if thresh < 0
            % falling threshold < 0
            idx = find(diff(dataVec <= thresh) == 1);
        else
            % rising threshold > 0
            idx = find(diff(dataVec >= thresh) == 1);
        end
    end

    % identify threshold crossings
    crossings = cell(nTrials, nChannels);
    if iscell(data)
        for iR = 1:nTrials
            for iC = 1:nChannels
                crossings{iR, iC} = findCrossings(data{iR}(:, iC), threshold(iR, iC));
            end
        end
    else
        for iR = 1:nTrials
            for iC = 1:nChannels
                crossings{iR, iC} = findCrossings(data(iR, :, iC), threshold(iR, iC));
            end
        end
    end

    % ensure they are spaced by at least 1 waveform
    % loop through and extract snippets
    [spikeTimes, waveforms] = deal(cell(nTrials, nChannels));
    
    if strcmp(p.Results.mode, 'leftToRight')
        % sweep left to right, thresholding as you go
        if ~quiet, prog = ERAASR.Utils.ProgressBar(nTrials, 'Extracting threhold crossings from continuous data'); end
        for iR = 1:nTrials
            if ~quiet, prog.update(iR); end
            for iC = 1:nChannels
                lastCrossing = -Inf;
                nCrossingsKept = 0;
                W = numel(crossings{iR, iC});
                mask = true(W, 1);
                waveMat = nan(W, snippetLength);

                if iscell(data)
                    nTimeThisTrial = size(data{iR},1);
                else
                    nTimeThisTrial = size(data, 2);
                end

                for iW = 1:W
                    thisCross = crossings{iR,iC}(iW);
                    if thisCross - lastCrossing < lockoutPost || thisCross <= snippetPre || thisCross + snippetPost - 1 > nTimeThisTrial
                        % too close to lastCrossing or to ends of data trace --> throw away
                        mask(iW) = false;
                    else
                        % keep this crossing, sample snippet
                        lastCrossing = thisCross;
                        nCrossingsKept = nCrossingsKept + 1;
                        if iscell(data)
                            waveMat(nCrossingsKept, :) = data{iR}(thisCross + (-snippetPre : snippetPost-1), iC)';
                        else
                            waveMat(nCrossingsKept, :) = data(iR, thisCross + (-snippetPre : snippetPost-1), iC)';
                        end
                    end
                end

                spikeTimes{iR, iC} = crossings{iR, iC}(mask);
                waveforms{iR, iC} = waveMat(1:nCrossingsKept, :);
            end
        end
        if ~quiet, prog.finish(); end

    elseif strcmp(p.Results.mode, 'largestFirst')
        % take biggest waveforms first, and ignore those within lockout
        % window
        if ~quiet, prog = ERAASR.Utils.ProgressBar(nTrials, 'Extracting threhold crossings from continuous data'); end
        for iR = 1:nTrials
            if ~quiet, prog.update(iR); end
            
            for iC = 1:nChannels
                thresh = threshold(iR, iC);
                cross = crossings{iR, iC};
                W = numel(cross);
                if iscell(data)
                    nTimeThisTrial = size(data{iR}, 1);
                else
                    nTimeThisTrial = size(data, 2);
                end
                waveMat = nan(W, snippetLength);

                % figure out the extreme value of each snipet
                waveExt = nan(W, 1);
                for iW = 1:W
                    thisCross = cross(iW);
                    if thisCross - snippetPre < 1, continue, end
                    if thisCross + snippetPost - 1 > nTimeThisTrial, continue; end
                    if iscell(data)
                        snip = data{iR}(thisCross+(0 : extremumWithin-1), iC);
                        waveMat(iW, :) = data{iR}(thisCross + (-snippetPre : snippetPost-1), iC)';
                    else
                        snip = data(iR, thisCross+(0 : extremumWithin-1), iC)';
                        waveMat(iW, :) = data(iR, thisCross + (-snippetPre : snippetPost-1), iC)';
                    end
                    if thresh > 0
                        waveExt(iW) = max(snip);
                    else
                        waveExt(iW) = min(snip);
                    end
                end

                maskPicked = false(W, 1);
                maskEligible = true(W, 1);

                while(any(maskEligible))
                    % waveExt(ineligible) will be NaN
                    [~, pick] = max(abs(waveExt));

                    if isnan(waveExt(pick))
                        break;
                    end
                    maskPicked(pick) = true;
                    maskEligible(pick) = false;

                    % mark ineligible other crossings within lockout window
                    tooClose = cross - cross(pick) >= -lockoutPre & cross-cross(pick) < lockoutPost;
                    maskEligible(tooClose) = false;
                    waveExt(tooClose) = NaN;
                end

                spikeTimes{iR, iC} = crossings{iR, iC}(maskPicked);
                waveforms{iR, iC} = waveMat(maskPicked, :);
            end
        end
        if ~quiet, prog.finish(); end
    else
        error('Unknown mode');
    end

end
