function [dataTensorCleaned, extract] = cleanArtifactTensor(dataTensor, varargin)
% cleanedTensor = cleanArtifactTensor(dataTensor, varargin)
% dataTensor is nTrials x extractWindowDuration x nChannels

    p = inputParser();
    p.addRequired('Fs', @isscalar);
    p.addParameter('hpCornerHz', 200, @isscalar); % light high pass filtering only
    p.addParameter('hpOrder', 4, @isscalar); % order of high-pass filter
    p.addParameter('upsampleBy', 1, @isscalar);
    p.addParameter('samplesPre', 0, @isscalar);
    p.addParameter('samplesPerPulse', 30, @isscalar);
    p.addParameter('nPulses', 20, @isscalar);

    p.addParameter('nPC_channels', 12, @isscalar);
    p.addParameter('nPC_trials', 2, @isscalar);
    p.addParameter('nPC_pulses', 6, @isscalar);

    p.addParameter('omit_bandwidth_channels', 3, @isscalar);
    p.addParameter('omit_bandwidth_trials', 1, @isscalar);
    p.addParameter('omit_bandwidth_pulses', 1, @isscalar);

    p.addParameter('pcaOnlyOmitted', true, @islogical); % thesis used this as false but true works better
    p.addParameter('cleanOverChannelsIndividualTrials', false, @islogical);
    p.addParameter('cleanOverPulsesIndividualChannels', false, @islogical);
    p.addParameter('cleanOverTrialsIndividualChannels', true, @islogical);
    
    p.addParameter('cleanPostStim', true, @islogical);

    p.addParameter('alignPulsesOverTrain', false, @islogical);

    p.addParameter('showFigures', true, @islogical);
    p.addParameter('plotTrials', 1, @isvector);
    p.addParameter('plotPulses', 1, @isvector);
    p.addParameter('figurePath', pwd, @ischar);
    p.addParameter('saveFigures', false, @islogical);
    p.addParameter('saveFigureCommand', @(filepath) save([filepath '.png']));
    
    p.addParameter('quiet', false, @islogical);
    p.parse(varargin{:});

    Fs = p.Results.Fs;
    nPC_channels = p.Results.nPC_channels;
    nPC_trials = p.Results.nPC_trials;
    nPC_pulses = p.Results.nPC_pulses;

    obw_channels = p.Results.omit_bandwidth_channels;
    obw_trials = p.Results.omit_bandwidth_trials;
    obw_pulses = p.Results.omit_bandwidth_pulses;
    
    plotTrialIdx = p.Results.plotTrials;
    plotPulseIdx = p.Results.plotPulses;
    
    figureIdx = 1;
    showFigures = p.Results.showFigures;
    figurePath = p.Results.figurePath;
    saveFigures = p.Results.saveFigures;
    saveFigureCommand = p.Results.saveFigureCommand;
    quiet = p.Results.quiet;
    upsample = p.Results.upsampleBy;
    samplesPre = p.Results.samplesPre;
    pulseLen = p.Results.samplesPerPulse * upsample;
    nPulses = p.Results.nPulses;
    totalTrainSamples = pulseLen * nPulses;
    

    dataTensorCleaned = dataTensor;
    
    if upsample > 1
        tvec = 1:size(dataTensor, 2);
        dataTensor = ERAASR.Utils.resampleTensorInTime(dataTensor, 2, tvec, 'timeDelta', 1/upsample);
    end

    % very light filtering of stim trials before artifact computation
    if p.Results.hpCornerHz > 0
        if ~quiet, fprintf('Filtering before performing artifact subtraction\n'); end
        dataTensor = ERAASR.highPassFilter(dataTensor, Fs, 'cornerHz', p.Results.hpCornerHz, ...
            'order', p.Results.hpOrder, 'filtfilt', true, 'subtractFirstSample', true);
    end

    % segment stim trials into individual pulses into segment tensor
    nTraces = size(dataTensor, 1);
    nChannels = size(dataTensor, 3);
    
    if samplesPre + totalTrainSamples > size(dataTensor, 2)
        error('dataTensor does not have enough timepoints for samplesPre + totalTrainSamples');
    end
    
    dataTensor_idxTime = samplesPre + (1:totalTrainSamples);
    segmentTensorRaw = reshape(dataTensor(:, dataTensor_idxTime, :), [nTraces, pulseLen, nPulses, nChannels]);

    if showFigures
        for iR = 1:numel(plotTrialIdx)
            for iP = 1:numel(plotPulseIdx)
                figSetName('raw_overChannels_trial%d_pulse%d', plotTrialIdx(iR), plotPulseIdx(iP));
                ERAASR.PlotUtils.plotOverChannels(segmentTensorRaw, 1, 1);
                saveFig();
            end
        end
    end

    %%%%%%%
    % Clean over channels using PCA technique
    %%%%%%%
   
    segmentTensorCached = segmentTensorRaw;

    % subtract mean over time points
    segmentTensor = bsxfun(@minus, segmentTensorRaw, nanmean(segmentTensorRaw, 2));
    extract.segmentTensor = segmentTensor; %#ok<*AGROW>

    if nPC_channels > 0
        if ~quiet, fprintf('Cleaning over channels using %d PCs\n', nPC_channels); end

        if p.Results.cleanOverChannelsIndividualTrials
            % originally R x T x P x C
            % TP x C x R --> each PC is TP x 1 for each R
            pcaMatOverTrials = ERAASR.TensorUtils.reshapeByConcatenatingDims(segmentTensor, {[2 3], 4, 1});
            pcaArt = nan(size(pcaMatOverTrials));

            pcaCleaned = nan(size(pcaMatOverTrials));
            
            if ~quiet, prog = ERAASR.Utils.ProgressBar(nTraces, 'Cleaning over channels for each individual trial'); end
            for r = 1:nTraces
                if ~quiet, prog.update(r); end
                pcaMat = pcaMatOverTrials(:, :, r); % TP x C
                [pcaCleaned(:, :, r), artPcs, pcaArt(:, :, r)] = ERAASR.cleanMatrixViaPCARegression(pcaMat, nPC_channels, ...
                    'omitAdjacentChannelsBandWidth', obw_channels, 'pcaOnlyOmitted', p.Results.pcaOnlyOmitted);
            end
            if ~quiet, prog.finish(); end
            segmentTensor = ERAASR.TensorUtils.undoReshapeByConcatenatingDims(pcaCleaned, {[2 3], 4, 1}, size(segmentTensor));

            % reconstruct artifact detection into segmentTensor shape
            tensorArtifact = ERAASR.TensorUtils.undoReshapeByConcatenatingDims(pcaArt, {[2 3], 4, 1}, size(segmentTensor));
        else
            % clean using PCA + regression
            % RTP x C --> each PC is RTP x 1
            pcaMat = ERAASR.TensorUtils.reshapeByConcatenatingDims(segmentTensor, {[2 3 1], 4});
            [pcaCleaned, artPcs, artMat] = ERAASR.cleanMatrixViaPCARegression(pcaMat, nPC_channels, ...
                'omitAdjacentChannelsBandWidth', obw_channels, 'pcaOnlyOmitted', p.Results.pcaOnlyOmitted);
            segmentTensor = ERAASR.TensorUtils.undoReshapeByConcatenatingDims(pcaCleaned, {[2 3 1], 4}, size(segmentTensor));

            % reconstruct artifact detection into segmentTensor shape
            tensorArtifact = ERAASR.TensorUtils.undoReshapeByConcatenatingDims(artMat, {[2 3 1], 4}, size(segmentTensor));
        end

        % subtract mean over time points
        segmentTensor = bsxfun(@minus, segmentTensor, nanmean(segmentTensor, 2));

        extract.pcaArt_channels = artPcs;
        extract.tensorArtifact_cleanedChannels = tensorArtifact;
        extract.segmentTensor_cleanedChannels = segmentTensor;

        if showFigures
            % show stacked pcs over channels
            figSetName('channels_artifact_PCs');
            pcsFirstTrial = artPcs(1:pulseLen*nPulses, :);
            ERAASR.PlotUtils.plotStackedPCTraces(pcsFirstTrial(:, 1:min(6, size(pcsFirstTrial,2))));
            saveFig();

            for iR = 1:numel(plotTrialIdx)
                figSetName('channels_prePost_overPulses_trial%d', plotTrialIdx(iR));
                ERAASR.PlotUtils.plotBeforeAfterOverTrainAllChannels(segmentTensorCached, segmentTensor, plotTrialIdx(iR));
                saveFig();

                figSetName('channels_cleaned_overPulses_trial%d', plotTrialIdx(iR))
                ERAASR.PlotUtils.plotOverTrainAllChannels(segmentTensor, plotTrialIdx(iR));
                saveFig();
            end

            for iP = 1:numel(plotPulseIdx)
                figSetName('channels_artifacts_overTrials', plotPulseIdx(iP))
                ERAASR.PlotUtils.plotRawWithArtifactOverTrialsAllChannels(segmentTensorCached, tensorArtifact, plotPulseIdx(iP));
                saveFig();
            end
        end
    end

    %%%%%%%
    % Clean over pulses using PCA
    %%%%%%%

    segmentTensorCache = segmentTensor;

    if p.Results.alignPulsesOverTrain
        % align pulses over the train to each other using delays computed
        % from the raw uncleaned segmentTensor
        delays = ERAASR.Utils.findDelaysTensor(segmentTensorRaw, ...
            'alignTo', 'first', 'quiet', quiet, ...
            'timeDimension', 2, 'simultaneousDimensions', 4, 'alignToEachOtherDimensions', 3, ...
            'upsample', 10, 'message', 'Aligning each pulse to others in train');

        segmentTensor = ERAASR.Utils.removeDelaysTensor(segmentTensor, delays, ...
            'fillMode', 'hold', 'quiet', quiet, ...
            'timeDimension', 2, 'simultaneousDimensions', 4, 'alignToEachOtherDimensions', 3, ...
            'upsample', 10, 'message',  'Applying alignments to cleaned pulses');
    end

    % clean by deleting top modes of PCA
    % we only do cleaning on some timepoints within the pulse and then
    % splice the cleaned data back in (these are the timepoints where
    % the pulse occurs)

    segmentTensorCached = segmentTensor;
    % in case we only want to look at the inital part of each pulse
    timeInPulseIdx = 1:pulseLen;
    segmentTensorPartial = segmentTensor(:, timeInPulseIdx, :, :);

    if nPC_pulses > 0
        if ~quiet, fprintf('Cleaning over pulses using %d PCs\n', nPC_pulses); end
        if p.Results.cleanOverPulsesIndividualChannels
            % originally R x T x P x C
            % to TR x P x C --> each PC is TR x 1 for each C
            pcaMatOverChannels = ERAASR.TensorUtils.reshapeByConcatenatingDims(segmentTensor, {[2 1], 3, 4});
            pcaArt = nan(size(pcaMatOverChannels));

            pcaCleaned = nan(size(pcaMatOverChannels));
            
            if ~quiet, prog = ERAASR.Utils.ProgressBar(nChannels, 'Cleaning over pulses for each individual channel'); end
            for c = 1:nChannels
                if ~quiet, prog.update(r); end
                pcaMat = pcaMatOverChannels(:, :, c); % TR x 1, weights are P x 1
                [pcaCleaned(:, :, c), ~, pcaArt(:, :, c)] = ERAASR.cleanMatrixViaPCARegression(pcaMat, nPC_pulses, ...
                    'omitAdjacentChannelsBandWidth', obw_pulses, 'pcaOnlyOmitted', p.Results.pcaOnlyOmitted);
            end
            if ~quiet, prog.finish(); end
            
            segmentTensorPartial = ERAASR.TensorUtils.undoReshapeByConcatenatingDims(pcaCleaned, {[2 1], 3, 4}, size(segmentTensorPartial));

            % reconstruct artifact detection into segmentTensor shape
            tensorArtifact = ERAASR.TensorUtils.undoReshapeByConcatenatingDims(pcaArt, {[2 1], 3, 4}, size(segmentTensorPartial));

        else
            % TRC x P --> each PC is TRC x 1, weights are P x 1,
            pcaMat = ERAASR.TensorUtils.reshapeByConcatenatingDims(segmentTensorPartial, {[2 1 4], 3});
            [pcaCleaned, ~, pcaArt] = ERAASR.cleanMatrixViaPCARegression(pcaMat, nPC_pulses, ...
                'omitAdjacentChannelsBandWidth', obw_pulses, 'pcaOnlyOmitted', p.Results.pcaOnlyOmitted);

            segmentTensorPartial = ERAASR.TensorUtils.undoReshapeByConcatenatingDims(pcaCleaned, {[2 1 4], 3}, size(segmentTensorPartial));

            tensorArtifact = ERAASR.TensorUtils.undoReshapeByConcatenatingDims(pcaArt, {[2 1 4], 3}, size(segmentTensorPartial));
        end

        % store preserving continuity to the right
        if max(timeInPulseIdx) == pulseLen
            segmentTensor = segmentTensorPartial;
        else
            spliceStart = timeInPulseIdx(end)+1;
            oldDelta = segmentTensor(:, spliceStart, :, :) - segmentTensor(:, spliceStart-1, :, :);
            segmentTensor(:, timeInPulseIdx, :, :) = segmentTensorPartial;
            currentDelta = segmentTensor(:, spliceStart, :, :) - segmentTensor(:, spliceStart-1, :, :);
            segmentTensor(:, spliceStart:end, :, :) = bsxfun(@minus, segmentTensor(:, spliceStart:end, :, :), currentDelta - oldDelta);
        end

        % transform back to pca mat
        % each trial is corrupted by a sum of artifacts
        % clean by deleting top 3 modes of PCA
        % TP x R x 1 for each C--> each PC is TPC x 1
        segmentTensor = bsxfun(@minus, segmentTensor, nanmean(segmentTensor, 2));

        extract.pcaArt_pulses = pcaArt;
        extract.tensorArtifact_cleanedPulses = tensorArtifact;
        extract.segmentTensor_cleanedPulses = segmentTensor;

        if showFigures
            for iR = 1:numel(plotTrialIdx)
                figSetName('pulses_artifacts_overPulses_trial%d', plotTrialIdx(iR));
                ERAASR.PlotUtils.plotRawWithArtifactOverPulsesAllChannels(segmentTensorCache, tensorArtifact, plotTrialIdx(iR));
                saveFig();

                figSetName('pulses_prePost_overPulses_trial%d', plotTrialIdx(iR));
                ERAASR.PlotUtils.plotBeforeAfterOverTrainAllChannels(segmentTensorCached, segmentTensor, plotTrialIdx(iR));
                saveFig();

                figSetName('pulses_cleaned_overPulses_trial%d', plotTrialIdx(iR));
                ERAASR.PlotUtils.plotOverTrainAllChannels(segmentTensor, plotTrialIdx(iR));
                saveFig();
            end

            for iP = 1:numel(plotPulseIdx)
                figSetName('pulses_overTrials_pulse%d', plotPulseIdx(iP));
                ERAASR.PlotUtils.plotOverTrialsAllChannels(segmentTensor, plotPulseIdx(iP));
                saveFig();
            end
        end
    end

    %%%%%%%
    % Clean over trials using PCA
    %%%%%%%

    if nPC_trials > 0
        if ~quiet, fprintf('Cleaning over trials using %d PCs\n', nPC_trials); end
        
        segmentTensorCached = segmentTensor;

        if p.Results.cleanOverTrialsIndividualChannels
            % TP x R x CPP
            pcaMatOverChannels = ERAASR.TensorUtils.reshapeByConcatenatingDims(segmentTensor, {[2 3], 1, 4});
            pcaArt = nan(size(pcaMatOverChannels));

            pcaCleaned = nan(size(pcaMatOverChannels));
            if ~quiet, prog = ERAASR.Utils.ProgressBar(nChannels, 'Cleaning over trials for each individual channel'); end
            for c = 1:nChannels
                if ~quiet, prog.update(c); end
                pcaMat = pcaMatOverChannels(:, :, c); % TP x R
                [pcaCleaned(:, :, c), ~, pcaArt(:, :, c)] = ERAASR.cleanMatrixViaPCARegression(pcaMat, nPC_trials, ...
                    'omitAdjacentChannelsBandWidth', obw_trials, 'pcaOnlyOmitted', p.Results.pcaOnlyOmitted);
            end
            if ~quiet, prog.finish(); end

            segmentTensor = ERAASR.TensorUtils.undoReshapeByConcatenatingDims(pcaCleaned, {[2 3], 1, 4}, size(segmentTensor));
            
            % reconstruct artifact into segmentTensor shape
            tensorArtifact = ERAASR.TensorUtils.undoReshapeByConcatenatingDims(pcaArt, {[2 3], 1, 4}, size(segmentTensor));

        else
            % originally R x T x P x C
            % TPC x R --> each PC is TPC x 1, weights are R x 1,
            pcaMat = ERAASR.TensorUtils.reshapeByConcatenatingDims(segmentTensor, {[2 3 4], 1});
            [pcaCleaned, ~, pcaArt] = ERAASR.cleanMatrixViaPCARegression(pcaMat, nPC_trials, ...
                'omitAdjacentChannelsBandWidth', obw_trials, 'pcaOnlyOmitted', p.Results.pcaOnlyOmitted);

            segmentTensor = ERAASR.TensorUtils.undoReshapeByConcatenatingDims(pcaCleaned, {[2 3 4], 1}, size(segmentTensor));
            
            % reconstruct artifact into segmentTensor shape
            tensorArtifact = ERAASR.TensorUtils.undoReshapeByConcatenatingDims(pcaArt, {[2 3 4], 1}, size(segmentTensor));
        end

        % subtract mean over timepoints again
        segmentTensor = bsxfun(@minus, segmentTensor, nanmean(segmentTensor, 2));
  
        extract.pcaArt_trials = pcaArt;
        extract.tensorArtifact_cleanedTrials = tensorArtifact;
        extract.segmentTensor_cleanedTrials = segmentTensor;
        
        extract.tensorArtifact_total = segmentTensorRaw - segmentTensor;

        if showFigures
            for iP = 1:numel(plotPulseIdx)
                figSetName('trials_artifacts_overTrials_pulse%d', plotPulseIdx(iP));
                ERAASR.PlotUtils.plotRawWithArtifactOverTrialsAllChannels(segmentTensorCached, tensorArtifact, 1);
                saveFig();

                figSetName('trials_beforeAfter_overTrials_pulse%d', plotPulseIdx(iP));
                ERAASR.PlotUtils.plotBeforeAfterOverTrialsAllChannels(segmentTensorCached, segmentTensor, plotPulseIdx(iP));
                saveFig();

                figSetName('trials_cleaned_overTrials_pulse%d', plotPulseIdx(iP));
                ERAASR.PlotUtils.plotOverTrialsAllChannels(segmentTensor, 1);
                saveFig();
            end

            for iR = 1:numel(plotTrialIdx(iR))
                figSetName('trials_cleaned_overPulses_trial%d', plotTrialIdx(iR));
                ERAASR.PlotUtils.plotOverTrainAllChannels(segmentTensor, plotTrialIdx(iR));
                saveFig();
            end
        end
    end

    % Rebuild timeseries

    % adjust adjacent pieces of the tensor to preserve continuity
    % N x T x P x C
    segmentTensorCont = segmentTensor;
    for iPulse = 1:nPulses
        if iPulse == 1
            % restore first timepoint to previous value
            segmentTensorCont(:, :, 1, :) = bsxfun(@plus, segmentTensorCont(:, :, 1, :), segmentTensorRaw(:, 1, 1, :) - segmentTensorCont(:, 1, 1, :));
        else
            currentDiffFromLast = segmentTensorCont(:, 1, iPulse, :) - segmentTensorCont(:, end, iPulse-1, :);
            segmentTensorCont(:, :, iPulse, :) = bsxfun(@minus, segmentTensorCont(:, :, iPulse, :), currentDiffFromLast);
        end
    end
    segmentTensorCleaned = segmentTensorCont;
    insert = reshape(segmentTensorCleaned, [nTraces, pulseLen * nPulses, nChannels]);
    
    % reinsert artifact subtracted pieces into original signal, bias post-stim data to preserve continuity
    dataTensor(:, dataTensor_idxTime, :) = insert;
    
    % clean the data post stim using PCA over channels
    if nPC_channels > 0 && p.Results.cleanPostStim
        if ~quiet
            fprintf('Cleaning post-stim transient over channels using %d PCs\n', nPC_channels);
        end
        
        dataPostStim = dataTensor(:, (dataTensor_idxTime(end)+1) : end, :);
        
        % clean the data following the last pulse using  PC again
        % subtract mean over time (R x T x C)
        meanEachTrialChannel = nanmean(dataPostStim, 2);
        dataPostStim = bsxfun(@minus, dataPostStim, meanEachTrialChannel);

        % RT x C
        pcaMat = ERAASR.TensorUtils.reshapeByConcatenatingDims(dataPostStim, {[1 2] 3});
        pcaCleaned = ERAASR.cleanMatrixViaPCARegression(pcaMat, nPC_channels, ...
            'omitAdjacentChannelsBandWidth', obw_channels, 'pcaOnlyOmitted', p.Results.pcaOnlyOmitted);
        dataPostStim = ERAASR.TensorUtils.undoReshapeByConcatenatingDims(pcaCleaned, {[1 2], 3}, size(dataPostStim));

        % subtract mean over trials
        dataPostStim = bsxfun(@minus, dataPostStim, nanmean(dataPostStim, 1));

        delta = bsxfun(@minus, dataPostStim(:, 1, :), dataTensor(:, dataTensor_idxTime(end), :));
        dataPostStim = bsxfun(@minus, dataPostStim, delta);

        dataTensor(:, (dataTensor_idxTime(end)+1) : end, :) = dataPostStim;
    end
    
    % downsample dataTensor if it was upsampled
    if upsample > 1
        if ~quiet
            fprintf('Downsampling dataTensor\n');
        end
        tvec = (1:size(dataTensor, 2)) / upsample;
        dataTensor = ERAASR.Utils.resampleTensorInTime(dataTensor, 2, tvec, 'timeDelta', 1);
    end
    
    % reinsert the cleaned data into the raw, unfiltered data tensor, since
    % acausal filtering can add a transient to the pre-cleaned data
    % maintain continuity at the insert point
    delta = bsxfun(@minus, dataTensor(:, samplesPre+1, :), dataTensorCleaned(:, samplesPre, :));
    dataTensor = bsxfun(@minus, dataTensor, delta);
    
    dataTensorRaw = dataTensorCleaned;
    dataTensorCleaned(:, (samplesPre+1) : end, :) = dataTensor(:, (samplesPre+1) : end, :);
    extract.artifact_total = dataTensorRaw - dataTensorCleaned;
    
    function figSetName(varargin)
        name = sprintf(varargin{:});
        figh = figure(figureIdx);
        figh.Name = name;
        figh.NumberTitle = 'off';
        figureIdx = figureIdx + 1;
    end

    function saveFig()
        if saveFigures
            figh = gcf;
            figfile = fullfile(figurePath, figh.Name);
            saveFigureCommand(figfile);
        end
    end

end




