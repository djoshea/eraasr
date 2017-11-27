classdef Parameters
    properties
        Fs % sampling rate in Hz
        
        thresholdChannel = 1; % threshold which channel
        thresholdValue % threshold channel at this value to find approximate start 
        thresholdHPCornerHz = 250; % High pass filter at this corner frequency before thresholding
        
        % For alignment
        alignChannel = 1;
        alignUpsampleBy = 10; % supersample by this ratio before alignment
        alignWindowPre
        alignWindowDuration
        
        %% For Cleaning
        
        extractWindowPre % number of samples to extract before the threshold crossing, should be sufficient to allow HP filtering 
        extractWindowDuration % total number of samples to extract, which includes the pre, during stim, and post stim windows
        cleanStartSamplesPreThreshold % number of samples before threshold crossing to include within the first pulse
        
        cleanHPCornerHz = 200 % light high pass filtering at the start of cleaning
        cleanHPOrder = 4; % high pass filter order 
        cleanUpsampleBy % upsample by this ratio during cleaning
        samplesPerPulse
        nPulses

        nPC_channels = 12;
        nPC_trials = 2;
        nPC_pulses = 6;

        % when reconstructing each, omit this number of *TOTAL*
        % channels/trials/pulses, including the one being reconstructed. 
        % So 3 means omit me and my immediate neighbors
        omit_bandwidth_channels = 3;
        omit_bandwidth_trials = 1;
        omit_bandwidth_pulses = 1;

        alignPulsesOverTrain logical = false; % do a secondary alignment within each train, in case you think there is pulse to pulse jitter. Works best with upsampling
        
        pcaOnlyOmitted logical = true; % if true, build PCs only from non-omitted channels/trials/pulses. if false, build PCs from all but set coefficients to zero post hoc
        
        cleanOverChannelsIndividualTrials logical = false;
        cleanOverPulsesIndividualChannels logical = false;
        cleanOverTrialsIndividualChannels logical = true;
    
        cleanPostStim logical = true; % clean the post stim window using a single PCR over channels

        showFigures logical = false; % useful for debugging and seeing well how the cleaning works
        plotTrials = 1; % which trials to plot in figures, can be vector
        plotPulses = 1; % which pulses to plot in figures, can be vector
        figurePath % folder to save the figures
        saveFigures logical = false; % whether to save the figures
        saveFigureCommand = @(filepath) print('-dpng', '-r300', [filepath '.png']); % specify a custom command to save the figure
        
        quiet logical = false; % if true, don't print anything
    end
   
    methods
        function opts = Parameters()
            opts.figurePath = pwd;
        end
        
        function check(opts)
            props = properties(opts);
            missing = cellfun(@(prop) isempty(opts.(prop)), props);
            
            if any(missing)
                error('Missing Parameter settings for %s', strjoin(props(missing), ', '));
            end
        end
    end
end