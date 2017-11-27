function dataHP = highPassFilter(data, Fs, varargin)
% data is either:
%   - nTrials x 1 cell of nTime x nChannels matrices
%   - nTrials x nTime x nChannels tensor

    p = inputParser();
    p.addParameter('cornerHz', 250, @isscalar);
    p.addParameter('order', 4, @isscalar);
    % if true subtracts first sample from each signal which can help reduce onset transients for IIR filters
    p.addParameter('subtractFirstSample', true, @islogical);
    % if true, the first sample will be added back on afterwards
    p.addParameter('addBackFirstSample', false, @islogical);
    p.addParameter('filtfilt', false, @islogical);
    p.addParameter('showProgress', false, @islogical);
    p.parse(varargin{:});
    
    if p.Results.cornerHz == 0
        dataHP = data;
        return;
    end
    
    hpCornerNormalized = p.Results.cornerHz / (Fs/2);
    [B, A] = butter(p.Results.order, hpCornerNormalized, 'high');
    
    if iscell(data)
        timeDim = 1;
    else
        timeDim = 2;
    end
    dataHP = ERAASR.Utils.filterIgnoreLeadingTrailingNaNs(B, A, data, ...
        'dim', timeDim, 'filtfilt', p.Results.filtfilt, 'showProgress', p.Results.showProgress, ...
        'subtractFirstSample', p.Results.subtractFirstSample, 'addBackFirstSample', p.Results.addBackFirstSample);
    
end