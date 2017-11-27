function rms = computeRMS(data, varargin)
% data is either:
%  - nTrials x 1 cell of time x channels matrices
%  - nTrials x nTime x channels tensor
% rms is nTrials x nChannels matrix of RMS values. if perTrial is set to
% false, all rows of rms will be identical

    p = inputParser();
    p.addParameter('perTrial', true, @islogical);
    p.addParameter('smoothOverTrials', 0, @isscalar);
    p.addParameter('clip', [], @(x) numel(x) <= 2); % +/- this value will be ignored in the computation
    p.parse(varargin{:});

    if ~isempty(p.Results.clip)
        data = TensorUtils.clip(data, p.Results.clip, NaN);
    end

    if iscell(data)
        temp = cellfun(@(x) nansum((x-nanmean(x, 1)).^2, 1), data, 'UniformOutput', false);
        ssqByTrial = cat(1, temp{:});
        temp = cellfun(@(x) sum(~isnan(x), 1), data, 'UniformOutput', false);
        countByTrial = cat(1, temp{:});
    else
        ssqByTrial = ERAASR.TensorUtils.squeezeDims(nansum((data-nanmean(data, 2)).^2, 2), 2);
        countByTrial = ERAASR.TensorUtils.squeezeDims(sum(~isnan(data), 2), 2);
    end

    if p.Results.perTrial
        rms = sqrt(ssqByTrial ./ countByTrial);

        if p.Results.smoothOverTrials > 1
            for iC = 1:size(rms, 2)
                rms(:, iC) = smooth(rms(:, iC), p.Results.smoothOverTrials);
                % smooth leaves the edges a bit rough
                rms(1, iC) = rms(2, iC);
                rms(end, iC) = rms(end-1, iC);
            end
        end
    else
        rms = sqrt(nansum(ssqByTrial, 1) ./ nansum(countByTrial, 1));
%         rms = repmat(rms, size(ssqByTrial, 1), 1);
    end

end