function [delays, aligned] = findDelaysTensor(data,varargin)
% Estimates delays between sets of simultaneously sampled signals.
%   delays = finddelaysMultiAnalogTensor(x,varargin)
%
%   inputs: 
%     data: numeric, can be arbitrary size and ordering of dimensions is
%     arbitrary. Use parameters to specify which dimensions are which.
%
%   parameters:
%     timeDimension [1]: size T, all delays will be calculated as shifts along this
%       dimension
%     simultaneousDimension [2]: size C, simultaneously recorded signals,
%       as in multiple channels to be considered simultaneous when
%       computing the cross-correlation vs time
%     alignToEachOtherDimensions [3]: size G, selecting all elements (:)
%       along these dimensions provides the set of T x C matrices that will
%       be aligned in time to each other 
%    any remaining dimensions will simply be looped over, and the delays
%      reported will completely independent of other sets of values.
%
%     alignTo: 'mean' or 'first'
%
%     maxLag: maximum shift to detect, default T/2
%
%   returns:
%     delays : will have same size as data, except along timeDimension and 
%       simultaneousDimension, where it has size 1. delays provides the time
%       delay in indices used to slide each T x C matrix
% 
%   also parameters for removeDelaysTensor:
%     fillMode: how to fill the "exposed" edges of the signal matrices,
%       either a scalar like NaN or 0 or 'hold' to extend the edge values
%
%   based off Matlab's finddelay

    p = inputParser();
    p.addOptional('maxLag', [], @isscalar);
    p.addParameter('timeDimension', 1, @isscalar);
    p.addParameter('simultaneousDimensions', 2, @isvector);
    p.addParameter('alignToEachOtherDimensions', 3, @isvector);
    
    p.addParameter('upsample', 1, @isscalar);
    p.addParameter('message', 'Aligning segments in time', @ischar);
    p.addParameter('quiet', false, @islogical);
    
    p.addParameter('alignTo', 'mean', @ischar);
    p.addParameter('periodicMode', false, @islogical);
    
    % for nargout == 2, passed to removeDelaysTensor
    p.addParameter('fillMode', NaN, @(x) isscalar(x) || ischar(x));
    p.parse(varargin{:});

    alignTo = p.Results.alignTo;
    maxLag = p.Results.maxLag;
    up = p.Results.upsample;
    quiet = p.Results.quiet;
    
    hold = false;
    if ischar(p.Results.fillMode)
        if strcmp(p.Results.fillMode, 'hold')
            hold = true;
            fill = NaN;
        else
            error('Unknown fillMode');
        end
    else
        fill = p.Results.fillMode;
    end
    
    % permute and reshape so that time is along dim 1, simultaneous channels are along
    % dim 2, each group to be aligned to each other span dim 3
    % we then loop over dim 4
    szDataOrig = size(data);
    timeDim = p.Results.timeDimension;
    channelDims = ERAASR.TensorUtils.makecol(p.Results.simultaneousDimensions);
    groupDims = ERAASR.TensorUtils.makecol(p.Results.alignToEachOtherDimensions);
    loopOverDims = setdiff((1:ndims(data))', [timeDim; channelDims; groupDims]);
    szDelays = szDataOrig;
    szDelays(timeDim) = 1;
    szDelays(channelDims) = 1;
    data = ERAASR.TensorUtils.reshapeByConcatenatingDims(data, {timeDim, channelDims, groupDims, loopOverDims});
    
    T = size(data, 1);
    C = size(data, 2);
    G = size(data, 3);
    L = size(data, 4);
   
    if isempty(maxLag)
        maxLag = floor(T/2) * up;
    end
    
    if nargout > 1
        aligned = nan(T, C, G, L);
    end
    
    switch alignTo
        case 'mean'
            ref = nanmean(data, 3);
        case 'first'
            ref = data(:, :, 1, :);
        otherwise
            error('Unknown alignTo mode %s', alignTo);
    end
   
    % The largest maximum window size determines the size of the 
    % cross-correlation vector/matrix c.
    % Preallocate normalized cross-correlation vector/matrix c.
    threshMaxC = 1e-8;
    c_normalized = zeros(2*maxLag+1,C,G);
    lags = -maxLag:maxLag;
    index_max = zeros(1,1,G);
    
    delays = nan(1, 1, G, L);
    maxC = nan(1, 1, G, L);
    spuriousFound = 0;
    
    if ~quiet, prog = ERAASR.Utils.ProgressBar(L*G, p.Results.message); end
    for iL = 1:L
        x = data(:, :, :, iL); % T x C x G
        y = ref(:, :, :, iL); % T x C
        
        if up > 1
            x = upsample(x, up);
            y = upsample(y, up);
        end

        % Compute absolute values of normalized cross-correlations between x and
        % all columns of y: function XCORR does not take into account special case
        % when either x or y is all zeros, so we don't use its normalization option
        % 'coeff'. Values of normalized cross-correlations computed for a lag of
        % zero are stored in the middle row of c at index i = max_maxlag+1 (c has
        % an odd number of rows).

        for g = 1:G
            if ~quiet, prog.update(G*(iL-1) + g); end
        
            for c = 1:C
                c_normalized(:,c,g) = absnanxcorrNorm(x(:,c,g), y(:,c), maxLag);
            end
        end

        % sum the cross correlations across channels
        c_normalized_sum = ERAASR.TensorUtils.squeezeDims(sum(c_normalized, 2), 2); % nLags x G;

        if p.Results.periodicMode
            % Find indices of lags resulting in the largest absolute values of
            % normalized cross-correlations: to deal with periodic signals, seek the
            % lowest (in absolute value) lag giving the largest absolute value of
            % normalized cross-correlation.
            % Find lowest positive or zero indices of lags (negative delays) giving the
            % largest absolute values of normalized cross-correlations. 
            [max_c_pos,index_max_pos] = max(c_normalized_sum(maxLag+1:end,:),[],1);    
            % Find lowest negative indices of lags (positive delays) giving the largest
            % absolute values of normalized cross-correlations. 
            [max_c_neg,index_max_neg] = max(flipud(c_normalized_sum(1:maxLag,:)),[],1);

            max_c = nan(G, 1);

            if isempty(max_c_neg)
                % Case where MAXLAG is all zeros.
                index_max = maxLag + index_max_pos;
            else
                for g=1:G
                    if max_c_pos(g)>max_c_neg(g)
                        % The estimated lag is positive or zero.
                        index_max(g) = maxLag + index_max_pos(g);
                        max_c(g) = max_c_pos(g);
                    elseif max_c_pos(g)<max_c_neg(g)
                        % The estimated lag is negative.
                        index_max(g) = maxLag + 1 - index_max_neg(g);
                        max_c(g) = max_c_neg(g);
                    elseif max_c_pos(g)==max_c_neg(g)
                        if index_max_pos(g)<=index_max_neg(g)
                            % The estimated lag is positive or zero.
                            index_max(g) = max_maxlag + index_max_pos(g);
                            max_c(g) = max_c_pos(g);
                        else
                            % The estimated lag is negative.
                            index_max(g) = max_maxlag + 1 - index_max_neg(g);
                            max_c(g) = max_c_neg(g);
                        end 
                    end   
                end
            end
        else
            [max_c, index_max] = max(c_normalized_sum, [], 1);
        end
        
        delays(1, 1, :, iL) = lags(index_max) / up;
        maxC(1, 1, :, iL) = max_c;
        
        if nargout > 1
            for g = 1:G
                % apply alignment
                if abs(delays(1, 1, g, iL)) > threshMaxC   
                    alignedSlice = shift(x(:, :, g), round(delays(1, 1, g, iL)*up), fill, hold);
                else
                    % not significant
                    delays(1, 1, g, iL) = 0;
                    spuriousFound = spuriousFound + 1;
                    alignedSlice = x(:, :, g);
                end
                % downsample
                if up > 1
                    aligned(:, :, g, iL) = downsample(alignedSlice, up);
                else
                    aligned(:, :, g, iL) = alignedSlice;
                end
            end
        end

       
       
    end
    if ~quiet, prog.finish(); end

    % Set to zeros estimated delays for which the normalized cross-correlation
    % values are below a given threshold (spurious peaks due to FFT roundoff
    % errors).
    if spuriousFound > 0
        warning('Non-significant cross-correlation found in %d alignments', spuriousFound);
    end

    if nargout > 1
        aligned = ERAASR.TensorUtils.undoReshapeByConcatenatingDims(aligned, {timeDim, channelDims, groupDims, loopOverDims}, szDataOrig);
    end
    
    % then reshape delays to match that of data
    delays = ERAASR.TensorUtils.undoReshapeByConcatenatingDims(delays, {timeDim, channelDims, groupDims, loopOverDims}, szDelays);
end

function c = absnanxcorrNorm(x, y, maxlag)
    cxx = nansum(abs(x).^2);
    cyy = nansum(abs(y).^2);

    x = (x - nanmean(x)) / nanstd(x);
    x(isnan(x)) = 0;
    
    c = xcorr(x, y, maxlag);
    
    if cxx==0 || cyy==0
        c(:) = 0;
    else
        c = c / sqrt(cxx*cyy);
    end
end

function matUp = upsample(mat, up)
    szMat = size(mat);
    T = szMat(1);
    P = prod(szMat) / T;
    mat = reshape(mat, [T P]);
    timeMat = 0:T-1;
    timeUpsample = 0: 1/up : T-1/up;
    szUp = [numel(timeUpsample) P];
    matUp = nan(szUp);
    for p = 1:P
        matUp(:, p) = interp1(timeMat, mat(:, p), timeUpsample, 'linear', 'extrap');
    end
    szUpFull = szMat;
    szUpFull(1) = numel(timeUpsample);
    matUp = reshape(matUp, szUpFull);
end

function mat = downsample(matUp, down)
    % mat is T x ... 
    szUp = size(matUp);
    Tup = szUp(1);
    P = prod(szUp) / Tup;
    matUp = reshape(matUp, [Tup P]);
    timeUp = 0: 1/down : (Tup-1)/down;
    timeMat = 0:Tup/down-1;
    szMat = [numel(timeMat) P];
    mat = nan(szMat);
    for p = 1:P
        mat(:, p) = interp1(timeUp, matUp(:, p), timeMat, 'linear', 'extrap');
    end
    szMatFull = szUp;
    szMatFull(1) = numel(timeMat);
    mat = reshape(mat, szMatFull);
end

function alignedMat = shift(mat, delay, fill, hold)
    % mat is T x C
    T = size(mat, 1);
    idxGrab = 1:T;
    idxInsert = (1:T) - delay;
    mask = idxInsert > 1 & idxInsert < T;
    idxInsert = idxInsert(mask);
    idxGrab = idxGrab(mask);

    alignedMat = repmat(fill, size(mat));
    alignedMat(idxInsert, :) = mat(idxGrab, :);
    
    if hold
        alignedMat(1:idxInsert(1)-1, :) = repmat(alignedMat(idxInsert(1), :), [numel(1:idxInsert(1)-1) 1]);
        alignedMat(idxInsert(end)+1:T, :) = repmat(alignedMat(idxInsert(end), :), [numel(idxInsert(end)+1:T) 1]);
    end
end
