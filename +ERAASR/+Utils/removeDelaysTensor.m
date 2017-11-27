function data = removeDelaysTensor(data, delays, varargin)
% delays is as given by delays = findDelaysTensor(data, varargin)
%
% parameters:
%   same as findDelaysTensor
%   fillMode: can be a scalar value like NaN or 0 to fill the edges of the
%     signals as they shift left or right. or can be 'hold' to hold the
%     edge values of the signals out to the end

    p = inputParser();
    p.addParameter('fillMode', NaN, @(x) isscalar(x) || ischar(x));
    p.addParameter('timeDimension', 1, @isscalar);
    p.addParameter('simultaneousDimensions', 2, @isvector);
    p.addParameter('alignToEachOtherDimensions', 3, @isvector);
    p.addParameter('upsample', 1, @isscalar);
    p.addParameter('message', 'Aligning segments in time', @ischar);
    p.addParameter('quiet', false, @islogical);
    p.parse(varargin{:});
    
    up = p.Results.upsample;
    quiet = p.Results.quiet;

    % permute and reshape so that time is along dim 1, simultaneous channels are along
    % dim 2, each group to be aligned to each other span dim 3
    % we then loop over dim 4
    szDataOrig = size(data);
    timeDim = p.Results.timeDimension;
    channelDims = ERAASR.TensorUtils.makecol(p.Results.simultaneousDimensions);
    groupDims = ERAASR.TensorUtils.makecol(p.Results.alignToEachOtherDimensions); % e.g. over trials or over discrete elements that are to be aligned
    loopOverDims = setdiff((1:ndims(data))', [timeDim; channelDims; groupDims]); % then we just repeat completely independently for each along this axis

    data = ERAASR.TensorUtils.reshapeByConcatenatingDims(data, {timeDim, channelDims, groupDims, loopOverDims});
    delays = ERAASR.TensorUtils.reshapeByConcatenatingDims(delays, {timeDim, channelDims, groupDims, loopOverDims});

    T = size(data, 1); % time dim
    %C = size(data, 2); % simultaneously recorded channels
    G = size(data, 3); 
    L = size(data, 4);
    
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

    if ~quiet, prog = ERAASR.Utils.ProgressBar(max(L,G), p.Results.message); end
    for iL = 1:L
        if L >= G && ~quiet, prog.update(iL); end
        for g = 1:G
            if G > L && ~quiet, prog.update(g); end
                
            slice = data(:, :, g, iL);
            
            % upsample
            if up > 1
                slice = upsample(slice, up);
            end
            % slice is Tup x C 
            
            % apply alignment
            alignedSlice = shift(slice, delays(:, :, g, iL)*up, fill, hold);
            
            % downsample
            if up > 1
                alignedSlice = downsample(alignedSlice, up);
            end
            % alignedSlice is T x C

            data(:, :, g, iL) = alignedSlice;
        end
    end
    if ~quiet, prog.finish(); end
    
    data = ERAASR.TensorUtils.undoReshapeByConcatenatingDims(data, {timeDim, channelDims, groupDims, loopOverDims}, szDataOrig);
end

function matUp = upsample(mat, up)
    % mat is T x ...
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
    idxInsert = (1:T) - round(delay);
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