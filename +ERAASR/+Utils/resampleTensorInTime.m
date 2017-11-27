function [data, timeNew] = resampleTensorInTime(data, timeDim, time, varargin)
% [mat, tvec] = embedTimeseriesInMatrix(dataCell, timeData, varargin)
% 
    p = inputParser();
    p.addParameter('interpolateMethod', 'linear', @ischar);
    p.addParameter('timeDelta', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('timeReference', 0, @isscalar);
    p.addParameter('binAlignmentMode', ERAASR.Utils.BinAlignmentMode.Centered, @(x) isa(x, 'ERAASR.Utils.BinAlignmentMode'));
    p.addParameter('resampleMethod', 'filter', @ischar); % valid modes are filter, average, repeat , interp   
    p.addParameter('origDelta', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('tMin', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('tMax', [], @(x) isempty(x) || isscalar(x));
    p.addParameter('uniformlySampled', false, @islogical); % can speed things up if you know it's arleady uniform
    p.parse(varargin{:});
    
    if isempty(data)
        timeNew = [];
        return;
    end
    
    assert(isvector(time));
    nTime = numel(time);
    assert(size(data, timeDim) == nTime);
    
    % figure out new time vector
    timeDelta = double(p.Results.timeDelta);
    timeReference = p.Results.timeReference;
    interpolateMethod = p.Results.interpolateMethod;
    binAlignmentMode = p.Results.binAlignmentMode;
   
    origDelta = p.Results.origDelta;
    if isempty(origDelta)
        origDelta = nanmedian(diff(time));
    end
    if isempty(timeDelta)
        timeDelta = origDelta;
    end
    
    time = ERAASR.Utils.removeSmallTimeErrors(time, timeDelta, timeReference);

    tMin = p.Results.tMin;
    tMax = p.Results.tMax;
    [tMinCalc, tMaxCalc] = p.Results.binAlignmentMode.getTimeLimitsForRebinning(min(time), max(time), origDelta, timeDelta, timeReference);
    if isempty(tMin), tMin = tMinCalc; end
    if isempty(tMax), tMax = tMaxCalc; end
    
    timeNew = (tMin:timeDelta:tMax)';
    nDimsOrig = ndims(data);
    
    if origDelta == 0 && tMin == tMax
        % single sample special case
        return;
    end
    
    data = ERAASR.TensorUtils.shiftdimToFirstDim(data, timeDim);
    deltaIsChanging = ~ERAASR.Utils.isequaltol(timeDelta, origDelta, origDelta / 1000);
    
    % do this to avoid off by one errors when generating time vectors
    if ~deltaIsChanging
        origDelta = timeDelta;
    end
    
    % build time vector for the original that starts at the appropriate
    % tMin so that we end up with the right samples
    timeUniform = (tMin:origDelta:tMax)';
    
    switch p.Results.resampleMethod
        case 'filter'
            if ~p.Results.uniformlySampled
                % sample to uniform grid
                data = interp1(time, data, timeUniform, interpolateMethod, 'extrap');
                time = timeUniform;
            end

            if deltaIsChanging
                % use resampling
                [data] = ERAASR.Utils.resamplePadEdges(data, time, timeNew, p.Results.binAlignmentMode, interpolateMethod);   
            end

            assert(size(data, 1) == numel(timeNew));
            
        case 'repeat'
            % TODO NEED TO TAKE BINALIGNMENTMODE INTO ACCOUNT
            
            % sample to uniform grid
            data = interp1(time, data, timeUniform, interpolateMethod);
            
            if timeDelta < origDelta
                % upsample via repelem
                P = origDelta/timeDelta;
                assert(ceil(P) == P, 'New sampling interval must be an integer multiple of old sampling interval');
                data = repelem(data, P, 1);
                
                % trim extra copies at the end
                data = data(1:numel(timeNew), :, :);
                
            elseif ~deltaIsChanging
                % fine as is
            else
                error('Cannot use repeat when downsampling');
            end
            
        case 'average'
            % TODO NEED TO TAKE BINALIGNMENTMODE INTO ACCOUNT
            
            % sample to uniform grid
            data = interp1(time, data, timeUniform, interpolateMethod);
            
            if timeDelta > origDelta
                % downsample via averaging
                P = timeDelta/origDelta;
                assert(ceil(P) == P, 'New sampling interval must be an integer multiple of old sampling interval');
                data = blockproc(data, [P 1], @(block) mean(block.data));
                
                % trim average from partial blocks at the end
                data = data(1:numel(timeNew), :, :);
                
            elseif ~deltaIsChanging
                % fine as is
            else
                error('Cannot use repeat when downsampling');
            end
                
        case 'interp'
            % first shift the original timestamps to their centers
            % if causal bins, e.g. 20 ms, want to interp data to centers at
            % +10, +30 and then label as 20, 40. offset is -10
            
            interpInputTimes =  time + binAlignmentMode.getOffsetToBinCenter(origDelta);
            interpOutputTimes = timeNew + binAlignmentMode.getOffsetToBinCenter(timeDelta);
            data = interp1(interpInputTimes, data, interpOutputTimes,  'linear', 'extrap');
            
%             interpToTime = timeNew + p.Results.binAlignmentMode.getOffsetToBinCenter(timeDelta);
%             data = interp1(time, data, interpToTime, interpolateMethod);
%             
        otherwise
            error('Unknown resampleMethod %s', p.Results.resampleMethod)
    end
    
    data = ERAASR.TensorUtils.unshiftdimToFirstDim(data, timeDim, nDimsOrig);
end