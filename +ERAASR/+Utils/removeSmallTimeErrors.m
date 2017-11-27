function time = removeSmallTimeErrors(time, timeDelta, timeReference, tol)
% remove small discrepancies from a time vector and close by time samples
    if nargin < 3
        timeReference = 0;
    end

    if nargin < 4
        tol = timeDelta / 1000;
    end
    tolScaled = tol / timeDelta;

    if iscell(time)
        % sometimes different delta for each column
        timeDelta = ERAASR.TensorUtils.singletonExpandToSize(timeDelta, size(time));
            
        for iT = 1:numel(time)
            time{iT} = innerFn(time{iT}, timeDelta(iT));
        end
    else
        time = innerFn(time, timeDelta);
    end
    
    function time = innerFn(time, timeDelta)
        if isnan(timeDelta) || timeDelta == 0
            return;
        end
        scaled = (time - timeReference) ./ timeDelta;

        fl = floor(scaled);
        mask = scaled - fl < tolScaled;
        scaled(mask) = fl(mask);

        cl = ceil(scaled);
        mask = cl - scaled < tolScaled;
        scaled(mask) = cl(mask);

        time = scaled .* timeDelta;
    end

end