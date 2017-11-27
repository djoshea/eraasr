function [traceCenters, hLines] = plotStackedTraces(tvec, data, varargin)
% [traceCenters] = plotStacked(tvec, data, varargin)
% Plots the columns of data stacked vertically on top of each other
%
% For data as numeric tensor, tvec is vector matching size(data, 1)
%   Dim 1 of data is time, dim 2 is traces to stack vertically [nTraces]
%   If dim 3 has length > 1, these traces will be superimposed at the same
%   location with the colors set by colorMap [nTracesSuperimposed]
%
% For data as cell matrix, tvec is cell with same size
%   Dim 1 of data is traces to be stacked vertically [nTraces]
%   Dim 2 is traces to be superimposed [nTracesSuperimposed]
%   Inside each cell time is dim 1, multiple traces with same color can be plotted over cols
%
% nTraces is 
%
% options include all options for plot(...) plus:
%   normalize: [false] scale each signal to the same range before plotting
%   spacingFraction: [1.02] space each trace by this fraction of the previous
%     trace's range

% TODO : fix intercalate - probably not working anymore

p = inputParser();
p.addParameter('evenSpacing', false, @islogical); % the vertical space allocated to each stacked trace is the same?
p.addParameter('normalize', false, @islogical); % the vertical height of each trace is normalized? or in original data units
p.addParameter('intercalate', false, @islogical); % the traces should be squished together as close as possible without touching
p.addParameter('spacingFraction', 1.2, @isscalar); % the gap between each trace / the height of those traces
p.addParameter('colormap', [], @(x) isempty(x) || isa(x, 'function_handle') || ismatrix(x)); % for superimposed traces 
p.addParameter('colormapStacked', [], @(x) isempty(x) || isa(x, 'function_handle') || ismatrix(x)); % for stacked traces, only used if colormap is empty
p.addParameter('maintainScaleSuperimposed', true, @islogical); % when superimposing multiple traces, keep the relative size and offset between the superimposed traces
p.addParameter('labels', {}, @(x) isempty(x) || isvector(x)); % labels over nTraces for the y axis
p.addParameter('labelRotation', 0, @isvector);
p.addParameter('labelsSuperimposed', {}, @iscell); % labels over the nSuperimposed traces, for clickable descriptions
p.addParameter('showLabels', 'auto', @(x) islogical(x) || ischar(x)); % show the labels on the left axis, slow if too many traces, 'auto' is true if nTraces < 25
p.addParameter('clickable', false, @islogical); % make each trace clickable and show a description
p.addParameter('timeUnits', '', @ischar); 
p.addParameter('timeScaleBar', false, @islogical); % use scale bar instead of tick bridge for time axis?
p.addParameter('dataUnits', [], @(x) ischar(x) || (isvector(x) && iscellstr(x))); % either a string describing units for all traces, or a nTraces x 1 cell of units for each set of traces running vertically
p.addParameter('verticalScaleBarHideLabel', false, @islogical);
p.addParameter('showVerticalScaleBars', false, @(x) islogical(x) || ischar(x)); % show intelligent y axis scale bars on the right hand side
p.addParameter('showDataRanges', false, @(x) islogical(x) || ischar(x)); % show intelligent y axis scale bars on the right hand side
p.addParameter('showSpanLines', true, @islogical);
p.addParameter('dataRangeFormat', '%.4g', @ischar);
% p.addParameter('lineStyle', '-', @ischar);
p.KeepUnmatched = true;
p.CaseSensitive = false;
p.parse(varargin{:});

if ~iscell(data)
    % all traces share common time vector
    nTraces = size(data, 2);
    nSuperimposed =size(data, 3);
    nTime = size(data, 1);
    if isempty(tvec)
        tvec = ERAASR.TensorUtils.makecol(0:nTime-1);
    end
    assert(numel(tvec) == nTime, 'Time vector must match size(data, 1)');
    
    if ~isfloat(data)
        data = single(data);
    end
else
    % everytraces has different time vector
    nTraces = size(data, 1);
    nSuperimposed = size(data, 2);
    emptyMask = cellfun(@isempty, data);
    data(emptyMask) = {NaN}; % prevents errors with cellfun later
end
   
% construct labels by traces
if isempty(p.Results.labels)
    labels = arrayfun(@num2str, 1:nTraces, 'UniformOutput', false);
elseif isnumeric(p.Results.labels)
    labels = arrayfun(@num2str, p.Results.labels, 'UniformOutput', false);
else
    labels = p.Results.labels;
end
labels = ERAASR.TensorUtils.makecol(labels);

% construct labels for each superimposed line within each trace
if isempty(p.Results.labelsSuperimposed)
    labelsLinesWithinEachTrace = arrayfun(@num2str, 1:nSuperimposed, 'UniformOutput', false);
else
    labelsLinesWithinEachTrace = p.Results.labelsSuperimposed;
end

if ~iscell(data)
    % invert the order of the traces so the first is plotted at the top
    if p.Results.maintainScaleSuperimposed
        % subtract the min so each group of traces has min == 0
        minEachGroup = ERAASR.TensorUtils.nanminMultiDim(data, [1 3]);
        dataLowOrig = minEachGroup;
        matShift = bsxfun(@minus, data, minEachGroup);

        rangesOrig = ERAASR.TensorUtils.nanmaxMultiDim(matShift, [1 3]);
        if p.Results.normalize
            norms = rangesOrig;
            matShift = bsxfun(@rdivide, matShift, norms);
        else
            norms = ones(nTraces, 1);
        end
    else
        % subtract the min so each trace has min at zero
        dataLowOrig = nanmin(data, [], 1);
        matShift = bsxfun(@minus, data, dataLowOrig);

        if p.Results.normalize
            maxEach = range(matShift, 1);
            matShift = bsxfun(@rdivide, matShift, maxEach);
        end
    end

    % compute the max range each row
    rangesNorm = ERAASR.TensorUtils.nanmaxMultiDim(matShift, [1 3]);
    
    rangesNorm(isnan(rangesNorm)) = 0;

    % figure out where each trace should start
    if p.Results.evenSpacing
        % trace(k) will be offset by spacing * k
        if p.Results.intercalate
            deltas = matShift(:, 2:end, :) * p.Results.spacingFraction - matShift(:, 1:end-1, :);
            maxDeltas = ERAASR.TensorUtils.nanmaxMultiDim(deltas, [1 3]); % max over time and superimposed traces
            traceOffsets = (nTraces-1:-1:0) * nanmax(maxDeltas);
        else
            rangesPadded = ERAASR.TensorUtils.makerow(rangesNorm * (p.Results.spacingFraction));
            traceOffsets = (nTraces-1:-1:0) * nanmax(rangesPadded);
        end
    else
        if p.Results.intercalate
            % the offsets should be set so that pairs of points across
            % all time points should satisfy:
            %   traceN+1(t) * spacingFraction <= traceN(t) + offset
            % or:
            %   offset = max_t (traceN+1(t) * spacingFraction - traceN(t))

            deltas = matShift(:, 2:end, :) * p.Results.spacingFraction - matShift(:, 1:end-1, :);
            maxDeltas = ERAASR.TensorUtils.nanmaxMultiDim(deltas, [1 3]); % max over time and superimposed traces
            cs = fliplr(cumsum(fliplr(maxDeltas)));
            traceOffsets = [cs, 0];

        else 
            rangesPadded = ERAASR.TensorUtils.makerow(rangesNorm * (p.Results.spacingFraction));
            cs = fliplr(cumsum(fliplr(rangesPadded)));
            traceOffsets = [cs(2:end), 0];
        end
    end  

    matShift = bsxfun(@plus, matShift, traceOffsets);

    % expand colormap to be exactly nSuperimposed long
    [map, colorByStack] = getColormap(p.Results.colormapStacked, nTraces, p.Results.colormap, nSuperimposed);

    if nSuperimposed == 1 && ~colorByStack
        % plot simultaneously
        hLines = plot(tvec, matShift, '-', 'Color', map(1, :), 'MarkerFaceColor', map(1, :), 'MarkerEdgeColor', map(1, :), p.Unmatched);
    else
%         set(gca, 'ColorOrder', map); % the map will superimpose automatically
        hold on;

        % here we arrange so that all traces are stacked vertically but that
        % all positions along dim 2 are grouped together 
        matShiftCat = ERAASR.TensorUtils.reshapeByConcatenatingDims(matShift, {1 [3 2]});

        hLines = plot(tvec, matShiftCat, '-', p.Unmatched);
        hLines = reshape(hLines, [nSuperimposed nTraces])';
        
        if colorByStack
            for iT = 1:nTraces
                set(hLines(iT, :), 'Color', map(iT, :), 'MarkerFaceColor', map(iT, :), 'MarkerEdgeColor', map(iT, :));
            end
        else
            for iS = 1:nSuperimposed
                set(hLines(:, iS), 'Color', map(iS, :), 'MarkerFaceColor', map(iS, :), 'MarkerEdgeColor', map(iS, :));
            end
        end
    end

else
    % cell mode - everyone has a different time vector
    if ~iscell(tvec)
        tvec = repmat({tvec}, nTraces, nSuperimposed);
    end
    
    [map, colorByStack] = getColormap(p.Results.colormapStacked, nTraces, p.Results.colormap, nSuperimposed);
        
    % invert the order of the traces so the first is plotted at the top
%     data = flipud(data);
%     labels = flipud(labels);
%     tvec = flipud(tvec);

    if p.Results.maintainScaleSuperimposed
        % subtract the min so each group trace has min at zero
        minEachRow = nanmin(cellfun(@(data) double(nanminNanEmpty(data(:))), data), [], 2);
        dataLowOrig = minEachRow;
        minCell = num2cell(repmat(minEachRow, 1, size(data, 2)));
        cellShift = cellfun(@(data, min) data - min, data, minCell, 'UniformOutput', false);

        rangesOrig = nanmax(cellfun(@(data) double(nanmax(data(:))), cellShift), [], 2);
        if p.Results.normalize
            maxCell = num2cell(repmat(rangesOrig, 1, size(data, 2)));
            cellShift = cellfun(@(matShift, max) double(matShift ./ max), cellShift, maxCell, 'UniformOutput', false);
            norms = rangesOrig;
        else
            norms = ones(nTraces, 1);
        end
    else
        % subtract the min so each trace has min at zero
        dataLowOrig = nanmin(data, [], 1);
        cellShift = cellfun(@(data) bsxfun(@minus, data, dataLowOrig), data, 'UniformOutput', false);

        if p.Results.normalize
            cellShift = cellfun(@(matShift) bsxfun(@rdivide, matShift, range(matShift, 1)), 'UniformOutput', false);
        end
    end

    % compute the max range each row
    rangesNorm = max(cellfun(@(matShift) double(nanmaxNanEmpty(matShift(:), [], 1)), cellShift), [], 2);

    rangesNorm(isnan(rangesNorm)) = 0;
    
    % can't intercalate without doing time consuming interpolation
    rangesPadded = rangesNorm * (p.Results.spacingFraction);
    cs = flipud(cumsum(flipud(rangesPadded)));
    traceOffsets = [cs(2:end); 0];

    traceOffsets = traceOffsets - min(traceOffsets);
%     traceCenters = (traceOffsets + rangesNorm / 2)';

    hLines = cell(nTraces, nSuperimposed);
    %lineDescriptionsCell = cell(nTraces, nSuperimposed);
    for iT = 1:nTraces
        for iS = 1:nSuperimposed
            matShift = cellShift{iT, iS} + traceOffsets(iT);
            if isempty(tvec{iT, iS})
                tvecThis = 1:numel(matShift);
            else
                tvecThis = tvec{iT, iS};
            end

            if colorByStack
                color = map(iT, :);
            else
                color = map(iS, :);
            end
            hLines{iT, iS} = plot(tvecThis, matShift, '-', 'Color', color, ...
                'MarkerFaceColor', color, 'MarkerEdgeColor', color, p.Unmatched);
            hold on;
        end
    end
end

traceCenters = (traceOffsets + rangesNorm / 2)';
traceHighs = traceOffsets + rangesNorm;
traceLows = traceOffsets;
dataHighOrig = dataLowOrig + rangesOrig;

ylim([nanmin(traceLows) nanmax(traceHighs)]);

if strcmp(p.Results.showLabels, 'auto')
    showLabels = nTraces < 25;
else
    showLabels = p.Results.showLabels;
end

if isempty(p.Results.timeUnits)
    xlab = 'Time';
else
    xlab = sprintf('Time (%s)', p.Results.timeUnits);
end

axis tight;
axis off;

end

function r = nanmaxNanEmpty(v1, varargin)
    if isempty(v1) && numel(varargin) >= 1 && isempty(varargin{1})
        r = NaN;
    else
        r = nanmax(v1, varargin{:});
    end
end

function [map, colorByStack] = getColormap(cmapStacked, nStacked, cmapSuperimposed, nSuperimposed)
    if ~isempty(cmapSuperimposed)
        cmapFn = cmapSuperimposed;
        colorByStack = false;
        N = nSuperimposed;
        map = expandWrapColormap(cmapFn, N);
    elseif ~isempty(cmapStacked)
        cmapFn = cmapStacked;
        colorByStack = true;
        N = nStacked;
        map = expandWrapColormap(cmapFn, N);
    else
        % no colormap specified
        if nSuperimposed == 1
            map = [0 0 0];
            colorByStack = false;
        else
            map =  expandWrapColormap(@jet, nSuperimposed);
            colorByStack = false;
        end
    end
end

function map = expandWrapColormap(map, n)
    if isa(map, 'function_handle')
        map = map(n);
    else
        map = convertColorsToMatrix(map);
        if size(map, 1) < n
            % make at least big enough
            map = repmat(map, ceil(n / size(map, 1)), 1);
            % cut it down to size (typically not necessary)
            map = map(1:n, :);
        end
    end
end

function map = convertColorsToMatrix(colors)

    hasAlpha = false;

    if ischar(colors)
        map = convStr(colors);
    elseif iscell(colors)
        N = numel(colors);
        map = nan(N, 4);

        for i = 1:N
            if ischar(colors{i})
                map(i, 1:3) = convStr(colors{i});
            elseif isvector(colors{i})
                n = numel(colors{i});
                if n == 3
                    map(i, 1:3) = colors{i};
                elseif n == 4
                    map(i, :) = colors{i};
                    hasAlpha = true;
                else
                    error('Numeric cell contents must be 3 or 4 vectors');
                end
            else
                error('Cell contents must be char or vectors');
            end
        end

        if ~hasAlpha
            map = map(:, 1:3);
        end

    elseif ismatrix(colors)
        assert(size(colors, 2) == 3 || size(colors, 2) == 4, 'Color matrix must have 3 or 4 columns'); 
        map = colors;
    else
        error('Colors must be cell array or matrix');
    end
end

function c = convStr(s)
    switch s
        case 'none'
            c = 'none';
        case 'k'
            c = [0 0 0];
        case 'b'
            c = [0 0 1];
        case 'g'
            c = [0 1 0];
        case 'c'
            c = [0 1 1];
        case 'r'
            c = [1 0 0];
        case 'm'
            c = [1 0 1];
        case 'y'
            c = [1 1 0];
        case 'w'
            c = [1 1 1];
        otherwise
            error('Unknown color %s', s);
    end
end

function r = nanminNanEmpty(v)
    if ~isempty(v)
        r = nanmin(v);
    else
        r = NaN; 
    end
end