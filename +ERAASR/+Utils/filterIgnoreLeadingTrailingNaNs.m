function out = filterIgnoreLeadingTrailingNaNs(B, A, mat, varargin)
    % matFilt = filterIgnoreLeadingTrailingNaNs(B, A, mat)
    
    p = inputParser;
    p.addParameter('dim', [], @(x) isempty(x) || isscalar(x));
    % if true subtracts first sample from each signal which can help reduce onset transients for IIR filters
    p.addParameter('subtractFirstSample', false, @islogical);
    % if true, the first sample will be added back on afterwards
    p.addParameter('addBackFirstSample', false, @islogical);
    p.addParameter('filtfilt', false, @islogical);
    p.addParameter('showProgress', false, @islogical);
    p.parse(varargin{:});
    
    if iscell(mat)
        if p.Results.showProgress
            prog = ERAASR.Utils.ProgressBar(numel(mat), 'Applying filter to data');
            out = cell(size(mat));
            for i = 1:numel(mat)
                prog.update(i);
                out{i} = filterMat(mat{i});
            end   
            prog.finish();
        else
            out = cellfun(@filterMat, mat, 'UniformOutput', false);
        end
    else
        out = filterMat(mat);
    end

    function mat = filterMat(mat)
        dim = p.Results.dim;
        if isempty(dim)
            % auto choose first non singular dimension
            dim = ERAASR.TensorUtils.firstNonSingletonDim(mat);
        end
        
        if dim ~= 1
            % place dim at second dimension slot
            ndimsOrig = ndims(mat);
            mat = ERAASR.TensorUtils.shiftdimToFirstDim(mat, dim);
        end
        
        nCol = size(mat(:, :), 2);
       
        % for each column, filter from the first to last non-NaN ind
        emptyCol = nan(size(mat, 1), 1);
        for iC = 1:nCol
            thisCol = mat(:, iC);
            filtCol = emptyCol;
            colMask = ~isnan(thisCol);
            start = find(colMask, 1, 'first');
            stop = find(colMask, 1, 'last');
            if p.Results.subtractFirstSample && ~isempty(start)
                firstSample = thisCol(start);
                thisCol = thisCol - firstSample;
            end
            if p.Results.filtfilt
                filtCol(start:stop) = filtfilt(B, A, double(thisCol(start:stop)));
            else
                filtCol(start:stop) = filter(B, A, double(thisCol(start:stop)));
            end
            if p.Results.subtractFirstSample && p.Results.addBackFirstSample && ~isempty(start)
                filtCol = filtCol + firstSample;
            end
            mat(:, iC) = filtCol;
        end

        if dim ~= 1
            % place dim at second dimension slot
            mat = ERAASR.TensorUtils.unshiftdimToFirstDim(mat, dim, ndimsOrig);
        end
    end
end