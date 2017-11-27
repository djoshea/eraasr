function [matClean, artPcs, artMat] = cleanMatrixViaPCARegression(mat, nPCs, varargin)

    p = inputParser();
    p.addParameter('omitAdjacentChannelsBandWidth', 1, @isscalar); % number of adjacent channels to omit from the regression
    p.addParameter('pcaOnlyOmitted', false, @islogical); % original method had this false
    p.parse(varargin{:});

    omit = p.Results.omitAdjacentChannelsBandWidth;

    % M x N (time x neurons)
    mat = bsxfun(@minus, mat, nanmean(mat, 1));

    N= size(mat, 2);
    coeff = pca(mat);

    cpart = coeff(:, 1:nPCs); % N x K (neurons --> PCs)
    artPcs = mat * cpart;
    
    matClean = nan(size(mat));
    artMat = nan(size(mat));
    
    for n = 1:N
        % build PCs but omit contributions from adjacent channels
        minOmit = n - floor((omit-1)/2);
        maxOmit = n + ceil((omit-1)/2);
        
        if ~p.Results.pcaOnlyOmitted
            % pcs from all but set coefficients to zero post hoc
            channelsOmit = intersect(1:N, minOmit:maxOmit);
            cpartOmit = cpart;
            cpartOmit(channelsOmit, :) = 0;
            pcs = mat * cpartOmit; % (MxN)*(NxK) == (MxK);

        else
            % new mode, build pcs only from omitted channels
            channelKeep = setdiff(1:N, minOmit:maxOmit);
            [~, pcs] = pca(mat(:, channelKeep));
            pcs = pcs(:, 1:min(nPCs, size(pcs, 2)));
        end

        % do regression to reconstruct this channel from pcs
        seg = mat(:, n);
        artMat(:, n) =  pcs * (pcs \ seg);
        seg = seg - artMat(:, n);
        matClean(:, n) = seg;
    end
end

% recon
% artEach = bsxfun(@times, pcs, (pcs\seg)');