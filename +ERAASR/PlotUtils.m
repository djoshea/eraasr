classdef PlotUtils
    methods(Static)
        function vec = lindelta(from, delta, N)
            % Generates vector with specified size and spacing
            %   vec = ERAASR.PlotUtils.lindelta(from, delta, N) generates a column vector of N
            %   linearly spaced points from from to delta.
            vec = from:delta:(from+(N-1)*delta);
        end
        
        function plotStackedPCTraces(pcMat)
            tvec = ERAASR.PlotUtils.lindelta(0, 1/30, size(pcMat, 1));
            nPCs = size(pcMat, 2);
            labels = arrayfun(@(n) sprintf('PC %d', n), 1:nPCs, 'UniformOutput', false);
            labels(2:end-1) = {''};
            ERAASR.Utils.plotStackedTraces(tvec, pcMat, ...
                'labels', labels, 'timeUnits', 'ms', 'dataUnits', 'uV', 'timeScaleBar', true, ...
                'showVerticalScaleBars', true, 'clickable', false);
        end

        function plotOverChannels(segmentTensor, trialIdx, pulseIdx)
            pulseLen = size(segmentTensor, 2);

            tvec = ERAASR.PlotUtils.lindelta(0, 1/30, pulseLen);
            ERAASR.PlotUtils.pt(2, tvec, segmentTensor(trialIdx, :, pulseIdx, :));
            axis tight;
            axis off;
        end

        function plotOverTrials(segmentTensor, pulseIdx, chIdx)
            pulseLen = size(segmentTensor, 2);
            tvec = ERAASR.PlotUtils.lindelta(0, 1/30, pulseLen);
            ERAASR.PlotUtils.pt(2, tvec, segmentTensor(:, :, pulseIdx, chIdx), 'alpha', 1);
            axis tight;
            axis off;
        end

        function plotBeforeAfterOverTrials(before, after, pulseIdx, chIdx)
            % before after figure
            pulseLen = size(before, 2);
            tvec = ERAASR.PlotUtils.lindelta(0, 1/30, pulseLen);
            ERAASR.PlotUtils.pt(2, tvec, before(:, :, pulseIdx, chIdx), 'alpha', 1);
            hold on
            pulseLen = size(after, 2);
            tvec = ERAASR.PlotUtils.lindelta(0, 1/30, pulseLen);
            ERAASR.PlotUtils.pt(2, tvec, after(:, :, pulseIdx, chIdx), 'alpha', 1, 'Color', 'b');
            axis tight;
            axis off;
        end

        function plotBeforeAfterOverTrialsAllChannels(before, after, pulseIdx)
            nChannels = size(before, 4);
            nCol = 4;
            nRow = ceil(nChannels) / nCol;

            for c = 1:nChannels
                ERAASR.PlotUtils.subtightplot(nRow, nCol, c);
                ERAASR.PlotUtils.plotBeforeAfterOverTrials(before, after, pulseIdx, c);
            end
        end

        function plotOverTrialsAllChannels(segmentTensor, pulseIdx)
            nChannels = size(segmentTensor, 4);
            nCol = 4;
            nRow = ceil(nChannels) / nCol;

            for c = 1:nChannels
                ERAASR.PlotUtils.subtightplot(nRow, nCol, c);
                ERAASR.PlotUtils.plotOverTrials(segmentTensor, pulseIdx, c);
            end
        end

        %%%%%%%%%% over train

        function plotOverTrain(segmentTensor, trialIdx, chIdx)
            alpha = 0.8;
            pulseLen = size(segmentTensor, 2);
            tvec = ERAASR.PlotUtils.lindelta(0, 1/30, pulseLen);
            ERAASR.PlotUtils.pt(2, tvec, segmentTensor(trialIdx, :, :, chIdx), 'alpha', alpha);
            axis tight;
            axis off;
        end

        function plotOverTrainAllChannels(segmentTensor, trialIdx)
            nChannels = size(segmentTensor, 4);
            nCol = 4;
            nRow = ceil(nChannels) / nCol;

            handles = gobjects(nChannels, 1);
            for c = 1:nChannels
                ERAASR.PlotUtils.subtightplot(nRow, nCol, c);
                ERAASR.PlotUtils.plotOverTrain(segmentTensor, trialIdx, c);
                ax = gca;
                ax.Units = 'centimeters';
                ax.LooseInset = [0.1 0.1 0.1 0.1];
                if c <= nCol
                    ax.LooseInset(4) = 0.3;
                end
                axis tight;
                axis off;
                h = title(sprintf('ch %d', c));
                h.BackgroundColor = 'none';
                
                handles(c) = ax;
            end
            
            ERAASR.PlotUtils.equalizeAxisLimits(handles);
        end

        %%%%%%%%%% before / after over train

        function plotBeforeAfterOverTrain(before, after, trialIdx, chIdx)
            % before after figure
            pulseLen = size(before, 2);
            tvec = ERAASR.PlotUtils.lindelta(0, 1/30, pulseLen);
            ERAASR.PlotUtils.pt(2, tvec, before(trialIdx, :, :, chIdx), 'alpha', 1, 'Color', 'k');
            hold on

            pulseLen = size(after, 2);
            tvec = ERAASR.PlotUtils.lindelta(0, 1/30, pulseLen);
            ERAASR.PlotUtils.pt(2, tvec, after(trialIdx, :, :, chIdx), 'alpha', 1, 'Color', 'b');
            axis tight
            axis off;
        end

        function plotBeforeAfterOverTrainAllChannels(before, after, trialIdx)
            nChannels = size(before, 4);
            nCol = 4;
            nRow = ceil(nChannels) / nCol;

            handles = gobjects(nChannels, 1);
            for c = 1:nChannels
                ERAASR.PlotUtils.subtightplot(nRow, nCol, c);
                ERAASR.PlotUtils.plotBeforeAfterOverTrain(before, after, trialIdx, c);

                ax = gca;
                ax.Units = 'centimeters';
                ax.LooseInset = [0.1 0.1 0.1 0.1];
                if c <= nCol
                    ax.LooseInset(4) = 0.3;
                end
                axis tight;
                axis off;
                h = title(sprintf('ch %d', c));
                h.BackgroundColor = 'none';
                
                handles(c) = ax;
            end

            ERAASR.PlotUtils.equalizeAxisLimits(handles);
        end

        %%%%%%%%%% raw with artifact over trials

        function plotRawWithArtifactOverTrials(before, after, pulseIdx, chIdx)
            % before after figure
            nTraces = size(before, 1);
            trialIdx = 1:10:nTraces;

            pulseLen = size(before, 2);
            tvec = ERAASR.PlotUtils.lindelta(0, 1/30, pulseLen);
            ERAASR.PlotUtils.pt(2, tvec, before(trialIdx, :, pulseIdx, chIdx), 'alpha', 1, 'Color', 'k');
            hold on
            pulseLen = size(after, 2);
            tvec = ERAASR.PlotUtils.lindelta(0, 1/30, pulseLen);
            ERAASR.PlotUtils.pt(2, tvec, after(trialIdx, :, pulseIdx, chIdx), 'alpha', 1, 'Color', 'r');
            axis tight;
            axis off;
        end

        function plotRawWithArtifactOverTrialsAllChannels(before, after, pulseIdx)
            nChannels = size(before, 4);
            nCol = 4;
            nRow = ceil(nChannels) / nCol;

            handles = gobjects(nChannels, 1);
            for c = 1:nChannels
                ERAASR.PlotUtils.subtightplot(nRow, nCol, c);
                ERAASR.PlotUtils.plotRawWithArtifactOverTrials(before, after, pulseIdx, c);

                ax = gca;
                ax.Units = 'centimeters';
                ax.LooseInset = [0.1 0.1 0.1 0.1];
                if c <= nCol
                    ax.LooseInset(4) = 0.3;
                end
                handles(c) = ax;
               
                axis tight;
                axis off;
                h = title(sprintf('ch %d', c));
                h.BackgroundColor = 'none';
                
            end

            ERAASR.PlotUtils.equalizeAxisLimits(handles);
        end
        
        function plotRawWithArtifactOverPulses(before, artifact, trialIdx, chIdx)
            % before after figure
            nPulses = size(before, 3);
            pulseIdx = 1:nPulses;

            pulseLen = size(before, 2);
            tvec = ERAASR.PlotUtils.lindelta(0, 1/30, pulseLen);
            ERAASR.PlotUtils.pt(2, tvec, before(trialIdx, :, pulseIdx, chIdx), 'alpha', 1, 'Color', 'k');
            hold on
            pulseLen = size(artifact, 2);
            tvec = ERAASR.PlotUtils.lindelta(0, 1/30, pulseLen);
            ERAASR.PlotUtils.pt(2, tvec, artifact(trialIdx, :, pulseIdx, chIdx), 'alpha', 1, 'Color', 'r');

            axis tight;
            axis off;
        end
        
        function plotRawWithArtifactOverPulsesAllChannels(before, artifact, trialIdx)
            nChannels = size(before, 4);
            nCol = 4;
            nRow = ceil(nChannels) / nCol;

            handles = gobjects(nChannels, 1);
            for c = 1:nChannels
                ERAASR.PlotUtils.subtightplot(nRow, nCol, c);
                ERAASR.PlotUtils.plotRawWithArtifactOverPulses(before, artifact, trialIdx, c);

                ax = gca;
                ax.Units = 'centimeters';
                ax.LooseInset = [0.1 0.1 0.1 0.1];
                if c <= nCol
                    ax.LooseInset(4) = 0.3;
                end
                
                axis tight;
                axis off;
                h = title(sprintf('ch %d', c));
                h.BackgroundColor = 'none';
                
                handles(c) = ax;
            end

            ERAASR.PlotUtils.equalizeAxisLimits(handles);
        end
        
        function [h, cmap] = pt(timeDim, varargin)
            % pt(timeDim, dataTensor, ...)
            % pt(timeDim, timeVec, dataTensor, ...)
            %
            % parameters:
            %   colormap
            %   other parameters will be passed thru to plot(...)
            %
            % like plot except treats timeDim as the timeDimension and moves everything
            % else to the second dim to be plotted on top of it

            narg = numel(varargin);
            if isvector(varargin{1}) && narg > 1 && isnumeric(varargin{2})
                x = varargin{2};
                tvec = ERAASR.TensorUtils.makecol(varargin{1});
                args = varargin(3:end);
            else
                x = varargin{1};
                tvec = (1:size(x, timeDim))';
                args = varargin(2:end);
            end

            % other dims taken care of automatically
            otherDims = ERAASR.TensorUtils.otherDims(size(x), timeDim);
            xr = ERAASR.TensorUtils.reshapeByConcatenatingDims(x, {timeDim, otherDims});
            nTraces = size(xr, 2);

            p = inputParser();
            p.addParameter('colormap', [], @(x) isempty(x) || (~ischar(x) && ismatrix(x)));
            p.addParameter('coloreval', [], @(x) isempty(x) || isvector(x));
            p.addParameter('alpha', 0.8, @isscalar);
            p.KeepUnmatched = true;
            p.PartialMatching = false;
            p.parse(args{:});

            cmap = p.Results.colormap;
            if isempty(cmap)
                cmap = parula(nTraces);
            end

            if isempty(p.Results.coloreval)
                set(gca, 'ColorOrder', cmap, 'ColorOrderIndex', 1);
                hold on;
                h = plot(tvec, xr, p.Unmatched);
                for iH = 1:numel(h)
                    h(iH).Color(4) = p.Results.alpha;
                end
                hold off;
            else
                % plot lines according to their value in cmap
                coloreval = p.Results.coloreval;
                colorevalLims = [nanmin(coloreval(:)), nanmax(coloreval(:))];
                coloreval = ERAASR.TensorUtils.rescaleIntervalToInterval(coloreval, colorevalLims, [0 1]);
                colors = ERAASR.PlotUtils.evalColorMapAt(cmap, coloreval);

                hold on;
                %     h = plot(tvec, xr, p.Unmatched);
                h = stairs(tvec, xr);

                for iH = 1:numel(h)
                    if any(isnan(colors(iH, :)))
                        delete(h(iH));
                    else
                        h(iH).Color = cat(2, colors(iH, :), p.Results.alpha);
                    end
                end
                hold off;

                ax = gca;
                ax.TickDir = 'out';
                ax.ColorSpace.Colormap = cmap;
                ax.CLim = colorevalLims;
                hc = colorbar;
                hc.TickDirection = 'out';
                
                axis off
            end
        end
        
        function ceval = evalColorMapAt(cmap, at)
            if isa(cmap, 'function_handle')
                cmap = cmap(1000);
            end
            N = size(cmap, 1);
            ceval = interp1((0:N-1) / (N-1), cmap, at);
        end
        
        function equalizeAxisLimits(hax, which)
            if nargin < 2
                which = 'xy';
            end
            doX = ismember('x', which);
            doY = ismember('y', which);
            doZ = ismember('z', which);

            argsX = {};
            argsY = {};
            argsZ = {};

            if doX
                xl = cell2mat(get(hax, 'XLim'));
                xl = [min(xl(:, 1)) max(xl(:, 2))];
                argsX = {'XLim', xl};
            end

            if doY
                yl = cell2mat(get(hax, 'YLim'));
                yl = [min(yl(:, 1)) max(yl(:, 2))];
                argsY = {'YLim', yl};
            end

            if doZ
                zl = cell2mat(get(hax, 'ZLim'));
                zl = [min(zl(:, 1)) max(zl(:, 2))];
                argsZ = {'ZLim', zl};
            end

            set(hax, argsX{:}, argsY{:}, argsZ{:});
        end
        
        function h=subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
            %function h=subtightplot(m,n,p,gap,marg_h,marg_w,varargin)
            %
            % Functional purpose: A wrapper function for Matlab function subplot. Adds the ability to define the gap between
            % neighbouring subplots. Unfotrtunately Matlab subplot function lacks this functionality, and the gap between
            % subplots can reach 40% of figure area, which is pretty lavish.  
            %
            % Input arguments (defaults exist):
            %   gap- two elements vector [vertical,horizontal] defining the gap between neighbouring axes. Default value
            %            is 0.01. Note this vale will cause titles legends and labels to collide with the subplots, while presenting
            %            relatively large axis. 
            %   marg_h  margins in height in normalized units (0...1)
            %            or [lower uppper] for different lower and upper margins 
            %   marg_w  margins in width in normalized units (0...1)
            %            or [left right] for different left and right margins 
            %
            % Output arguments: same as subplot- none, or axes handle according to function call.
            %
            % Issues & Comments: Note that if additional elements are used in order to be passed to subplot, gap parameter must
            %       be defined. For default gap value use empty element- [].      
            %
            % Usage example: h=subtightplot((2,3,1:2,[0.5,0.2])

            if (nargin<4) || isempty(gap),    gap=0.01;  end
            if (nargin<5) || isempty(marg_h),  marg_h=0.05;  end
            if (nargin<5) || isempty(marg_w),  marg_w=marg_h;  end
            if isscalar(gap),   gap(2)=gap;  end
            if isscalar(marg_h),  marg_h(2)=marg_h;  end
            if isscalar(marg_w),  marg_w(2)=marg_w;  end
            gap_vert   = gap(1);
            gap_horz   = gap(2);
            marg_lower = marg_h(1);
            marg_upper = marg_h(2);
            marg_left  = marg_w(1);
            marg_right = marg_w(2);

            %note n and m are switched as Matlab indexing is column-wise, while subplot indexing is row-wise :(
            [subplot_col,subplot_row]=ind2sub([n,m],p);  

            % note subplot suppors vector p inputs- so a merged subplot of higher dimentions will be created
            subplot_cols=1+max(subplot_col)-min(subplot_col); % number of column elements in merged subplot 
            subplot_rows=1+max(subplot_row)-min(subplot_row); % number of row elements in merged subplot   

            % single subplot dimensions:
            %height=(1-(m+1)*gap_vert)/m;
            %axh = (1-sum(marg_h)-(Nh-1)*gap(1))/Nh; 
            height=(1-(marg_lower+marg_upper)-(m-1)*gap_vert)/m;
            %width =(1-(n+1)*gap_horz)/n;
            %axw = (1-sum(marg_w)-(Nw-1)*gap(2))/Nw;
            width =(1-(marg_left+marg_right)-(n-1)*gap_horz)/n;

            % merged subplot dimensions:
            merged_height=subplot_rows*( height+gap_vert )- gap_vert;
            merged_width= subplot_cols*( width +gap_horz )- gap_horz;

            % merged subplot position:
            merged_bottom=(m-max(subplot_row))*(height+gap_vert) +marg_lower;
            merged_left=(min(subplot_col)-1)*(width+gap_horz) +marg_left;
            pos_vec=[merged_left merged_bottom merged_width merged_height];

            % h_subplot=subplot(m,n,p,varargin{:},'Position',pos_vec);
            % Above line doesn't work as subplot tends to ignore 'position' when same mnp is utilized
            h=subplot('Position',pos_vec,varargin{:});
        end
    end
end