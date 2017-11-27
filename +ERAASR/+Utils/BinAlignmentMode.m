classdef BinAlignmentMode < int32
    enumeration
        Centered (1)
        Causal (2) % align bin so that all spikes came before the timepoint
        Acausal (3) % align bin so that all spikes came before the timepoint
    end
   
    methods
        function offsetStart = getBinStartOffsetForBinWidth(mode, binWidth)
            switch mode
                case ERAASR.Utils.BinAlignmentMode.Centered
                    offsetStart = -binWidth/2;
                case ERAASR.Utils.BinAlignmentMode.Causal
                    offsetStart = -binWidth;
                case ERAASR.Utils.BinAlignmentMode.Acausal
                    offsetStart = 0;
            end
        end
        
        function offset = getOffsetToBinCenter(mode, binWidth)
            % add this to a timestamp to get the bin center
            switch mode
                case ERAASR.Utils.BinAlignmentMode.Centered
                    offset = 0;
                case ERAASR.Utils.BinAlignmentMode.Causal
                    offset = -binWidth/2;
                case ERAASR.Utils.BinAlignmentMode.Acausal
                    offset = binWidth/2;
            end
        end
        
        function offsetStop = getBinStopOffsetForBinWidth(mode, binWidth)
            offsetStop = mode.getBinStartOffsetForBinWidth(binWidth) + binWidth;
        end
        
        function [tMinNew, tMaxNew] = getTimeLimitsForRebinning(mode, tMinOrig, tMaxOrig, origDelta, newDelta, timeReference)
            if nargin < 6
                timeReference = 0;
            end
            
            tol = min(origDelta, newDelta) / 1000;
            
            if origDelta == 0
                assert(isequaln(tMinOrig, tMaxOrig));
                % single time sample
                tMinNew = tMinOrig;
                tMaxNew = tMaxOrig;
            else
            
                switch mode
                    % the assumption here is that we treat the existing time
                    % bins (with width origDelta) as having the same
                    % binAlignmentMode as mode.
                    case ERAASR.Utils.BinAlignmentMode.Causal
                        tMinNew = timeReference + ERAASR.Utils.ceiltol((tMinOrig - origDelta - timeReference + newDelta) ./ newDelta, tol) .* newDelta;
                        tMaxNew = timeReference + ERAASR.Utils.floortol((tMaxOrig - timeReference) / newDelta, tol) * newDelta;
                    case ERAASR.Utils.BinAlignmentMode.Acausal
                        tMinNew = timeReference + ERAASR.Utils.ceiltol((tMinOrig - timeReference) / newDelta, tol) .* newDelta;
                        tMaxNew = timeReference + ERAASR.Utils.floortol((tMaxOrig + origDelta - timeReference - newDelta) ./ newDelta, tol) .* newDelta;
                    case ERAASR.Utils.BinAlignmentMode.Centered
                        tMinNew = timeReference + ERAASR.Utils.ceiltol((tMinOrig - origDelta/2 - timeReference + newDelta/2) ./ newDelta, tol) .* newDelta;
                        tMaxNew = timeReference + ERAASR.Utils.floortol((tMaxOrig + origDelta/2 - timeReference - newDelta/2) ./ newDelta, tol) .* newDelta;
                end
            end
                    
            mask = tMinNew > tMaxNew;
            tMinNew(mask) = NaN;
            tMaxNew(mask) = NaN;
        end
        
        function [tvecLabelCell, tbinsForHistcCell] = generateMultipleBinnedTimeVectors(mode, starts, stops, binWidth, finalDelta, timeReference)
            % tvecLabel will be the time that each time bin is labeled
            % with, tbinsForHistc is the argument to histc that will select
            % times in the bin so that the output of histc matches
            % tvecLabel. The last bin output by histc should be thrown
            % away. The time limits are chosen that the first bin is
            % located entirely within starts:stops, but the last bin is
            % extended to ensure that stops is included within the last
            % bin. This is to ensure consistent behavior when starts ==
            % stops, when we would then want a single bin extending from
            % starts : (starts + binWidth).
            %
            % actual spike zero/start/stop will be at 0/start/stop + binAlignmentMode.getBinStartOffsetForBinWidth
                
            assert(numel(starts) == numel(stops));
            if nargin < 5
                finalDelta = binWidth;
            end
            if nargin < 6
                timeReference = 0;
            end
            
            [tvecLabelCell, tbinsForHistcCell] = deal(cell(numel(starts), 1));
            
            binOffset = mode.getBinStartOffsetForBinWidth(binWidth);
            
            [tMin, tMax] = mode.getTimeLimitsForRebinning(starts, stops, binWidth, finalDelta, timeReference);
            
            for i = 1:numel(starts)
                if isnan(tMin(i)) || isnan(tMax(i))
                    continue;
                end
                
                tvecLabel = (tMin(i):binWidth:tMax(i))';
                if isempty(tvecLabel)
                    continue;
                end
                
                % don't include the offset from the time labels
                tvecLabelCell{i} = tvecLabel;
                
                % add the end of the last bin onto the histc bins
                tbinsForHistcCell{i} = [tvecLabel+binOffset; tvecLabel(end)+binOffset+binWidth];
            end
        end
        
        function [tvecLabel, tbinsForHistc, tbinsValidMat] = generateCommonBinnedTimeVector(mode, starts, stops, binWidth)
            % note that tbinsForHistc will have one additional entry than
            % tbinsValidMat and tvecLabel. This is for histc. Since histc returns 
            % an additional bin which captures values equal to the final
            % value, you should throw away the final bin returned by histc.
            start = nanmin(starts);
            stop = nanmax(stops);
            [tvecLabelCell, tbinsForHistcCell] = mode.generateMultipleBinnedTimeVectors(start, stop, binWidth);
            tvecLabel = tvecLabelCell{1};
            tbinsForHistc = tbinsForHistcCell{1};
            
            binStartOffset = mode.getBinStartOffsetForBinWidth(binWidth);
            
            nTrials = numel(starts);
            tbinsValidMat = false(nTrials, numel(tvecLabel));
            for iTrial = 1:nTrials
                tbinsValidMat(iTrial, :) = tbinsForHistc(1:end-1) >= starts(iTrial)+binStartOffset & ...
                    tbinsForHistc(2:end) <= stops(iTrial) + binWidth + binStartOffset; % ensure that we take the last bin labeled at exactly stop(iTrial)
            end
        end
    end
    
end