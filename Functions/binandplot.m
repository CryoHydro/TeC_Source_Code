function [binEdges,y25,yMed,y75,binMeds]= binandplot(x,y,nBins)
% function [binEdges,y25,yMed,y75]= binandplot(x,y,nBins,xMin,xMax)
% 
% binEdges = xMin:(xMax-xMin)/nBins:xMax;
% [yMed,y25,y75] = deal(nan(1,nBins));
% 
% for iBin = 1:nBins
%     yMed(iBin) = nanmedian(y(x > binEdges(iBin) & x <= binEdges(iBin+1)));
%     y25(iBin) = prctile(y(x > binEdges(iBin) & x <= binEdges(iBin+1)),25);
%     y75(iBin) = prctile(y(x > binEdges(iBin) & x <= binEdges(iBin+1)),75);
%     nSamps = sum(x > binEdges(iBin) & x <= binEdges(iBin+1));
%     if nSamps > 30
%         fill([binEdges(iBin) binEdges(iBin+1) binEdges(iBin+1) ...
%             binEdges(iBin)],[y25(iBin) y25(iBin) y75(iBin) y75(iBin)],...
%             'k','FaceAlpha',0.1); hold on
%         plot([binEdges(iBin) binEdges(iBin+1)],[yMed(iBin) yMed(iBin)],...
%             '-k','LineWidth',1)
%     end
% end
% 
% end
%% alternative uses quantiles to split up
binEdges=[min(x),quantile(x,nBins)];
% binEdges = xMin:(xMax-xMin)/nBins:xMax;
[yMed,y25,y75,binMeds] = deal(nan(1,nBins));

for iBin = 1:nBins
    yMed(iBin) = nanmedian(y(x > binEdges(iBin) & x <= binEdges(iBin+1)));
    y25(iBin) = prctile(y(x > binEdges(iBin) & x <= binEdges(iBin+1)),25);
    y75(iBin) = prctile(y(x > binEdges(iBin) & x <= binEdges(iBin+1)),75);
    nSamps = sum(x > binEdges(iBin) & x <= binEdges(iBin+1));
%     if nSamps > 10
        fill([binEdges(iBin) binEdges(iBin+1) binEdges(iBin+1) ...
            binEdges(iBin)],[y25(iBin) y25(iBin) y75(iBin) y75(iBin)],...
            'k','FaceAlpha',0.1); hold on
        plot([binEdges(iBin) binEdges(iBin+1)],[yMed(iBin) yMed(iBin)],...
            '-k','LineWidth',1)
%     end
    binMeds(iBin)=nanmedian(x(x > binEdges(iBin) & x <= binEdges(iBin+1)));
end

end