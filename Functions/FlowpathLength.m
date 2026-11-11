%% Glacier flow-path distance from headwalls/ridges
% Requires TopoToolbox 3
%

function [FPL] = FlowpathLength(DEM,glacierID)
%% Glacier headwall-to-cell flow-path distance
% Requires TopoToolbox 3
%
% Inputs:
%   DEM                     GRIDobj DEM (assumed projected in metres)
%   glacierID      GRIDobj; glacier IDs, NaN/0 outside glacier
%
% Output:
%   FPL                        GRIDobj; flow distance [m] from the most
%                              distant upstream headwall/divide.
%
% Notes:
% - A single FLOWobj is derived for speed.
% - Smoothing is performed separately within each glacier ID, avoiding
%   elevation leakage across glacier margins.
% - 'upstream' flowdistance returns 0 at headwall/source cells and values
%   increasing down-glacier.

assert(isequal(size(DEM.Z), size(glacierID.Z)), ...
    'DEM and glacier-ID GRIDobj must have identical dimensions.');


f = fspecial('gaussian',60,2);
    DEMs = DEM;
    DEMs.Z = filter2(f, filter2(f, DEMs.Z));
    clear f

    glacierID.Z(isnan(glacierID.Z)) = 0;

    DEMg = DEMs; DEMg.Z(glacierID.Z==0)=NaN;
    FDg  = FLOWobj(DEMg,'preprocess','carve');
    FPL  = flowdistance(FDg,'upstream');
    clear DEMs DEMg FDg

    FPL.Z(glacierID.Z==0) = NaN;
    uG = unique(glacierID.Z(:)); uG(uG == 0) = [];
    % iterate an inverse 'upstream' flowpath for each glacier
    % individually using it's maximum achieved distance
    % [NOTE: Needs improvements to account fo artefacts on coarser DEMs]
    for xx = 1:length(uG)
        fpl = FPL.Z(glacierID.Z == uG(xx)); z = DEM.Z(glacierID.Z == uG(xx));
        fpl = max(fpl)-fpl;
        minVal = min(fpl(fpl>-1));
        fpl = fpl + minVal;
        fpl(fpl < 0) = 0;
         
        Outlier = z > prctile(z,20) & fpl > prctile(fpl,75);
        if sum(Outlier) > 0
            ft = polyfit(z(~Outlier),fpl(~Outlier),1);
            fpl(Outlier) = polyval(ft,z(Outlier));
        end
        fpl(fpl < 0) = 0;

        FPL.Z(glacierID.Z == uG(xx)) = fpl;
    end

    FPL.Z(FPL.Z <0) = 0;
    sigma_px = 3.5;

    valid = isfinite(FPL.Z);
    Z = FPL.Z;
    Z(~valid) = 0;
    
    FPL_smooth = FPL;
    FPL_smooth.Z = imgaussfilt(Z, sigma_px) ./ ...
                   imgaussfilt(double(valid), sigma_px);
    
    FPL_smooth.Z(~valid) = NaN;