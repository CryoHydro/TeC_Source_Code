function[SLE,SLnoise]=Snowline(DTM,SND,In_SWE,LAI_H,Aspect,OPT_thru_vis)
% function to compute the snow line elevation, based on the minimization of
% snow pixels below/land pixels above a given elevation. 
% based on Krajci et al.(2016)

% Possibility to extend to calculate aspect-dependent snowlines


demmin = min(DTM,[],'all'); %minimum elevation
demmax = max(DTM,[],'all'); %maximum elevation
SLE = nan;
SLnoise = 10000000000;
SNDadj = SND; SNDadj(isnan(DTM)) = nan;
if OPT_thru_vis == 1
%assumption that satellites cannot see 'throughfall-snow':
%adjust 'visible' snow cover in cases where high vegetation is present but no snow intercepted
SNDadj(LAI_H>0 & In_SWE < 0.5) = 0; 
end
nl = length(find(SNDadj<=0.005 & SNDadj >= 0)); %number of land pixels without snow

if ~isempty(SNDadj>0)
for elev = demmin:10:demmax %loop through elevation profile by 10m
tmp = SNDadj; tmp(isnan(DTM))=nan; 
tmp((DTM>elev)==1) = nan;

nsp = length(find(tmp>0.005)); %number of snow pixels below elev
nlp = nl-length(find(tmp<=0.005 & tmp >= 0)); %number of land pixels above elev

if(nsp+nlp <= SLnoise) %minimize sum of nsp and nlp
    SLnoise = nlp+nsp;
    if elev == max(demmin:10:demmax)
      SLE = nan;
    else
      SLE = elev;
    end
end
end

end