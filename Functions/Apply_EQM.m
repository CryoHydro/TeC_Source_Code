function [BCmain] = Apply_EQM(DATE,EST,coeff,VAR)

% BCmain = nan(size(EST));
nTimesteps = size(EST);

qs = size(coeff,1);

%% DAILY AGGREGATE
    if strcmp(VAR,'PP')
          for dd = 1:length(DATE)/24
                    DAY_EXTRACT = (24*dd - 23):(24*dd);
                     EST_DAY(dd) = sum(EST(DAY_EXTRACT));
                        for hh = 1:24
                          HourDist(DAY_EXTRACT(hh)) = EST(DAY_EXTRACT(hh)) / sum(EST(DAY_EXTRACT));
                        end
          end
          DATE_DAY = DATE(1:24:end-23);
    end

%% APPLY BC
if size(coeff,2) > 1
    for MM = 1:12

        M_IDX = month(DATE) == MM;
        if sum(M_IDX) > 0 
            MLOC = find(M_IDX);
            ESTm = EST(MLOC);

            ESTQm(:,MM) = quantile(ESTm,qs);    

        end
    end
    
    
        if strcmp(VAR,'PP')
            for iTimestep = 1:length(EST_DAY)
                [~,indQ] = min(abs(EST_DAY(iTimestep)-ESTQm(:,month(DATE_DAY(iTimestep)))));
                BCmain(iTimestep) = EST_DAY(iTimestep)*coeff(indQ,month(DATE_DAY(iTimestep)));
            end
        else
            for iTimestep = 1:nTimesteps
                [~,indQ] = min(abs(EST(iTimestep)-ESTQm(:,month(DATE(iTimestep)))));
                BCmain(iTimestep) = EST(iTimestep)-coeff(indQ,month(DATE(iTimestep)));    
            end
        end   
    
    
else
    
ESTQ = quantile(EST,qs);     
        if strcmp(VAR,'PP')
            for iTimestep = 1:length(EST_DAY)
                [~,indQ] = min(abs(EST_DAY(iTimestep)-ESTQ));
                BCmain(iTimestep) = EST_DAY(iTimestep)*coeff(indQ);
            end
        else
            for iTimestep = 1:nTimesteps
                [~,indQ] = min(abs(EST(iTimestep)-ESTQ));
                BCmain(iTimestep) = EST(iTimestep)-coeff(indQ);    
            end
        end   


end

if strcmp(VAR,'PP')
  for dd = 1:length(EST)/24    
      DAY_EXTRACT = (24*dd - 23):(24*dd);
        for hh = 1:24
            BCmain_HR(DAY_EXTRACT(hh)) = BCmain(dd) *  HourDist(DAY_EXTRACT(hh));
        end
  end  
   
  BCmain = BCmain_HR;
end   



