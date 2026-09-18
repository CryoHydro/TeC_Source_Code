%Code to pull out the Ta lapse rate for the Shallap catchment

%Shallap is in catchment 406 with main station 144.

addpath('C:\Users\CFyffe\Dropbox\Northumbria_University\PeruGROWS\TOPKAPI\MODEL_Cat\FORCINGS');

%load data
Ta_lapseT = readtable('Ta_lapse_RS.csv');

SH_Ta_lapse = Ta_lapseT(:,[1 8]); %Column 8 is for station 144

SH_Ta_lapse.hour = hour(SH_Ta_lapse.Time);
SH_Ta_lapse.month = month(SH_Ta_lapse.Time);
SH_Ta_lapse=movevars(SH_Ta_lapse,'hour','After','Time');
SH_Ta_lapse=movevars(SH_Ta_lapse,'month','After','hour');

SH_hm_lapse = grpstats(SH_Ta_lapse,["hour","month"],"mean"); %We don't really need the average since the values are just repeated over time, but it allows us to extract them

%Export table
save('SH_hm_Ta_lapse.mat',"SH_hm_lapse");