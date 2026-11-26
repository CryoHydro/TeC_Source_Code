%Quick script to test the radiation partition

addpath('C:\Users\CFyffe\Dropbox\Northumbria_University\PeruGROWS\Modelling\T_and_C_HIMAL\TC_PrePost\PREprocessing\Cat');
%To try other scripts

x1 = datetime(2014,11,1,0,0,0);
x2 = datetime(2014,11,1,23,0,0);
Datetime=x1:hours(1):x2;

Date = datenum(Datetime);
Lat = -9.4900;
Lon = -77.3436;
Zbas = 4767;
DeltaGMT = -5;
Pr = [0 1 0 1 1 1 1 1 1 1 1.2 1.5 1 1 0 0 1.2 1.5 1 1 0 0 0 0]';
Tdew = [2 3 4 2 3 4 2 3 4 2 3 4 2 3 4 2 3 4 2 3 4 2 3 4]';
Tdew =-Tdew;
Rsw = [0 0 0 0 0 0 100 200 300 400 400 500 700 700 900 900 800 600 500 300 200 100 0 0]';
t_bef = 1;
t_aft = 0;
GRAPH_VAR = 0;

%Note should find t_bef and t_aft once first (with previous version of
%automatic partition) and then can use the same values for other times

 [SAD1,SAD2,SAB1,SAB2,PARB,PARD, A,Rsws,Nsim,Prw,Tatm1,Tatm2] = SH_Automatic_Radiation_Partition_fast(Date,Lat,Lon,Zbas,...
     DeltaGMT,Pr,Tdew,Rsw,t_bef,t_aft);

%[SD,SB,SAD1,SAD2,SAB1,SAB2,PARB,PARD,N,Rsws,t_bef,t_aft]=Peru_Automatic_Radiation_Partition_I(Date,Lat,Lon,Zbas,DeltaGMT,Pr,Tdew,Rsw,GRAPH_VAR);

%Ok so it works fine with a time series, but get stuck on a scalar on
% Nn=interp1(xnt,Nsim(not(isnan(Nsim))),xn,'linear'); line 258 / 147
% depending on routine