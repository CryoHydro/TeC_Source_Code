%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%  GRAPHICAL OUTPUT FUNCTION     %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% TETHYS-CHLORIS(T&C)ADVANCED HYDROLOGICAL MODEL%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% TIME SERIES OUTPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load OUTPUT_SIM_Rietholzbach_AVG.dat
%A=OUTPUT_SIM_Rietholzbach_AVG;


Pr = A(:,1); 
Ta = A(:,2); 
Ws= A(:,3);
Ds = A(:,4); 
ea= A(:,5); 
N= A(:,6); 
Pre = A(:,7); 
Tdew = A(:,8); 
Rsw = A(:,9); 
PAR = A(:,10);
Ca = A(:,11);
%%%
Ts =A(:,12); 
Tdamp = A(:,13); 
Csno = A(:,14); 
Cice = A(:,15); 
Csnow = A(:,16); 
Cicew= A(:,17); 
Pr_sno = A(:,18); 
Pr_liq = A(:,19); 
%%%%
Rn = A(:,20); 
H = A(:,21);
G = A(:,22); 
Gfin = A(:,23); 
QE = A(:,24); 
Qv = A(:,25); 
Qfm = A(:,26); 
%%%
SWE = A(:,27); 
SND = A(:,28); 
WR_SP = A(:,29); 
dw_SNO = A(:,30); 
ros = A(:,31); 
In_SWE = A(:,32); 
SP_wc= A(:,33);
ICE= A(:,34);
ICE_D= A(:,35);
WR_IP= A(:,36);
NIce = A(:,37);
IP_wc =  A(:,38);
%%%
T_H= A(:,39); 
T_L= A(:,40);
EIn_H= A(:,41);
EIn_L= A(:,42);
EG= A(:,43);
ESN= A(:,44);
EWAT= A(:,45);
EICE = A(:,46);
Dr_L= A(:,47);
Dr_H=  A(:,48);
%%%%
EIn_urb =  A(:,49);
EIn_rock =  A(:,50);
In= A(:,51);
Inveg =  A(:,52);
WAT =  A(:,53);
FROCK =  A(:,54);
%%%
f= A(:,55);
WIS= A(:,56);
Rd = A(:,57);
Rh = A(:,58);
Lk = A(:,59);
Lk_wat = A(:,60);
Lk_rock = A(:,61);
OF = A(:,62);
OS = A(:,63);
ZWT = A(:,64);
Fract_sat= A(:,65);
%%%
Qlat_in = A(:,66);
Qlat_out= A(:,67);
q_runon= A(:,68);
Q_channel= A(:,69);
V= A(:,70);
O= A(:,71);
Q_exit= A(:,72);
Qsub_exit = A(:,73);
Swe_exit = A(:,74);
%%%
r_soil= A(:,75); 
alp_soil = A(:,76); 
ra = A(:,77);
Tdp = A(:,78);
er= A(:,79);
TsVEG = A(:,80);
DQ_S = A(:,81);
DT_S = A(:,82);
dQ_S = A(:,83);
%%%%
CK1 = A(:,84);
t = A(:,85);
CKt = A(:,86);
%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NN=1:length(CKt); NN=NN';

ET= EG +ESN +T_L +T_H + EIn_H + EIn_L ; 
PLOT_TS=1; 

if PLOT_TS==1 

figure(2)
subplot(2,1,1);
set(gca,'FontSize',11);
plot(NN,OF,'r','LineWidth', 1.5);
hold on; grid on;
plot(NN,OS,'b','LineWidth', 1.5);
%plot(NN,OH,'g','LineWidth', 1.5);
%plot(NN,OL,'--k','LineWidth', 1.5);
title(' \theta  Soil Moisture')
ylabel('[ ]')
legend('Infiltration Layer','Evaporation Layer')
subplot(2,1,2);
hold on; grid on;
%plot(NN,Epot,':y','LineWidth', 1.5);
hold on; grid on;
plot(NN,EG,'b','LineWidth', 1.5);
plot(NN,T_H,'r','LineWidth', 1.5);
plot(NN,T_L,'m','LineWidth', 1.5);
title('Evaporation and Transpiration')
xlabel('Hour'); ylabel('[mm/h]')
%legend('E-Virtual','Ground','H-Vegetation','L-Vegetation')
legend('Ground','H-Vegetation','L-Vegetation')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
subplot(1,1,1);
set(gca,'FontSize',11);
plot(NN,T_H,'r','LineWidth', 1.5);
hold on; grid on;
plot(NN,T_L,'m','LineWidth', 1.5);
title('Moisture Outfluxes')
ylabel('[mm/h]')
plot(NN,EG,':g','LineWidth', 1.5);
ylabel('[mm/h]')
plot(NN,EIn_H+EIn_L,'--y','LineWidth', 1.5);
legend('H-Vegetation','L-Vegetation','Ground Evaporation','Interception Evaporation')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
subplot(2,1,1);
set(gca,'FontSize',11);
plot(NN,Ts,'r','LineWidth', 1.5);
hold on; grid on;
plot(NN,Tdp,'b','LineWidth', 1.5);
plot(NN,Ta(1:length(NN)),':g','LineWidth', 1.5);
title('Temperature')
ylabel('[°C]')
legend('Surface Radiative','Soil at Zdep','Air Temperature')
subplot(2,1,2);
hold on; grid on;
plot(NN,Rn,':k','LineWidth', 1.5);
hold on; grid on;
plot(NN,QE,'g','LineWidth', 1.5);
plot(NN,H,'b','LineWidth', 1.5);
plot(NN,G,'m','LineWidth', 1.5);
plot(NN,Qv,'y','LineWidth', 1.5);
legend('Observed','Simulated')
title('HEAT FLUXES ')
xlabel('Hour'); ylabel('[W/m^2]')
legend('Rn','LE','H','G','Qv')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(11)
set(gca,'FontSize',11);
subplot(3,1,1);
plot(NN,SND,'r','LineWidth', 1.5);
hold on; grid on;
title('Snow Depth');
ylabel('[m]')
subplot(3,1,2);
plot(NN,SWE,'g','LineWidth', 1.5);
hold on; grid on;
title('SWE');
ylabel('[mm]')
subplot(3,1,3);
plot(NN,ros,'g','LineWidth', 1.5);
hold on; grid on;
title('Snow Density');
xlabel('Hour'); ylabel('[kg/m^3]')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(12)
set(gca,'FontSize',11);
subplot(1,2,1);
plot(NN,Q_exit,'k','LineWidth', 1.5);
hold on; grid on;
plot(NN,Qsub_exit,'--g','LineWidth', 1.5);
title('Discharge');
ylabel('[mm/h]')
legend('SUR.','SUB.')
figure(12)
set(gca,'FontSize',11);
subplot(1,2,2);
plot(NN,Fract_sat,'k','LineWidth', 1.5);
hold on; grid on;
title('Saturated Fraction ');
ylabel('[-]')

end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('D:\TIFU\KEMMEexp\MeteoandTopo\Data_Q_Kleine_Emme.mat')
Qobs= Q(:,3)/(100*100*47707)*3600*1000; %%[mm/h]
Qobs=Qobs(6578:end,1); 

figure(105)
subplot(1,1,1)
plot(Q_exit,'b','LineWidth',1.5)
hold on ; grid on
plot(Qobs,'--g','LineWidth',1.5)
legend('Qsim','Qmeas')
xlim([100 30150])

set(gca,'FontSize',11);
subplot(1,2,1);
plot(NN,Q_exit,'k','LineWidth', 1.5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fff=0; 
if fff == 1
    VRtm1=0;
    fr=0.3; 
    %Kres_Rock=30
    Qrock = dth*FROCK/Kres_Rock; %%% [mm]
    %rr=Q_exit;
    rr = Q_exit-Qrock-Rd+Lk_rock*(1-fr);
    %rrs = rr - WR_IP;
    
    k=1/6;
    for i = 1:length(rr); Qoout(i)=k*VRtm1; VR= VRtm1 + rr(i) - k*VRtm1; VR(VR<0)=0; VRtm1=VR  ; end
    Qoout1= Qoout; 
    
    rr=Lk_rock*fr; 
    k=1/8760;
    for i = 1:length(rr); Qoout(i)=k*VRtm1; VR= VRtm1 + rr(i) - k*VRtm1; VR(VR<0)=0; VRtm1=VR  ; end
   
    Qoout2 =Qoout; 
    Qoout = Qoout2 + Rd+ Qoout1;
end

load('D:\T&C_setups\DATA_FOLDER\Discharge_Gletsch.mat') %% 
%01-Oct-1989 01:00:00
Qobs=Qobs(138049:(138049+length(Q_exit)));  
Qobs= Qobs/(50*50*15866)*3600*1000; %%[mm/h]

figure(105)
subplot(1,1,1)
plot(Q_exit,'b','LineWidth',1.5)
hold on ; grid on
plot(Qobs,'--g','LineWidth',1.5)
legend('Qsim','Qmeas')
xlim([1 length(Q_exit)])
plot(WR_SP,'--m','LineWidth',1.5)
plot(WR_IP,'--c','LineWidth',1.5)

n=length(Q_exit);
fr=24;
m=floor(n/fr);
Qs_day=reshape(Q_exit(1:m*fr),fr,m);
Qobs_day=reshape(Qobs(1:m*fr),fr,m);
Qs_day=sum(Qs_day);
Qobs_day=sum(Qobs_day);

Rp=corrcoef(Qs_day,Qobs_day);
R2_Q_day=Rp(1,2)^2;

figure(65)
subplot(1,1,1)
set(gca,'FontSize',12);
plot(1:length(Qs_day),Qobs_day,'xb','LineWidth',2.0);
hold on; grid on
plot(1:length(Qs_day), Qs_day,'r','LineWidth',2.0);
legend('Q_{obs}','Q_{sim}')
xlabel('Time');
title('DISCHARGE')
ylabel('[mm/day]')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Qs= Q_exit'; 
Qs = Qs(61368:end);
load('D:\T&C_setups\DATA_FOLDER\Data_Rietholzbach.mat'); 
x1=61368;1;%219168;
x2=262922;
Qmeas=Q1(x1+1:x2); 
t=Date(x1+1:x2); 



figure(5)
subplot(3,1,1)
plot(t,Qs,'b','LineWidth',1.5)
hold on ; grid on
plot(t,Qmeas,'--g','LineWidth',1.5)
legend('Qsim','Qmeas')
datetick('x',11)
subplot(3,1,2)
semilogy(t,Qs,'b','LineWidth',1.5)
hold on ; grid on
semilogy(t,Qmeas,'--g','LineWidth',1.5)
legend('Qsim','Qmeas')
datetick('x',11)
subplot(3,1,3)
plot(t,cumsum(Qs),'b','LineWidth',1.5)
hold on ; grid on
plot(t,cumsum(Qmeas),'--g','LineWidth',1.5)
legend('Qsim','Qmeas')
datetick('x',11)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[yr]=year(t);
Anno=min(yr):max(yr);
Pyr=zeros(length(Anno),1);
Qobs_yr=zeros(length(Anno),1);
Qsim_yr=zeros(length(Anno),1);
ET_yr=zeros(length(Anno),1);
r=0;
for i=min(Anno):max(Anno)
    r=r+1;
    Pyr(r)= nansum(Pr(find(yr==i)));
    ET_yr(r)= nansum(ET(find(yr==i)));
    Qobs_yr(r)= nansum(Qmeas(find(yr==i)));
    Qsim_yr(r) = nansum(Qs(find(yr==i)));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(6)
subplot(2,1,1)
set(gca,'FontSize',12);
plot(Anno,Pyr,'b','LineWidth',2.0);
hold on; grid on
plot(Anno,ET_yr,'xr','LineWidth',2.0);
plot(Anno,Pyr-Qobs_yr,'xg','LineWidth',2.0);
legend('Pr.','ET_{sim}','ET_{bal}')
xlabel('Water Year');
title('ET')
ylabel('[mm/yr]')
subplot(2,1,2)
set(gca,'FontSize',12);
plot(Anno,Qobs_yr,'xb','LineWidth',2.0);
hold on; grid on
plot(Anno,Qsim_yr,'r','LineWidth',2.0);
legend('Q_{obs}','Q_{sim}')
xlabel('Water Year');
title('DISCHARGE')
ylabel('[mm/yr]')

a=find(ET_yr>0);
Rp=corrcoef(Qobs_yr(a),Qsim_yr(a));
R2_Q_yr=Rp(1,2)^2;


%%%%%%%%%%%%%%%%%%%
n=length(Qs);
fr=24;
m=floor(n/fr);
Qs_day=reshape(Qs(1:m*fr),fr,m);
Qmeas_day=reshape(Qmeas(1:m*fr),fr,m);
Qs_day=sum(Qs_day);
Qmeas_day=sum( Qmeas_day);

Rp=corrcoef(Qs_day,Qmeas_day);
R2_Q_day=Rp(1,2)^2;

figure(65)
subplot(1,1,1)
set(gca,'FontSize',12);
plot(1:length(Qs_day), Qmeas_day,'xb','LineWidth',2.0);
hold on; grid on
plot(1:length(Qs_day), Qs_day,'r','LineWidth',2.0);
legend('Q_{obs}','Q_{sim}')
xlabel('Water Year');
title('DISCHARGE')
ylabel('[mm/day]')


[yr]=year(t);[mo]=month(t);
Anno=min(yr):max(yr);
Pyr=zeros(length(Anno)*12,1);
Qobs_yr=zeros(length(Anno)*12,1);
Qsim_yr=zeros(length(Anno)*12,1);
ET_yr=zeros(length(Anno)*12,1);
r=0;
for i=min(Anno):max(Anno)
    for j=1:12
        r=r+1;
        Pyr(r)= nansum(Pr(intersect(find(yr==i), find(mo==j))));
        ET_yr(r)= nansum(ET(intersect(find(yr==i), find(mo==j))));
        Qobs_yr(r)= nansum(Qmeas(intersect(find(yr==i), find(mo==j))));
        Qsim_yr(r) = nansum(Qs(intersect(find(yr==i), find(mo==j))));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a=find(ET_yr>0);
Rp=corrcoef(Qobs_yr(a),Qsim_yr(a));
R2_Q_mo=Rp(1,2)^2;

figure(7)
subplot(2,1,1)
set(gca,'FontSize',12);
plot(1:length(Anno)*12,Pyr,'b','LineWidth',2.0);
hold on; grid on
plot(1:length(Anno)*12,ET_yr,'xr','LineWidth',2.0);
plot(1:length(Anno)*12,Pyr-Qobs_yr,'xg','LineWidth',2.0);
for ij=1:12*length(Anno)
    if mod(ij,12)==0
        plot([ij ij],[0,max(Pyr)],'--y')
    end
end
legend('Pr.','ET_{sim}','ET_{bal}')
xlabel('Month');
title('ET')
ylabel('[mm/mo]')
subplot(2,1,2)
set(gca,'FontSize',12);
plot(1:length(Anno)*12,Qobs_yr,'xb','LineWidth',2.0);
hold on; grid on
plot(1:length(Anno)*12,Qsim_yr,'r','LineWidth',2.0);
for ij=1:12*length(Anno)
    if mod(ij,12)==0
        plot([ij ij],[0,max(Qobs_yr)],'--y')
    end
end
legend('Q_{obs}','Q_{sim}')
xlabel('Month');
title('DISCHARGE')
ylabel('[mm/mo]')


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%CNS=(1-sum((Qs-Qmeas).^2)/sum((Qs-mean(Qmeas)).^2));
%Rp=corrcoef(Qs,Qmeas); R2=Rp(1,2)^2;
%RMSE=sqrt(sum((Qs-Qmeas).^2)/length(Qmeas)); %% Root mean square error

CNS_h=(1-sum((Qs(not(isnan(Qmeas)))-Qmeas(not(isnan(Qmeas)))).^2)/nansum((Qs(not(isnan(Qmeas)))-mean(Qmeas(not(isnan(Qmeas))))).^2));
Rp=corrcoef(Qs(not(isnan(Qmeas))),Qmeas(not(isnan(Qmeas)))); R2=Rp(1,2)^2;
RMSE=sqrt(sum((Qs(not(isnan(Qmeas)))-Qmeas(not(isnan(Qmeas)))).^2)/length(Qmeas(not(isnan(Qmeas))))); %% Root mean square error

CNS_day=(1-sum((Qs_day(not(isnan(Qmeas_day)))-Qmeas_day(not(isnan(Qmeas_day)))).^2)/nansum((Qs_day(not(isnan(Qmeas_day)))-mean(Qmeas_day(not(isnan(Qmeas_day))))).^2));

%%%%%%%%%%%%%%%
qsim= Qs; %% [mmm/h]
nn=length(qsim);
t=1:nn; t=t/24;
qobs1=Qmeas;
figure(8)
subplot(3,1,1)
set(gca,'FontSize',9);
plot(t,qsim,'r','LineWidth',2.0);
hold on
plot(t,qobs1,'--g','LineWidth',2.0);
grid on
title('Discharge')
legend('SIM.','OBS.')
xlabel('Day');
ylabel('[mm h^{-1}]')
subplot(3,1,2)
set(gca,'FontSize',9);
semilogy(t,qsim,'r','LineWidth',2.0);
hold on
semilogy(t,qobs1,'--g','LineWidth',2.0);
grid on
title('Discharge')
legend('SIM.','OBS.')
xlabel('Day');
ylabel('[mm h^{-1}]')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
I=isnan(qobs1);
qsim(I)=NaN;
q_obs_c = qobs1; q_obs_c(isnan(q_obs_c))=0;
q_sim_c = qsim;  q_sim_c(isnan(q_sim_c))=0;
subplot(3,1,3)
set(gca,'FontSize',9);
plot(t,cumsum(q_sim_c),'r','LineWidth',2.0);
hold on
plot(t,cumsum(q_obs_c),'--b','LineWidth',2.0);
grid on
title('Cumulative Discharge')
legend('SIM.','OBS.')
xlabel('Day');
ylabel('[mm]')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
th=1e-5;
cf=1;  %% [mm/h] to [lt/s]
nv=length(qsim) - sum(isnan(qsim));
dur= 365*((1:nv)/(nv+1));
p1 = -cf*sort(-qsim);
p2 = -cf*sort(-qobs1);
p1(p1<th)=NaN;
p2(p2<th)=NaN;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(9)
subplot(2,1,1)
semilogy(dur,p1(1:nv),'r','LineWidth',2.0);
hold on; grid on
semilogy(dur,p2(1:nv),'--b','LineWidth',2.0);
legend('SIM.','OBS.')
xlabel('Duration');
title('Duration Curve')
ylabel('[mm/h]')
xlim([0 366])
%ylim([0.001 100000])
figure(9)
subplot(2,1,2)
loglog(dur,p1(1:nv),'r','LineWidth',2.0);
hold on; grid on
loglog(dur,p2(1:nv),'--b','LineWidth',2.0);
legend('SIM.','OBS.')
xlabel('Duration');
title('Duration Curve')
ylabel('[mm/h]')
xlim([0 366])
%ylim([0.001 100000])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Qsim_F = nansum(qsim)/( sum(not(isnan(qobs1)))/8765);
Qobs_F = nansum(qobs1)/( sum(not(isnan(qobs1)))/8765);


figure(101)
subplot(3,1,3)
set(gca,'FontSize',11.5);
semilogy(dur,p1(1:nv),'r','LineWidth',2.0);
hold on; grid on
semilogy(dur,p2(1:nv),'--b','LineWidth',2.0);
legend('SIM.','OBS.')
xlabel('Duration');
title('Duration Curve')
ylabel('[mm/h]')
xlim([0 366])
ylim([1e-5 1e+2])

subplot(3,1,2)
set(gca,'FontSize',11.5);
plot(1:length(Anno)*12,Qobs_yr,'xb','LineWidth',2.0);
hold on; grid on
plot(1:length(Anno)*12,Qsim_yr,'r','LineWidth',2.0);
for ij=1:12*length(Anno)
    if mod(ij,12)==0
        plot([ij ij],[0,max(Qobs_yr)],'--y')
    end
end
legend('Q_{obs}','Q_{sim}')
xlabel('Month');
title('Monthly Q')
ylabel('[mm/mo]')
%xlim([1 275])

subplot(3,1,1)
set(gca,'FontSize',11.5);
plot(t,qsim,'r','LineWidth',2.0);
hold on
plot(t,qobs1,'--g','LineWidth',2.0);
grid on
title('Hourly Q')
legend('SIM.','OBS.')
xlabel('Day');
ylabel('[mm h^{-1}]')
%xlim([5800 6000])






LAI_H=A(:,9);
LAI_L=A(:,10);
NPP_H=A(:,11);
NPP_L=A(:,12);
ANPP_H=A(:,13);
ANPP_L=A(:,14);
RA_H=A(:,15);
RA_L=A(:,16);