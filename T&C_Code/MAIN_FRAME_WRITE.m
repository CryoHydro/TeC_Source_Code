%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% MAIN_FRAME OF HBM-VEG %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INITIZIALIZATION VARIABLES
%%%%%--->>  INIZIALIZATION VARIABLES  <<--- %%%%%%%%%%%%%%%%%%%%%%%
%%% j time dt = 1 day  %%%%%%%%%%%%%%%%
%%% i time dt = 1h
%%% ms soil layer
%%% cc Crown Area number
%%% NN time step
dtd = 1; %% [day]
dth = dt/3600; %% [hour]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NCP = 7; %% Number of Carbon Pool
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% CALL PARAMETERS AND INITIAL CONDITION
In_H=zeros(1,cc);In_L=zeros(1,cc);
Ci_sunH=zeros(1,cc);Ci_sunL=zeros(1,cc);
Ci_shdH=zeros(1,cc);Ci_shdL=zeros(1,cc);
%%%%
LAI_L=zeros(1,cc); B_L=zeros(1,cc,NCP);NPP_L=zeros(1,cc);Rg_L=zeros(1,cc);
RA_L=zeros(1,cc);Rmc_L=zeros(1,cc);ANPP_L=zeros(1,cc);
Rms_L=zeros(1,cc);Rmr_L=zeros(1,cc); PHE_S_L=zeros(1,cc); dflo_L=zeros(1,cc);
AgeL_L=zeros(1,cc);e_rel_L=zeros(1,cc);SAI_L=zeros(1,cc); hc_L=zeros(1,cc);
LAI_H=zeros(1,cc); B_H=zeros(1,cc,NCP);NPP_H=zeros(1,cc);Rg_H=zeros(1,cc);
RA_H=zeros(1,cc);Rms_H=zeros(1,cc);Rmr_H=zeros(1,cc); ANPP_H=zeros(1,cc);
PHE_S_H=zeros(1,cc); dflo_H=zeros(1,cc);Rmc_H=zeros(1,cc);
AgeL_H=zeros(1,cc);e_rel_H=zeros(1,cc);SAI_H=zeros(1,cc); hc_H=zeros(1,cc);
%%%
Sr_H=zeros(1,cc);  Slf_H=zeros(1,cc);
Sfr_H=zeros(1,cc); Swm_H=zeros(1,cc);
Sr_L=zeros(1,cc); Slf_L=zeros(1,cc);
Sfr_L=zeros(1,cc); Swm_L=zeros(1,cc);
LAIdead_H=zeros(1,cc); LAIdead_L=zeros(1,cc);
Llitter=zeros(1,cc);
%%%%%%%%%%%%%%%%%
run(PARAM_IC);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Check_Land_Cover_Fractions(Crock,Curb,Cwat,Cbare,Ccrown);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B_H=squeeze(B_H);
B_L=squeeze(B_L);
if cc == 1 
    B_H=B_H';
    B_L=B_L';
end
%%%% Lateral Contribution
q_runon=0; %%[mm/h]
Qi_in=zeros(1,ms); %%[mm/h]
%%%%%%%%%%%%%%%%%%%%%%%%%
tic ;
%profile on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bau = waitbar(0,'Waiting...');
j=1; %%% time dt = 1 day  %%%%%%%%%%%%%%%%
Tstm0= Ts;
%%% i time dt = 1h
%%% ms soil layer
%%% cc Crown Area present
%%%%%%%%%%%%%%%% NUMERICAL METHODS OPTIONS
%Opt_CR = optimset('TolX',3);
%Opt_ST = optimset('TolX',0.1);
%Opt_ST = optimset('TolX',0.1,'Display','iter');
Opt_CR = optimset('TolFun',1);%,'UseParallel','always');
Opt_ST = optimset('TolFun',0.1);%,'UseParallel','always');
OPT_SM=  odeset('AbsTol',0.05,'MaxStep',dth);
OPT_VD=  odeset('AbsTol',0.05);
OPT_STh = odeset('AbsTol',5e+3);
OPT_VegSnow = 0;
OPT_SoilTemp = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Vtm1= V;
Tstm1 = Ts;
SNDtm1 = SND;
snow_albtm1 = snow_alb;
Citm1_sunH = Ci_sunH;
Citm1_shdH = Ci_shdH;
Citm1_sunL = Ci_sunL;
Citm1_shdL = Ci_shdL;
e_snotm1 = e_sno;
In_Htm1 =In_H;
In_Ltm1 =In_L;
In_urbtm1 = In_urb;
In_rocktm1 = In_rock;
SWEtm1= SWE;
In_SWEtm1 = In_SWE;
Sdptm1 = Sdp;
Tdptm1 = Tdp;
Tdamptm1 = Tdamp;
WATtm1=WAT;
ICEtm1=ICE;
IP_wctm1=IP_wc;
ICE_Dtm1=ICE_D;
Cicewtm1= 0 ;%Cicew;
FROCKtm1=FROCK;
t_slstm1= t_sls;
rostm1=ros;
SP_wctm1=SP_wc ;
tau_snotm1 = tau_sno;
EKtm1 = EK;
LAI_Htm1=LAI_H;
B_Htm1= B_H;
PHE_S_Htm1=PHE_S_H;
dflo_Htm1= dflo_H;
AgeL_Htm1= AgeL_H;
SAI_Htm1= SAI_H;
hc_Htm1= hc_H;
LAI_Ltm1= LAI_L;
B_Ltm1= B_L;
PHE_S_Ltm1= PHE_S_L;
dflo_Ltm1= dflo_L ;
AgeL_Ltm1= AgeL_L ;
SAI_Ltm1= SAI_L;
hc_Ltm1= hc_L;
Llittertm1 = Llitter;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZATION
Ta_t=zeros(1,24);
Tdp_H_t=zeros(cc,720);
Tdp_L_t=zeros(cc,720);
OS_t=zeros(1,24);
OH_t=zeros(cc,168);
OL_t=zeros(cc,168);
An_H_t=zeros(cc,24);
An_L_t=zeros(cc,24);
Rdark_H_t=zeros(cc,24);
Rdark_L_t=zeros(cc,24);
NPP_H_t=zeros(cc,7);
NPP_L_t=zeros(cc,7);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:NN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if  (mod(i,1000) == 0) || (i == 2)
        %waitbar(i/NN,bau)
        disp('Iter:'); disp(i);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Datam(i,4)==1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        j=j+1; [jDay]= JULIAN_DAY(Datam(i,:));
        [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day]= SetSunVariables(Datam(i,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
        clear h_S delta_S zeta_S T_sunrise T_sunset
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for cc=1:length(Ccrown)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ZR_H(cc) > 0
                OsatH = sum(RfH_Zs(cc,:).*Osat);
                %%%
                [LAI_H(cc),B_H(cc,:),NPP_H(cc),ANPP_H(cc),Rg_H(cc),RA_H(cc),Rms_H(cc),Rmr_H(cc),Rmc_H(cc),PHE_S_H(cc),...
                    dflo_H(cc),AgeL_H(cc),e_rel_H(cc),LAIdead_H(cc),Sr_H(cc),Slf_H(cc),Sfr_H(cc),Swm_H(cc),Rexmy_H(cc)]= VEGGIE_UNIT_S(LAI_Htm1(cc),...
                    B_Htm1(cc,:),PHE_S_Htm1(cc),dflo_Htm1(cc),AgeL_Htm1(cc),...
                    Ta_t,Tdp_H_t(cc,:),OH_t(cc,:),An_H_t(cc,:),Rdark_H_t(cc,:),OS_t,Oss_H(cc),Owp_H(cc),OsatH,i,j,NPP_H_t(cc,:),jDay,...
                    Sl_H(cc),mSl_H(cc),Ns_H(cc),Nr_H(cc),r_H(cc),gR_H(cc),LtR_H(cc),eps_ac_H(cc),aSE_H(cc),Trr_H(cc),...
                    dd_max_H(cc),dc_C_H(cc),Tcold_H(cc),drn_H(cc),dsn_H(cc),age_cr_H(cc),Bfac_lo_H(cc),Bfac_ls_H(cc),...
                    Tlo_H(cc),Tls_H(cc),mjDay_H(cc),LDay_min_H(cc),dmg_H(cc),...
                    Mf_H(cc),Wm_H(cc),LAI_min_H(cc),L_day,LDay_cr_H(cc),...
                    fab_H(cc),fbe_H(cc),Klf_H(cc),ff_r_H(cc),dexmy_H(cc),jDay_cut_H(cc,:),LAI_cut_H(cc),Lat,OPT_VD);
                %%%%%%%%%%%%%%%%%%%
                SAI_H(cc) =SAI_Htm1(cc);
                if aSE_H(cc) == 2
                    [hc_H(cc)] = GrassHeight(LAI_H(cc),LAIdead_H(cc));
                else
                    hc_H(cc)= hc_Htm1(cc); %%%[m]
                end
            else
                LAI_H(cc)=0;B_H(cc,:)=0;NPP_H(cc)=0;ANPP_H(cc)=0;Rg_H(cc)=0;RA_H(cc)=0;Rms_H(cc)=0;
                Rmr_H(cc)=0;PHE_S_H(cc)=0;dflo_H(cc)=0;AgeL_H(cc)=0;e_rel_H(cc)=0;
                SAI_H(cc)=0; hc_H(cc)=0;
                LAIdead_H(cc) =0;
                Sr_H(cc)=0;Slf_H(cc)=0;Sfr_H(cc)=0;Swm_H(cc)=0;
            end
            %%%%%%%%%%%%%%%%%
            if ZR_L(cc) > 0
                OsatL = sum(RfL_Zs(cc,:).*Osat);
                %%%
                [LAI_L(cc),B_L(cc,:),NPP_L(cc),ANPP_L(cc),Rg_L(cc),RA_L(cc),Rms_L(cc),Rmr_L(cc),Rmc_L(cc),PHE_S_L(cc),...
                    dflo_L(cc),AgeL_L(cc),e_rel_L(cc),LAIdead_L(cc),Sr_L(cc),Slf_L(cc),Sfr_L(cc),Swm_L(cc),Rexmy_L(cc)]= VEGGIE_UNIT_S(LAI_Ltm1(cc),...
                    B_Ltm1(cc,:),PHE_S_Ltm1(cc),dflo_Ltm1(cc),AgeL_Ltm1(cc),...
                    Ta_t,Tdp_L_t(cc,:),OL_t(cc,:),An_L_t(cc,:),Rdark_L_t(cc,:),OS_t,Oss_L(cc),Owp_L(cc),OsatL,i,j,NPP_L_t(cc,:),jDay,...
                    Sl_L(cc),mSl_L(cc),Ns_L(cc),Nr_L(cc),r_L(cc),gR_L(cc),LtR_L(cc),eps_ac_L(cc),aSE_L(cc),Trr_L(cc),...
                    dd_max_L(cc),dc_C_L(cc),Tcold_L(cc),drn_L(cc),dsn_L(cc),age_cr_L(cc),Bfac_lo_L(cc),Bfac_ls_L(cc),...
                    Tlo_L(cc),Tls_L(cc),mjDay_L(cc),LDay_min_L(cc),dmg_L(cc),...
                    Mf_L(cc),Wm_L(cc),LAI_min_L(cc),L_day,LDay_cr_L(cc),...
                    fab_L(cc),fbe_L(cc),Klf_L(cc),ff_r_L(cc),dexmy_L(cc),jDay_cut_L(cc,:),LAI_cut_L(cc),Lat,OPT_VD);
                %%%%%%%%%%%%%%%%%%%
                SAI_L(cc) =SAI_Ltm1(cc);
                if aSE_L(cc) == 2
                    [hc_L(cc)] = GrassHeight(LAI_L(cc),LAIdead_L(cc));
                else
                    hc_L(cc)= hc_Ltm1(cc); %%%[m]
                end
            else
                LAI_L(cc)=0;B_L(cc,:)=0;NPP_L(cc)=0;ANPP_L(cc)=0;Rg_L(cc)=0;RA_L(cc)=0;Rms_L(cc)=0;
                Rmr_L(cc)=0;PHE_S_L(cc)=0;dflo_L(cc)=0;AgeL_L(cc)=0;e_rel_L(cc)=0;
                SAI_L(cc)=0; hc_L(cc)=0;
                LAIdead_L(cc) =0;
                Sr_L(cc)=0;Slf_L(cc)=0;Sfr_L(cc)=0;Swm_L(cc)=0;
            end
            %%%%%%%%%%%%%%%%%
            Llitter(cc)=Llittertm1(cc);
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [V,O,ZWT,OF,OS,OH,OL,Rd,Qi_out,WTR,...
    Rh,Lk,f,WIS,Ts,Pr_sno,Pr_liq,Csno,Cice,NDVI,rb_H,rb_L,rs_sunH,...
    rs_sunL,rs_shdH,rs_shdL,r_litter,...
    An_L,An_H,Rdark_L,Rdark_H,Ci_sunH,Ci_sunL,Ci_shdH,Ci_shdL,...
    rap_H,rap_L,r_soil,b_soil,alp_soil,ra,Rn,...
    H,QE,Qv,T_H,T_L,EIn_H,EIn_L,EG,ESN,ESN_In,ELitter,EWAT,EICE,EIn_urb,EIn_rock,dw_SNO,...
    G,Gfin,Tdp,Sdp,Tdamp,Tdp_H,Tdp_L,SWE,SND,ros,In_SWE,SP_wc,WR_SP,U_SWE,NIn_SWE,dQ,Qfm,t_sls,DQ,DT,...
    WAT,ICE,ICE_D,IP_wc,WR_IP,NIce,Cicew,Csnow,FROCK,...
    In_H,In_L,In_Litter,In_urb,In_rock,Dr_H,Dr_L,SE_rock,SE_urb,Lk_wat,Lk_rock,er,...
    gsr_H,Psi_x_H,Psi_l_H,Jsx_H,Jxl_H,Kleaf_H,Kx_H,Vx_H,Vl_H,...
    gsr_L,Psi_x_L,Psi_l_L,Jsx_L,Jxl_L,Kleaf_L,Kx_L,Vx_L,Vl_L,...
    snow_alb,tau_sno,e_sno,Ws_under,dQVEG,TsV,EK,POT]=HYDROLOGIC_UNIT(Vtm1,aR,Zs,...
    EvL_Zs,Inf_Zs,Zinf,RfH_Zs,RfL_Zs,dz,Ks_Zs,Dz,ms,Kbot,Pr(i),Ta(i),Ds(i),Ws(i),zatm,Tstm1,dt,dth,ea(i),N(i),Pre(i),Tstm0,...
    LAI_H,SAI_H,LAI_L,SAI_L,LAIdead_H,LAIdead_L,BLit,Sllit,Kct,...
    Sp_SN_In,Sp_LAI_L_In,Sp_LAI_H_In,Datam(i,:),DeltaGMT,Lon,Lat,t_bef,t_aft,...
    Ccrown,Cbare,Crock,Curb,Cwat,...
    SAB1(i),SAB2(i),SAD1(i),SAD2(i),PARB(i),PARD(i),SvF,SNDtm1,snow_albtm1,Color_Class,OM_H,OM_L,...
    PFT_opt_H,PFT_opt_L,hc_H,hc_L,d_leaf_H,d_leaf_L,...
    Osat,Ohy,L,Pe,O33,alpVG,nVG,...
    gcI,KcI,TminS,TmaxS,...
    KnitH,KnitL,mSl_H,Sl_H,mSl_L,Sl_L,Ca(i),Oa,Citm1_sunH,Citm1_shdH,CT_H,Citm1_sunL,Citm1_shdL,CT_L,Vmax_H,Vmax_L,FI_H,FI_L,a1_H,go_H,a1_L,go_L,...
    DSE_H,Ha_H,Do_H,e_rel_H,DSE_L,Ha_L,Do_L,e_rel_L,gmes_H,rjv_H,gmes_L,rjv_L,...
    e_snotm1,In_Htm1,In_Ltm1,In_Littertm1,In_urbtm1,In_rocktm1,SWEtm1,In_SWEtm1,....
    Sdptm1,Tdptm1,Tdamptm1,rsd,lan_dry,lan_s,cv_s,...
    WATtm1,ICEtm1,IP_wctm1,ICE_Dtm1,Ice_wc_sp,ros_Ice_thr,Aice,Cicewtm1,...
    Vxtm1_H,Vltm1_H,Vxtm1_L,Vltm1_L,Psi_xtm1_H,Psi_ltm1_H,Psi_xtm1_L,Psi_ltm1_L,...
    Psi_sto_50_H,Psi_sto_99_H,Rrootl_H,rcyl_H,rroot_H,ZR95_H,...
    Psi_sto_50_L,Psi_sto_99_L,Rrootl_L,rcyl_L,rroot_L,ZR95_L,...
    Axyl_H,PsiL50_H,PsiL99_H,Kleaf_max_H,Cl_H,Kx_max_H,PsiX50_H,Cx_H,...
    Axyl_L,PsiL50_L,PsiL99_L,Kleaf_max_L,Cl_L,Kx_max_L,PsiX50_L,Cx_L,...
    FROCKtm1,Krock,Ws_undertm1,...
    Tdew(i),t_slstm1,rostm1,SP_wctm1,fpr,WatFreez_Th,dz_ice,...
    Th_Pr_sno,ros_max1,ros_max2,...
    In_max_urb,In_max_rock,K_usle,tau_snotm1,Slo_top,Slo_pot,Asur,Ared,aTop,EKtm1,q_runon,Qi_in,...
    pow_dis,a_dis,...
    SPAR,SN,OPT_VegSnow,OPT_SoilTemp,Opt_CR,Opt_ST,OPT_SM,OPT_STh)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    STOT= SAB1(i)+SAB2(i)+SAD1(i)+SAD2(i);
    ALB = SAB1(i)/STOT*snow_alb.dir_vis + SAD1(i)/STOT*snow_alb.dif_vis + ...
        SAB2(i)/STOT*snow_alb.dir_nir + SAD2(i)/STOT*snow_alb.dif_nir;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% v-coordinate
    CK1 = f*dth*Asur*Ared + sum(Vtm1 - V)*Asur*Ared - EG*dth - Lk*dth ...
        - sum(Qi_out)*dth -Rd -sum(T_L).*dth -sum(T_H).*dth  + sum(Qi_in)*dth  ;
    CK2 = Pr_sno*dth + Pr_liq*dth*max(Csno,Csnow) ...
        - ESN*dth - ESN_In*dth - WR_SP*dth*Asur + SWEtm1-SWE ...
        + (In_SWEtm1 - In_SWE) + (SP_wctm1 -SP_wc) -sum(Dr_H)*Csno -sum(Dr_L)*Csno  + ...
        sum(In_Htm1 - In_H).*Csno + sum(In_Ltm1 - In_L).*Csno;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% INITIAL CONDITION FOR THE NEXT STEP        %%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Tstm0 =2*Ts-Tstm1;
    Vtm1=V;
    Tstm1 =Ts;
    SNDtm1 = SND ;
    snow_albtm1 = snow_alb ;
    Citm1_sunH = Ci_sunH;
    Citm1_sunL = Ci_sunL;
    Citm1_shdH = Ci_shdH;
    Citm1_shdL = Ci_shdL;
    e_sontm1 = e_sno;
    In_Htm1 = In_H ;
    In_Ltm1 = In_L;
    In_urbtm1 = In_urb;
    In_rocktm1 = In_rock;
    SWEtm1 = SWE ;
    In_SWEtm1 = In_SWE;
    Tdamptm1 = Tdamp;
    Tdptm1 = Tdp;
    Sdptm1 = Sdp;
    WATtm1=WAT;
    ICEtm1=ICE;
    IP_wctm1=IP_wc;
    ICE_Dtm1=ICE_D;
    Cicewtm1=Cicew;
    FROCKtm1=FROCK;
    t_slstm1= t_sls;
    rostm1=ros;
    SP_wctm1=SP_wc ;
    tau_snotm1 = tau_sno;
    EKtm1 = EK;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Datam(i,4)==1)
        LAI_Htm1=LAI_H;
        B_Htm1= B_H;
        PHE_S_Htm1=PHE_S_H;
        dflo_Htm1= dflo_H;
        AgeL_Htm1= AgeL_H;
        SAI_Htm1= SAI_H;
        hc_Htm1= hc_H;
        LAI_Ltm1= LAI_L;
        B_Ltm1= B_L;
        PHE_S_Ltm1= PHE_S_L;
        dflo_Ltm1= dflo_L ;
        AgeL_Ltm1= AgeL_L ;
        SAI_Ltm1= SAI_L;
        hc_Ltm1= hc_L;
        %%%
        Llitter_tm1 = Llitter;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% MEMORY CONDITION FOR VEGETATION  MODEL    %%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% 1 day
    if i > 24
        Ta_t(:,1:23)=Ta_t(:,2:24);
        Ta_t(:,24)=Ta(i);
        An_H_t(:,1:23)=An_H_t(:,2:24);
        An_H_t(:,24)=An_H;
        An_L_t(:,1:23)=An_L_t(:,2:24);
        An_L_t(:,24)=An_L;
        Rdark_H_t(:,1:23)=Rdark_H_t(:,2:24);
        Rdark_H_t(:,24)=Rdark_H;
        Rdark_L_t(:,1:23)=Rdark_L_t(:,2:24);
        Rdark_L_t(:,24)=Rdark_L;
        OS_t(:,1:23)=OS_t(:,2:24);
        OS_t(:,24)=OS;
    else
        Ta_t(:,i)=Ta(i);
        An_H_t(:,i)=An_H;
        An_L_t(:,i)=An_L;
        Rdark_H_t(:,i)=Rdark_H;
        Rdark_L_t(:,i)=Rdark_L;
        OS_t(:,i)=OS;
    end
    %%%% 7 day
    if i > 168
        OH_t(:,1:167)= OH_t(:,2:168);
        OH_t(:,168)= OH;
        OL_t(:,1:167)= OL_t(:,2:168);
        OL_t(:,168)= OL;
    else
        OH_t(:,i)=OH;
        OL_t(:,i)=OL;
    end
    %%%% 30 day
    if i > 720
        Tdp_H_t(:,1:719)= Tdp_H_t(:,2:720);
        Tdp_H_t(:,720)= Tdp_H;
        Tdp_L_t(:,1:719)= Tdp_L_t(:,2:720);
        Tdp_L_t(:,720)= Tdp_L;
    else
        Tdp_H_t(:,i)=Tdp_H;
        Tdp_L_t(:,i)=Tdp_L;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Datam(i,4)==1)
        if j > 7
            NPP_H_t(:,1:6)= NPP_H_t(:,2:7);
            NPP_L_t(:,1:6)= NPP_L_t(:,2:7);
            NPP_H_t(:,7) = NPP_H;
            NPP_L_t(:,7) = NPP_L;
        else
            NPP_H_t(:,j)=NPP_H;  %%% 7 day
            NPP_L_t(:,j)=NPP_L;  %%% 7 day
        end
    end
    %%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%% WRITING  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if i==2
        tit5=strcat('OUTPUT_',TITLE_SAVE,'_VAR.dat');
        fid5=fopen(tit5,'a');
    end
    %%%%
    fprintf(fid5,'%g\t',i);%%11
    fprintf(fid5,'%g\t',Ts);%12
    fprintf(fid5,'%g\t',Tdamp);%13
    fprintf(fid5,'%g\t',Csno);%14
    fprintf(fid5,'%g\t',Cice);%15
    fprintf(fid5,'%g\t',Csnow);%16
    fprintf(fid5,'%g\t',Cicew);%17
    fprintf(fid5,'%g\t',Pr_sno);%18
    fprintf(fid5,'%g\t',Pr_liq);%19
    %%%%
    fprintf(fid5,'%g\t',Rn);%20
    fprintf(fid5,'%g\t',H);%21
    fprintf(fid5,'%g\t',G);%22
    fprintf(fid5,'%g\t',Gfin);%23
    fprintf(fid5,'%g\t',QE);%24
    fprintf(fid5,'%g\t',Qv);%25
    fprintf(fid5,'%g\t',Qfm);%26
    %%%%
    fprintf(fid5,'%g\t',SWE);%27
    fprintf(fid5,'%g\t',SND);%28
    fprintf(fid5,'%g\t',WR_SP);%29
    fprintf(fid5,'%g\t',U_SWE);%30
    fprintf(fid5,'%g\t',NIn_SWE);%31
    fprintf(fid5,'%g\t',dw_SNO);%32
    fprintf(fid5,'%g\t',ros);%33
    fprintf(fid5,'%g\t',In_SWE);%34
    fprintf(fid5,'%g\t',SP_wc);%35
    fprintf(fid5,'%g\t',ICE);%36
    fprintf(fid5,'%g\t',ICE_D);%37
    fprintf(fid5,'%g\t',WR_IP);%38
    fprintf(fid5,'%g\t',IP_wc);%39
    fprintf(fid5,'%g\t',NIce);%40
    %%%
    fprintf(fid5,'%g\t',EG);%41
    fprintf(fid5,'%g\t',ESN);%42
    fprintf(fid5,'%g\t',ESN_In);%43
    fprintf(fid5,'%g\t',EWAT);%44
    fprintf(fid5,'%g\t',EICE);%45
    fprintf(fid5,'%g\t',EIn_urb);%46
    fprintf(fid5,'%g\t',EIn_rock);%47
    fprintf(fid5,'%g\t',SE_rock);%48
    fprintf(fid5,'%g\t',SE_urb); %49
    fprintf(fid5,'%g\t',In_urb);%50
    fprintf(fid5,'%g\t',In_rock); %51
    fprintf(fid5,'%g\t',WAT); %52
    fprintf(fid5,'%g\t',FROCK); %53
    %%%
    fprintf(fid5,'%g\t',f);%54
    fprintf(fid5,'%g\t',WIS);%55
    fprintf(fid5,'%g\t',Rd);%56
    fprintf(fid5,'%g\t',Rh);%57
    fprintf(fid5,'%g\t',Lk);%58
    fprintf(fid5,'%g\t',Lk_wat);%59
    fprintf(fid5,'%g\t',Lk_rock);%60
    %%%
    fprintf(fid5,'%g\t',OF);%61
    fprintf(fid5,'%g\t',OS);%62
    fprintf(fid5,'%g\t',ZWT);%63
    %%%
    fprintf(fid5,'%g\t',DQ);%64
    fprintf(fid5,'%g\t',DT);%%65
    fprintf(fid5,'%g\t',dQ);%66
    %%%
    %%%
    fprintf(fid5,'%g\t',CK1);%67
    fprintf(fid5,'%g\t\n',CK2);%68
    if i==NN
        fclose(fid5);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%
    if i==2
        tit6=strcat('OUTPUT_',TITLE_SAVE,'_SOIL.dat');
        fid6=fopen(tit6,'a');
    end
    for kjk = 1:ms
        fprintf(fid6,'%g\t',V(1,kjk));
        fprintf(fid6,'%g\t',O(1,kjk));
        fprintf(fid6,'%g\t',Tdp(1,kjk));
    end
    fprintf(fid6,'%g\t\n',i);
    if i==NN
        fclose(fid6);
    end
    %%%%%%%%%%%%%%
    if i==2
        for ijki=1:cc
            tit7{1+ijki}=strcat('OUTPUT_',TITLE_SAVE,'_PFT_',num2str(ijki),'.dat');
            fid7(1+ijki)=fopen(tit7{1+ijki},'a');
        end
    end
    for ijki = 1:cc
        fprintf(fid7(1+ijki),'%g\t',T_H(ijki)); %1
        fprintf(fid7(1+ijki),'%g\t',T_L(ijki)); %2
        fprintf(fid7(1+ijki),'%g\t',EIn_H(ijki)); %3
        fprintf(fid7(1+ijki),'%g\t',EIn_L(ijki)); %4
        fprintf(fid7(1+ijki),'%g\t',Dr_H(ijki));  %5
        fprintf(fid7(1+ijki),'%g\t',Dr_L(ijki));  %6
        fprintf(fid7(1+ijki),'%g\t',In_H(ijki));  %7
        fprintf(fid7(1+ijki),'%g\t',In_L(ijki));  %8
        fprintf(fid7(1+ijki),'%g\t',Tdp_H(ijki));   %9
        fprintf(fid7(1+ijki),'%g\t',Tdp_L(ijki));  %10
        fprintf(fid7(1+ijki),'%g\t',OH(ijki));   %11
        fprintf(fid7(1+ijki),'%g\t',OL(ijki));  %12
        fprintf(fid7(1+ijki),'%g\t',rb_H(ijki)); %13
        fprintf(fid7(1+ijki),'%g\t',rb_L(ijki)); %14
        %%%%
        fprintf(fid7(1+ijki),'%g\t',rs_sunH(ijki)); %15
        fprintf(fid7(1+ijki),'%g\t',rs_sunL(ijki)); %16
        fprintf(fid7(1+ijki),'%g\t',rs_shdH(ijki)); %17
        fprintf(fid7(1+ijki),'%g\t',rs_shdL(ijki)); %18
        fprintf(fid7(1+ijki),'%g\t',rap_H(ijki)); %19
        fprintf(fid7(1+ijki),'%g\t',rap_L(ijki)); %20
        %%%
        fprintf(fid7(1+ijki),'%g\t',An_H(ijki)); %21
        fprintf(fid7(1+ijki),'%g\t',An_L(ijki)); %22
        fprintf(fid7(1+ijki),'%g\t',Rdark_H(ijki)); %23
        fprintf(fid7(1+ijki),'%g\t',Rdark_L(ijki)); %24
        %%%
        fprintf(fid7(1+ijki),'%g\t',Ci_sunH(ijki)); %25
        fprintf(fid7(1+ijki),'%g\t',Ci_sunL(ijki)); %26
        fprintf(fid7(1+ijki),'%g\t',Ci_shdH(ijki)); %27
        fprintf(fid7(1+ijki),'%g\t',Ci_shdL(ijki)); %28
        %%%
        fprintf(fid7(1+ijki),'%g\t',LAI_H(ijki)); %29
        fprintf(fid7(1+ijki),'%g\t',LAI_L(ijki)); %30
        fprintf(fid7(1+ijki),'%g\t',NPP_H(ijki)); %31
        fprintf(fid7(1+ijki),'%g\t',NPP_L(ijki)); %32
        fprintf(fid7(1+ijki),'%g\t',ANPP_H(ijki)); %33
        fprintf(fid7(1+ijki),'%g\t',ANPP_L(ijki)); %34
        fprintf(fid7(1+ijki),'%g\t',RA_H(ijki)); %35
        fprintf(fid7(1+ijki),'%g\t',RA_L(ijki)); %36
        fprintf(fid7(1+ijki),'%g\t',Rg_H(ijki)); %37
        fprintf(fid7(1+ijki),'%g\t',Rg_L(ijki)); %38
        fprintf(fid7(1+ijki),'%g\t',Rms_H(ijki)); %39
        fprintf(fid7(1+ijki),'%g\t',Rms_L(ijki)); %40
        fprintf(fid7(1+ijki),'%g\t',Rmr_H(ijki)); %41
        fprintf(fid7(1+ijki),'%g\t',Rmr_L(ijki)); %42
        fprintf(fid7(1+ijki),'%g\t',Rmc_H(ijki)); %43
        fprintf(fid7(1+ijki),'%g\t',Rmc_L(ijki)); %44
        fprintf(fid7(1+ijki),'%g\t',PHE_S_H(ijki)); %45
        fprintf(fid7(1+ijki),'%g\t',PHE_S_L(ijki)); %46
        fprintf(fid7(1+ijki),'%g\t',dflo_H(ijki)); %47
        fprintf(fid7(1+ijki),'%g\t',dflo_L(ijki)); %48
        fprintf(fid7(1+ijki),'%g\t',AgeL_H(ijki)); %49
        fprintf(fid7(1+ijki),'%g\t',AgeL_L(ijki)); %50
        fprintf(fid7(1+ijki),'%g\t',SAI_H(ijki));  %51
        fprintf(fid7(1+ijki),'%g\t',SAI_L(ijki));  %52
        fprintf(fid7(1+ijki),'%g\t',LAIdead_H(ijki));  %53
        fprintf(fid7(1+ijki),'%g\t',LAIdead_L(ijki));  %54
        fprintf(fid7(1+ijki),'%g\t',hc_H(ijki)); %55
        fprintf(fid7(1+ijki),'%g\t',hc_L(ijki)); %56
        fprintf(fid7(1+ijki),'%g\t',e_rel_L(ijki)); %57
        fprintf(fid7(1+ijki),'%g\t',e_rel_H(ijki)); %58
        fprintf(fid7(1+ijki),'%g\t',B_H(ijki,1)); %59
        fprintf(fid7(1+ijki),'%g\t',B_H(ijki,2)); %60
        fprintf(fid7(1+ijki),'%g\t',B_H(ijki,3)); %61
        fprintf(fid7(1+ijki),'%g\t',B_H(ijki,4)); %62
        fprintf(fid7(1+ijki),'%g\t',B_H(ijki,5)); %63
        fprintf(fid7(1+ijki),'%g\t',B_H(ijki,6)); %64
        fprintf(fid7(1+ijki),'%g\t',B_H(ijki,7)); %65
        fprintf(fid7(1+ijki),'%g\t',B_L(ijki,1)); %66
        fprintf(fid7(1+ijki),'%g\t',B_L(ijki,2)); %67
        fprintf(fid7(1+ijki),'%g\t',B_L(ijki,3)); %68
        fprintf(fid7(1+ijki),'%g\t',B_L(ijki,4)); %69
        fprintf(fid7(1+ijki),'%g\t',B_L(ijki,5)); %70
        fprintf(fid7(1+ijki),'%g\t',B_L(ijki,6)); %71
        fprintf(fid7(1+ijki),'%g\t',B_L(ijki,7)); %72
        fprintf(fid7(1+ijki),'%g\t\n',Llitter(ijki)); %73
    end
    if i==NN
        for ijki=1:cc
            fclose(fid7(1+ijki));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
Computational_Time =toc;
%profile off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('COMPUTATIONAL TIME [s] ')
disp(Computational_Time)
disp(' COMPUTATIONAL TIME [ms/cycle] ')
disp(1000*Computational_Time/NN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%







