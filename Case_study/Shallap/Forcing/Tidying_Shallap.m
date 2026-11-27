%Quick code to tidy Shallap meteo when using as test site

Meteo_data = load('SM_Data_cl_and_fl_20260202.mat','SM_fl4_hourly');
Meteo_data = Meteo_data.SM_fl4_hourly;
idForc = isbetween(Meteo_data.DateTime,Meteo_data.DateTime(1),"2024-09-30 23:00"); %This is to prevent Pr issues where NaNs
Meteo_data_short = Meteo_data(idForc,:);
To_keep_vars = {'DateTime','Ta','RH','u','Pr_CorRH','SWinCor','SWout','LWin','LWout'};
All_vars = Meteo_data_short.Properties.VariableNames;
To_keep = matches(All_vars,To_keep_vars);
Meteo_data_short = Meteo_data_short(:,To_keep);
Meteo_data_short = renamevars(Meteo_data_short,'u','WS');
Meteo_data_short = renamevars(Meteo_data_short,'Pr_CorRH','Pr');
Meteo_data_short = renamevars(Meteo_data_short,'SWinCor','SWin');

save('SM_Data_cl_and_fl_20260202_tidy.mat','Meteo_data_short');

%Also sort lapse rates
hm_Ta_lapse = load('SH_hm_Ta_lapse.mat');
hm_Ta_lapse = hm_Ta_lapse.SH_hm_lapse;
hm_Ta_lapse.Properties.VariableNames{5}='Ta_lapse'; 
save('SH_hm_Ta_lapse_tidy.mat',"hm_Ta_lapse");