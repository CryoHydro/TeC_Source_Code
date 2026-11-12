%==========================================================================
% OUTPUT CONFIGURATION - Kyzylsu
%==========================================================================
% User preferences for the distributed model outputs. This file is read by
% T&C_Code/Output_Manager.m (the same for all users) at the start of a run.
% It is a plain MATLAB script: edit it in any text editor, use % to comment
% out a line or switch off a variable.
%
% Units and descriptions are defined in Output_Manager.m (OM_catalogue) and
% are written next to the outputs (Outputs_info.csv / file attributes).
%==========================================================================

%% Global settings
%--------------------------------------------------------------------------
OUT_set.map_format   = {'mat','netcdf'};  % Maps: 'mat', 'netcdf' or both
OUT_set.point_format = {'csv','mat'};     % Points: 'csv', 'mat' or both
OUT_set.checkpoint_every_days = 1;       % Model days between two restart checkpoints
OUT_set.keep_n_checkpoints    = 2;        % Number of checkpoints kept in Store/ (older ones deleted)
OUT_set.single_precision      = true;     % Store outputs as single (half the size of double)

%% Variables
%--------------------------------------------------------------------------
% Name         : variable name in the model workspace
% Map period   : 'day', 'month', 'year' or '' (no map)
% Map func     : 'mean', 'sum', 'min', 'max' or 'last' (value at the end of the period)
% Point period : 'hour', 'day', 'month', 'year' or '' (no point output)
% Point func   : as Map func (ignored for 'hour')
% Crown        : only for variables with one column per crown (cc_max):
%                'weighted' = sum over crowns weighted by Ccrown (one value)
%                'each'     = one column per crown
%                ''         = not a crown variable
% Variables with one column per soil layer (O, Tdp, ...) are exported with
% one column per layer (Tdp_1 ... Tdp_<ms_max>), whatever ms_max is.
%
%  Name              Map period  Map func  Point period  Point func  Crown
OUT_var = {
 % --- Vegetation -------------------------------------------------------
  'ANPP_H'          'month'     'mean'    ''            ''          'weighted'
  'ANPP_L'          'month'     'mean'    ''            ''          'weighted'
  'LAI_H'           'month'     'mean'    'hour'        ''          'weighted'
  'LAI_L'           'month'     'mean'    'hour'        ''          'weighted'
  'NDVI'            'month'     'mean'    'hour'        ''          ''
  'NPP_H'           'month'     'mean'    'hour'        ''          'weighted'
  'NPP_L'           'month'     'mean'    'hour'        ''          'weighted'
  'hc_H'            'month'     'mean'    ''            ''          'weighted'
  'hc_L'            'month'     'mean'    ''            ''          'weighted'
  'An_H'            'month'     'mean'    'hour'        ''          'weighted'
  'An_L'            'month'     'mean'    'hour'        ''          'weighted'
  'RA_H'            'month'     'mean'    'hour'        ''          'weighted'
  'RA_L'            'month'     'mean'    'hour'        ''          'weighted'
  'Rg_H'            'month'     'mean'    'hour'        ''          'weighted'
  'Rg_L'            'month'     'mean'    'hour'        ''          'weighted'
  'SAI_H'           'month'     'mean'    'hour'        ''          'weighted'
  'SAI_L'           'month'     'mean'    'hour'        ''          'weighted'
  'PHE_S_H'         'month'     'last'    'hour'        ''          'weighted'
  'PHE_S_L'         'month'     'last'    'hour'        ''          'weighted'
 % --- Evaporation and transpiration --------------------------------------
  'EG'              'month'     'sum'     'hour'        ''          ''
  'EICE'            'month'     'sum'     'hour'        ''          ''
  'EIn_H'           'month'     'sum'     'hour'        ''          'weighted'
  'EIn_L'           'month'     'sum'     'hour'        ''          'weighted'
  'EIn_rock'        'month'     'sum'     'hour'        ''          ''
  'EIn_urb'         'month'     'sum'     ''            ''          ''
  'EWAT'            'month'     'sum'     ''            ''          ''
  'T_H'             'month'     'sum'     'hour'        ''          'weighted'
  'T_L'             'month'     'sum'     'hour'        ''          'weighted'
  'ESN'             'month'     'sum'     'hour'        ''          ''
  'ELitter'         'month'     'sum'     ''            ''          ''
  'ESN_In'          'month'     'sum'     ''            ''          ''
 % --- Snow and ice ---------------------------------------------------------
  'SND'             'day'       'mean'    'hour'        ''          ''   % daily map
  'SWE'             'month'     'mean'    'hour'        ''          ''
  'ros'             'month'     'mean'    'hour'        ''          ''
  'snow_albedo'     'month'     'mean'    'hour'        ''          ''
  'Smelt'           'month'     'sum'     'hour'        ''          ''
  'Imelt'           'month'     'sum'     'hour'        ''          ''
  'SP_wc'           'month'     'mean'    ''            ''          ''
  'SWE_avalanched'  'month'     'sum'     ''            ''          ''
  'Cicew'           'month'     'mean'    ''            ''          ''
  'Cice'            'month'     'mean'    ''            ''          ''
  'IP_wc'           'month'     'mean'    'hour'        ''          ''
  'ICE_D'           'month'     'mean'    'hour'        ''          ''
  'ICE'             'month'     'mean'    'hour'        ''          ''
  'U_SWE'           'month'     'sum'     ''            ''          ''
 % --- Groundwater -----------------------------------------------------------
  'FROCK'           'month'     'mean'    'hour'        ''          ''
  'Gfin'            'month'     'mean'    ''            ''          ''
  'Lk_rock'         'month'     'sum'     ''            ''          ''
  'Lk'              'month'     'sum'     ''            ''          ''
  'Lk_wat'          'month'     'sum'     ''            ''          ''
  'ZWT'             'month'     'mean'    'hour'        ''          ''
 % --- Soil moisture and soil temperature -----------------------------------
  'O'               'month'     'mean'    'hour'        ''          ''
  'OF'              'month'     'mean'    'hour'        ''          ''
  'OH'              'month'     'mean'    'hour'        ''          ''
  'OL'              'month'     'mean'    'hour'        ''          ''
  'OS'              'month'     'mean'    'hour'        ''          ''
  'Tdp'             'month'     'mean'    'hour'        ''          ''
  'Tdp_L'           'month'     'mean'    'hour'        ''          ''
  'Tdp_H'           'month'     'mean'    'hour'        ''          ''
 % --- Runoff and interception ----------------------------------------------
  % Routing outputs (QpointC, QpointH, UpointC, UpointH), computed at the outlet and the points:
  % Point period -> columns of the point files (Points/), Map period -> time series (Series/).
  % Catchment values (Q_exit, Qsub_exit): time series (Series/), period from the Map or Point columns.
  'QpointC'         ''          ''        'hour'        ''          ''   % discharge at the outlet + points
  'QpointH'         ''          ''        'hour'        ''          ''   % hillslope outflow at the outlet + points
  'UpointC'         ''          ''        'hour'        ''          ''   % channel flow velocity at the outlet + points
  'Q_exit'          ''          ''        'hour'        ''          ''   % surface water leaving the domain
  'Qsub_exit'       ''          ''        'hour'        ''          ''   % subsurface water leaving the domain
  'Q_channel'       ''          ''        'hour'        ''          ''   % channel water at the points (stream gauges): water level
  'q_runon'         'month'     'sum'     ''            ''          ''
  'Rd'              'month'     'sum'     ''            ''          ''
  'Rh'              'month'     'sum'     'hour'        ''          ''
  'f'               'month'     'sum'     ''            ''          ''
  'In_rock'         'month'     'mean'    'hour'        ''          ''
  'In_SWE'          'month'     'mean'    'hour'        ''          ''
  'In_L'            'month'     'mean'    'hour'        ''          'weighted'
  'In_H'            'month'     'mean'    'hour'        ''          'weighted'
 % --- Meteorology (forcing) -------------------------------------------------
  'Pr_liq'          'month'     'sum'     'hour'        ''          ''
  'Pr_sno'          'month'     'sum'     'hour'        ''          ''
  'Pr_St'           'month'     'sum'     'hour'        ''          ''
  'Ta_St'           'day'       'mean'    'hour'        ''          ''   % daily map
  'Tdew_St'         'month'     'mean'    'hour'        ''          ''
  'Ws_St'           'month'     'mean'    'hour'        ''          ''
  'U_St'            'month'     'mean'    'hour'        ''          ''
  'N_St'            'month'     'mean'    'hour'        ''          ''
  'ea_St'           'month'     'mean'    'hour'        ''          ''
  'Pre_St'          'month'     'mean'    'hour'        ''          ''
  'SAD1_St'         'month'     'mean'    'hour'        ''          ''
  'SAD2_St'         'month'     'mean'    'hour'        ''          ''
  'SAB1_St'         'month'     'mean'    'hour'        ''          ''
  'SAB2_St'         'month'     'mean'    'hour'        ''          ''
 % --- Energy balance --------------------------------------------------------
  'Ts'              'month'     'mean'    'hour'        ''          ''
  'Rn'              'month'     'mean'    'hour'        ''          ''
  'QE'              'month'     'mean'    'hour'        ''          ''
  'H'               'month'     'mean'    'hour'        ''          ''
  'HV'              'month'     'mean'    'hour'        ''          ''
  'G'               'month'     'mean'    'hour'        ''          ''
  'Qfm'             'month'     'mean'    'hour'        ''          ''
  'Qv'              'month'     'mean'    'hour'        ''          ''
  'dQ_S'            'month'     'mean'    'hour'        ''          ''
  'dQVEG'           'month'     'mean'    'hour'        ''          ''
  'Csno'            'month'     'mean'    'hour'        ''          ''
};
