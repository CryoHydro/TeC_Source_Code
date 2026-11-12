%==========================================================================
% OUTPUT CONFIGURATION - Kyzylsu
%==========================================================================
% User preferences for the distributed model outputs. This file is read by
% T&C_Code/Output_Manager.m (the same for all users) at the start of a run.
% It is a plain MATLAB script: edit it in any text editor, use % to comment
% out a line or switch off a variable.
%
% Units and descriptions are defined in Output_Manager.m (OM_catalogue) and
% are written next to the outputs (variables_info.csv / file attributes).
%==========================================================================

%% Global settings
%--------------------------------------------------------------------------
OUT_set.map_format   = {'mat','netcdf'};  % Maps: 'mat', 'netcdf' or both
OUT_set.point_format = {'csv','mat'};     % Points: 'csv', 'mat' or both
OUT_set.checkpoint_every_days = 30;       % Model days between two restart checkpoints
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
  'ANPP_H'          'day'       'mean'    ''            ''          'weighted'
  'ANPP_L'          'day'       'mean'    ''            ''          'weighted'
  'LAI_H'           'day'       'mean'    'hour'        ''          'weighted'
  'LAI_L'           'day'       'mean'    'hour'        ''          'weighted'
  'NDVI'            'day'       'mean'    'hour'        ''          ''
  'NPP_H'           'day'       'mean'    'hour'        ''          'weighted'
  'NPP_L'           'day'       'mean'    'hour'        ''          'weighted'
  'hc_H'            'day'       'mean'    ''            ''          'weighted'
  'hc_L'            'day'       'mean'    ''            ''          'weighted'
  'An_H'            'day'       'mean'    'hour'        ''          'weighted'
  'An_L'            'day'       'mean'    'hour'        ''          'weighted'
  'RA_H'            'day'       'mean'    'hour'        ''          'weighted'
  'RA_L'            'day'       'mean'    'hour'        ''          'weighted'
  'Rg_H'            'day'       'mean'    'hour'        ''          'weighted'
  'Rg_L'            'day'       'mean'    'hour'        ''          'weighted'
  'SAI_H'           'day'       'mean'    'hour'        ''          'weighted'
  'SAI_L'           'day'       'mean'    'hour'        ''          'weighted'
  'PHE_S_H'         'day'       'last'    'hour'        ''          'weighted'
  'PHE_S_L'         'day'       'last'    'hour'        ''          'weighted'
 % --- Evaporation and transpiration --------------------------------------
  'EG'              'day'       'sum'     'hour'        ''          ''
  'EICE'            'day'       'sum'     'hour'        ''          ''
  'EIn_H'           'day'       'sum'     'hour'        ''          'weighted'
  'EIn_L'           'day'       'sum'     'hour'        ''          'weighted'
  'EIn_rock'        'day'       'sum'     'hour'        ''          ''
  'EIn_urb'         'day'       'sum'     ''            ''          ''
  'EWAT'            'day'       'sum'     ''            ''          ''
  'T_H'             'day'       'sum'     'hour'        ''          'weighted'
  'T_L'             'day'       'sum'     'hour'        ''          'weighted'
  'ESN'             'day'       'sum'     'hour'        ''          ''
  'ELitter'         'day'       'sum'     ''            ''          ''
  'ESN_In'          'day'       'sum'     ''            ''          ''
 % --- Snow and ice ---------------------------------------------------------
  'SND'             'day'       'mean'    'hour'        ''          ''
  'SWE'             'day'       'mean'    'hour'        ''          ''
  'ros'             'day'       'mean'    'hour'        ''          ''
  'snow_albedo'     'day'       'mean'    'hour'        ''          ''
  'Smelt'           'day'       'sum'     'hour'        ''          ''
  'Imelt'           'day'       'sum'     'hour'        ''          ''
  'SP_wc'           'day'       'mean'    ''            ''          ''
  'SWE_avalanched'  'day'       'sum'     ''            ''          ''
  'Cicew'           'day'       'mean'    ''            ''          ''
  'Cice'            'day'       'mean'    ''            ''          ''
  'IP_wc'           'day'       'mean'    'hour'        ''          ''
  'ICE_D'           'day'       'mean'    'hour'        ''          ''
  'ICE'             'day'       'mean'    'hour'        ''          ''
  'U_SWE'           'day'       'sum'     ''            ''          ''
 % --- Groundwater -----------------------------------------------------------
  'FROCK'           'day'       'mean'    'hour'        ''          ''
  'Gfin'            'day'       'mean'    ''            ''          ''
  'Lk_rock'         'day'       'sum'     ''            ''          ''
  'Lk'              'day'       'sum'     ''            ''          ''
  'Lk_wat'          'day'       'sum'     ''            ''          ''
  'ZWT'             'day'       'mean'    'hour'        ''          ''
 % --- Soil moisture and soil temperature -----------------------------------
  'O'               'day'       'mean'    'hour'        ''          ''
  'OF'              'day'       'mean'    'hour'        ''          ''
  'OH'              'day'       'mean'    'hour'        ''          ''
  'OL'              'day'       'mean'    'hour'        ''          ''
  'OS'              'day'       'mean'    'hour'        ''          ''
  'Tdp'             'day'       'mean'    'hour'        ''          ''
  'Tdp_L'           'day'       'mean'    'hour'        ''          ''
  'Tdp_H'           'day'       'mean'    'hour'        ''          ''
 % --- Runoff and interception ----------------------------------------------
  'QpointC'         'day'       'mean'    ''            ''          ''   % outlets: written as time series, not maps
  'q_runon'         'day'       'sum'     ''            ''          ''
  'Rd'              'day'       'sum'     ''            ''          ''
  'Rh'              'day'       'sum'     'hour'        ''          ''
  'f'               'day'       'sum'     ''            ''          ''
  'In_rock'         'day'       'mean'    'hour'        ''          ''
  'In_SWE'          'day'       'mean'    'hour'        ''          ''
  'In_L'            'day'       'mean'    'hour'        ''          'weighted'
  'In_H'            'day'       'mean'    'hour'        ''          'weighted'
 % --- Meteorology (forcing) -------------------------------------------------
  'Pr_liq'          'day'       'sum'     'hour'        ''          ''
  'Pr_sno'          'day'       'sum'     'hour'        ''          ''
  'Pr_St'           'day'       'sum'     'hour'        ''          ''
  'Ta_St'           'day'       'mean'    'hour'        ''          ''
  'Tdew_St'         'day'       'mean'    'hour'        ''          ''
  'Ws_St'           'day'       'mean'    'hour'        ''          ''
  'U_St'            'day'       'mean'    'hour'        ''          ''
  'N_St'            'day'       'mean'    'hour'        ''          ''
  'ea_St'           'day'       'mean'    'hour'        ''          ''
  'Pre_St'          'day'       'mean'    'hour'        ''          ''
  'SAD1_St'         'day'       'mean'    'hour'        ''          ''
  'SAD2_St'         'day'       'mean'    'hour'        ''          ''
  'SAB1_St'         'day'       'mean'    'hour'        ''          ''
  'SAB2_St'         'day'       'mean'    'hour'        ''          ''
 % --- Energy balance --------------------------------------------------------
  'Ts'              'day'       'mean'    'hour'        ''          ''
  'Rn'              'day'       'mean'    'hour'        ''          ''
  'QE'              'day'       'mean'    'hour'        ''          ''
  'H'               'day'       'mean'    'hour'        ''          ''
  'HV'              'day'       'mean'    'hour'        ''          ''
  'G'               'day'       'mean'    'hour'        ''          ''
  'Qfm'             'day'       'mean'    'hour'        ''          ''
  'Qv'              'day'       'mean'    'hour'        ''          ''
  'dQ_S'            'day'       'mean'    'hour'        ''          ''
  'dQVEG'           'day'       'mean'    'hour'        ''          ''
  'Csno'            'day'       'mean'    'hour'        ''          ''
};
