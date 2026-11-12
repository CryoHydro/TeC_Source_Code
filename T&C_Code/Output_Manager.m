%==========================================================================
% OUTPUT MANAGER - T&C distributed model
%==========================================================================
% Same script for all users. User preferences (which variables, how often,
% which aggregation, file formats, checkpoints) are read from:
%       <site_name>/Parameters/Output_Config.m
%
% Called from the launcher at three points:
%   OM_phase = 'init';  Output_Manager   % just before the time loop. Sets
%                                        % t_start (= fts, or the step after
%                                        % the checkpoint when restarting)
%   OM_phase = 'step';  Output_Manager   % last line inside the time loop
%   OM_phase = 'close'; Output_Manager   % after the time loop
%
% What it does
%   - Aggregates the chosen variables on the fly (running sum/mean/min/max/
%     last), so no hourly cubes are kept in memory.
%   - Writes every record as soon as its period is complete (maps: daily,
%     monthly, yearly; points: hourly, daily, monthly, yearly). All writing
%     happens at the end of each model day, so outputs can be inspected while
%     the run is going.
%   - Maps:   Maps/<period>/<Var>/<Var>_<period>_<yyyy>.mat / .nc  (one folder
%             per variable, one file per calendar year, records indexed by date)
%   - Points: Points/<period>/<Point>_<period>.csv            (whole run)
%             Points/<period>/<Point>_<period>_<yyyy>.mat     (per year)
%   - Series: Series/<period>/<Var>_<period>.csv / _<yyyy>.mat (outlets, e.g. QpointC)
%   - Outputs_info.csv: units, description and aggregation of every column.
%   - run_status.txt:   progress of the run (updated every model day).
%   - Store/Checkpoint_<date>.mat: restart checkpoints.
%
% Restart
%   A checkpoint is the whole model workspace (state variables, memory of
%   the previous 24 h, the output manager state OM, the time step t) except
%   the forcing, the run settings and the paths, which are re-read by the
%   launcher. With Restart_ID = 1 the launcher sets everything up as usual,
%   then this script loads the latest checkpoint (or the one of Restart_Date),
%   removes from the output files any record written after the checkpoint,
%   and the time loop continues from the next step. The model parameters
%   come from the checkpoint, so the continued run is identical to an
%   uninterrupted one. Run_EndDate can be changed to extend a run.
%
% All manager state is in the structure OM. Temporary variables start with
% OM_ and are cleared at the end of each call.
%==========================================================================

switch OM_phase

    %======================================================================
    case 'init'
    %======================================================================
        OM_cfg = [site_name '/Parameters/Output_Config.m'];
        if ~isfile(OM_cfg)
            error('Output_Manager: configuration file not found: %s', OM_cfg);
        end
        clear OUT_set OUT_var
        run(OM_cfg); % defines OUT_set and OUT_var

        OM_grid.mask = (MASK == 1);
        OM_grid.nrow = size(MASK,1);
        OM_grid.ncol = size(MASK,2);
        OM_grid.x = 1:OM_grid.ncol; OM_grid.y = (1:OM_grid.nrow)'; OM_grid.xy_units = 'cell index';
        if exist('x','var') && exist('y','var') && numel(x) == OM_grid.ncol && numel(y) == OM_grid.nrow
            OM_grid.x = x(:)'; OM_grid.y = y(:); OM_grid.xy_units = 'm (UTM)';
        end
        OM_new = om_setup(OUT_set, OUT_var, Directories.save, OM_grid, Date(1), site_name);
        clear OUT_set OUT_var

        if restart.id == 1
            OM_file = om_find_checkpoint(OM_new.dir.store, restart.date);
            fprintf('Output_Manager: restarting from %s\n', OM_file);
            OM_load = setdiff(who('-file', OM_file), OM_new.ckpt_exclude); % never restore forcing, run settings or run-length series
            load(OM_file, OM_load{:}); % model state, t and OM
            if OM.run_start ~= Date(1)
                error('Output_Manager: the checkpoint belongs to a run starting %s, this run starts %s.', ...
                    char(OM.run_start), char(Date(1)));
            end
            if ~isequal(OM.cfg_sig, OM_new.cfg_sig)
                warning('Output_Manager: Output_Config.m changed since the checkpoint. The configuration of the checkpoint is kept for this run.');
            end
            OM.dir = OM_new.dir;
            om_trim_outputs(OM);
            t_start = t + 1;
            fprintf('Output_Manager: outputs complete until %s, continuing from step %d\n', char(OM.tc), t_start);
        else
            OM = OM_new;
            t_start = fts;
        end
        OM.job_clock = tic;
        OM.job_t_start = t_start;
        OM.N_time_step = N_time_step;

    %======================================================================
    case 'step'
    %======================================================================
        OM_date = datetime(Datam(t-1,1), Datam(t-1,2), Datam(t-1,3), Datam(t-1,4), 0, 0); % date of this step

        % Values of the output variables at this step
        OM_vals = cell(OM.nvar,1);
        OM_missing = false(OM.nvar,1);
        for OM_i = 1:OM.nvar
            if exist(OM.var(OM_i).name, 'var')
                OM_vals{OM_i} = eval(OM.var(OM_i).name);
                if ismatrix(OM_vals{OM_i}) && ~isvector(OM_vals{OM_i}) && numel(OM_vals{OM_i}) == num_cell
                    OM_vals{OM_i} = OM_vals{OM_i}(:); % 2-D map (m_cell x n_cell, e.g. Q_channel): one value per cell
                end
            else
                OM_missing(OM_i) = true;
            end
        end
        if any(OM_missing)
            error('Output_Manager: these variables of Output_Config.m do not exist in the model: %s', ...
                strjoin({OM.var(OM_missing).name}, ', '));
        end

        if ~OM.ready % first step: sizes of the variables
            % Cells of the routing outputs (QpointC, QpointH, UpointC, UpointH): Xout, Yout
            if exist('Xout','var') && exist('Yout','var')
                OM_idx_routed = sub2ind([m_cell n_cell], Yout(:), Xout(:));
            else
                OM_idx_routed = [];
            end
            OM = om_build(OM, OM_vals, num_cell, cc_max, Points, NAMEout, OM_idx_routed);
        end

        [OM, OM_day_end] = om_step(OM, OM_vals, Ccrown_OUTv, OM_date);

        OM_last = (t == N_time_step);
        if OM_day_end || OM_last
            OM.days_since_ckpt = OM.days_since_ckpt + OM_day_end;
            if OM.days_since_ckpt >= OM.set.checkpoint_every_days || OM_last
                % Save the checkpoint (before closing incomplete periods at the end of the run)
                OM.tc = OM_date + hours(1);   % outputs are complete up to this time
                OM.days_since_ckpt = 0;
                OM_file = [OM.dir.store 'Checkpoint_' char(datetime(OM.tc,'Format','yyyy-MM-dd_HHmm')) '.mat'];
                OM.last_ckpt = OM_file;
                OM_list = who;
                OM_list = OM_list(~ismember(OM_list, OM.ckpt_exclude) & ~startsWith(OM_list, 'OM_'));
                OM_tclock = tic;
                save([OM_file '.tmp'], OM_list{:}, '-v7.3');
                movefile([OM_file '.tmp'], OM_file, 'f');
                om_clean_checkpoints(OM.dir.store, OM.set.keep_n_checkpoints);
                fprintf('Output_Manager: checkpoint saved (%s, %.1f s)\n', OM_file, toc(OM_tclock));
            end
            if OM_last
                OM = om_close_all(OM); % write incomplete periods (end of the run)
            end
            om_write_status(OM, OM_date, t, false);
        end

    %======================================================================
    case 'close'
    %======================================================================
        om_write_status(OM, datetime(Datam(t-1,1),Datam(t-1,2),Datam(t-1,3),Datam(t-1,4),0,0), t, true);
        fprintf('Output_Manager: run finished, outputs in %s\n', OM.dir.root);

    otherwise
        error('Output_Manager: unknown OM_phase "%s" (init, step or close)', OM_phase);
end

clear -regexp ^OM_

%##########################################################################
%                           LOCAL FUNCTIONS
%##########################################################################

%--------------------------------------------------------------------------
function OM = om_setup(S, V, save_dir, grid, run_start, site_name)
% Reads and checks the user configuration
if ~iscell(V) || size(V,2) ~= 6
    error('Output_Manager: OUT_var must be a cell array with 6 columns (Name, Map period, Map func, Point period, Point func, Crown).');
end
def = struct('map_format',{{'mat'}}, 'point_format',{{'csv'}}, 'checkpoint_every_days',30, ...
             'keep_n_checkpoints',2, 'single_precision',true);
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(S, fn{i}); S.(fn{i}) = def.(fn{i}); end
end
S.map_format = cellstr(lower(S.map_format)); S.point_format = cellstr(lower(S.point_format));
if any(~ismember(S.map_format, {'mat','netcdf'})); error('Output_Manager: map_format must be ''mat'' and/or ''netcdf''.'); end
if any(~ismember(S.point_format, {'mat','csv'})); error('Output_Manager: point_format must be ''mat'' and/or ''csv''.'); end

V(cellfun(@isempty, V)) = {''};
names = V(:,1);
[u, ~, j] = unique(names);
if numel(u) < numel(names)
    error('Output_Manager: variables listed twice in Output_Config.m: %s', strjoin(u(accumarray(j,1) > 1), ', '));
end
funcs = {'mean','sum','min','max','last'};
err = {};
for i = 1:size(V,1)
    nm = V{i,1}; mp = lower(V{i,2}); mf = lower(V{i,3}); pp = lower(V{i,4}); pf = lower(V{i,5}); cr = lower(V{i,6});
    if ~isvarname(nm); err{end+1} = sprintf('"%s" is not a valid variable name', nm); end %#ok<*AGROW>
    if strcmp(mp,'hour'); err{end+1} = sprintf('%s: hourly maps are not allowed (daily is the highest map frequency)', nm); end
    if ~ismember(mp, {'','day','month','year'}); err{end+1} = sprintf('%s: unknown map period "%s"', nm, mp); end
    if ~isempty(mp) && ~ismember(mf, funcs); err{end+1} = sprintf('%s: unknown map func "%s"', nm, mf); end
    if ~ismember(pp, {'','hour','day','month','year'}); err{end+1} = sprintf('%s: unknown point period "%s"', nm, pp); end
    if ismember(pp, {'day','month','year'}) && ~ismember(pf, funcs); err{end+1} = sprintf('%s: unknown point func "%s"', nm, pf); end
    if ~ismember(cr, {'','weighted','each'}); err{end+1} = sprintf('%s: Crown must be '''', ''weighted'' or ''each''', nm); end
    if strcmp(pp,'hour'); pf = 'last'; end
    V(i,2:6) = {mp, mf, pp, pf, cr};
end
if ~isempty(err)
    error('Output_Manager: errors in Output_Config.m:\n  %s', strjoin(err, '\n  '));
end
V = V(~(cellfun(@isempty,V(:,2)) & cellfun(@isempty,V(:,4))), :); % variables not exported

cat = om_catalogue();
OM.var = struct('name',V(:,1), 'map_period',V(:,2), 'map_func',V(:,3), 'pt_period',V(:,4), ...
                'pt_func',V(:,5), 'crown',V(:,6), 'units','', 'desc','', 'kind','', 'ncol_in',0, 'cols',{{}}, ...
                'routed',false, 'series_cols',{{}});
for i = 1:numel(OM.var)
    k = find(strcmp(cat(:,1), OM.var(i).name), 1);
    if isempty(k)
        OM.var(i).units = '?'; OM.var(i).desc = '(not in the Output_Manager catalogue)';
    else
        OM.var(i).units = cat{k,2}; OM.var(i).desc = cat{k,3};
    end
end
OM.nvar = numel(OM.var);
OM.set = S;
OM.cfg_sig = {S, V};
OM.grid = grid;
OM.site = site_name;
OM.run_start = run_start;
OM.ready = false;
OM.days_since_ckpt = 0;
OM.tc = NaT;
OM.last_ckpt = '';
OM.class = 'single'; if ~S.single_precision; OM.class = 'double'; end

OM.dir.root   = save_dir;
OM.dir.maps   = [save_dir 'Maps/'];
OM.dir.points = [save_dir 'Points/'];
OM.dir.series = [save_dir 'Series/'];
OM.dir.store  = [save_dir 'Store/'];
OM.dir.status = [save_dir 'run_status.txt'];

% Variables not saved in the checkpoints: forcing (re-loaded each month),
% run settings and paths (re-read by the launcher at restart)
OM.ckpt_exclude = {'Ta_S','Pr_S','N_S','Tdew_S','SAB1_S','SAB2_S','PARB_S','Pr_ratio','Pre_S','U_S', ...
    'Ws_S','SAD1_S','SAD2_S','PARD_S','esat_S','ea_S','Ds_S','Forc_S','Forc_S_dc','Forc_S_dc_r', ...
    'idcli_r','idTa','idapply','Meteo_data','Ta','U','Ws','Pr','SWin','SWout','Nin','LWout','Uf', ...
    'esat','ea','Ds','Tdew','SD','SB','SAD1','SAD2','SAB1','SAB2','PARB','PARD','dsERA','dsGRID', ...
    'dsCoeff','demCurv','hm_Ta_lapse','Pr_lapse','Forcing_Date','Forcing_Date_num','year_loaded', ...
    'month_loaded','date_forMonth','id_date_for_month','restart','Directories','CONFIG','CONFIG_vals', ...
    'dateRun','x1','x2','Date','N_time_step','Nd_time_step','Datam','outlocation','Func_Dir', ...
    'OUT_set','OUT_var','ans','Ca'}; % Ca: CO2 series cut to the run period
end

%--------------------------------------------------------------------------
function OM = om_build(OM, vals, num_cell, cc_max, Points, NAMEout, idx_routed)
% Sizes and column names of the variables (done at the first step)
% idx_routed: cells (linear index) of the routing outputs, one per value of
% QpointC, QpointH, UpointC, UpointH ([] if unknown)
OM.idx_routed = idx_routed(:);
OM.num_cell = num_cell;
missing_pts = ~ismember(Points.idx(:), OM.idx_routed);
for i = 1:OM.nvar
    v = OM.var(i); X = vals{i};
    if isrow(X) && numel(X) == num_cell; X = X'; end
    v.routed = false;
    if numel(X) > 1 && numel(X) == numel(OM.idx_routed) && size(X,1) ~= num_cell && ~isempty(v.pt_period)
        % Routing output with a Point period: written in the point files
        % (one value per point, at its cell); also as a series if it has a Map period
        if any(missing_pts)
            error('Output_Manager: %s is computed only at the routed cells (Xout, Yout), not at the points %s.', ...
                v.name, strjoin(cellstr(string(Points.Name(missing_pts))), ', '));
        end
        v.kind = 'grid'; v.routed = true; v.ncol_in = 1; v.cols = {v.name};
        if ~isempty(v.map_period)   % series of the routed cells, as before
            onames = cellstr(string(NAMEout));
            if numel(onames) == numel(X); v.series_cols = matlab.lang.makeValidName(onames(:)');
            else; v.series_cols = compose([v.name '_%d'], 1:numel(X)); end
        end
        OM.var(i) = v;
        continue
    end
    if size(X,1) == num_cell
        v.kind = 'grid'; v.ncol_in = size(X,2);
        if ~isempty(v.crown)
            if v.ncol_in ~= cc_max
                error('Output_Manager: %s has %d columns, not cc_max = %d, so Crown must be ''''.', v.name, v.ncol_in, cc_max);
            end
            if strcmp(v.crown,'weighted'); v.cols = {v.name};
            else; v.cols = compose([v.name '_c%d'], 1:cc_max); end
        elseif v.ncol_in == 1
            v.cols = {v.name};
        else
            v.cols = compose([v.name '_%d'], 1:v.ncol_in); % e.g. soil layers 1..ms_max
        end
    else % outlet variables (e.g. QpointC): written as time series
        v.kind = 'series'; v.ncol_in = numel(X);
        onames = cellstr(string(NAMEout)); % NAMEout can be a char (one outlet) or a list
        if numel(X) == 1 && numel(onames) ~= 1; v.cols = {v.name};     % catchment value (e.g. Q_exit)
        elseif numel(onames) == numel(X); v.cols = matlab.lang.makeValidName(onames(:)');
        else; v.cols = compose([v.name '_%d'], 1:numel(X)); end
        % Period of the time series: Map period/func, or Point period/func if no
        % map period is given (needed for hourly series: maps are at most daily)
        if isempty(v.map_period)
            v.map_period = v.pt_period; v.map_func = v.pt_func;
        elseif ~isempty(v.pt_period)
            warning('Output_Manager: %s is not a gridded variable: its point settings are ignored, its map settings define its time series.', v.name);
        end
    end
    OM.var(i) = v;
end

per = {'hour','day','month','year'};
% Map accumulators (one per variable)
im = find(strcmp({OM.var.kind},'grid') & ~cellfun(@isempty,{OM.var.map_period}) & ~[OM.var.routed]);
OM.map = struct('vi',num2cell(im), 'period',{OM.var(im).map_period}, 'func',{OM.var(im).map_func}, ...
                'A',[], 'n',0, 't0',NaT);
% Table groups: points (one group per period) and series (one per variable)
OM.tab = struct('type',{}, 'period',{}, 'vi',{}, 'funcs',{}, 'files',{}, 'rows',{}, 'cols',{}, ...
                'A',{}, 'n',{}, 't0',{}, 'buf',{});
pnames = matlab.lang.makeValidName(cellstr(string(Points.Name)));
for p = 1:numel(per)
    iv = find(strcmp({OM.var.kind},'grid') & strcmp({OM.var.pt_period}, per{p}));
    if isempty(iv); continue; end
    g = struct('type','points', 'period',per{p}, 'vi',iv, 'funcs',{{OM.var(iv).pt_func}}, ...
               'files',{pnames(:)'}, 'rows',Points.idx(:), 'cols',{[OM.var(iv).cols]}, ...
               'A',{cell(1,numel(iv))}, 'n',0, 't0',NaT, 'buf',struct('t0',{},'n',{},'V',{}));
    OM.tab(end+1) = g;
end
for i = find(strcmp({OM.var.kind},'series') & ~cellfun(@isempty,{OM.var.map_period}))
    g = struct('type','series', 'period',OM.var(i).map_period, 'vi',i, 'funcs',{{OM.var(i).map_func}}, ...
               'files',{{OM.var(i).name}}, 'rows',1, 'cols',{OM.var(i).cols}, ...
               'A',{cell(1,1)}, 'n',0, 't0',NaT, 'buf',struct('t0',{},'n',{},'V',{}));
    OM.tab(end+1) = g;
end
for i = find([OM.var.routed] & ~cellfun(@isempty,{OM.var.map_period}))   % routing outputs also as series
    g = struct('type','series', 'period',OM.var(i).map_period, 'vi',i, 'funcs',{{OM.var(i).map_func}}, ...
               'files',{{OM.var(i).name}}, 'rows',OM.idx_routed, 'cols',{OM.var(i).series_cols}, ...
               'A',{cell(1,1)}, 'n',0, 't0',NaT, 'buf',struct('t0',{},'n',{},'V',{}));
    OM.tab(end+1) = g;
end

% Folders and overview of all outputs
dirs = {OM.dir.root, OM.dir.store};
for k = 1:numel(OM.map); dirs{end+1} = [OM.dir.maps OM.map(k).period '/' OM.var(OM.map(k).vi).name]; end
for k = 1:numel(OM.tab)
    if strcmp(OM.tab(k).type,'points'); dirs{end+1} = [OM.dir.points OM.tab(k).period];
    else; dirs{end+1} = [OM.dir.series OM.tab(k).period]; end
end
dirs = unique(dirs);
for k = 1:numel(dirs); if ~isfolder(dirs{k}); mkdir(dirs{k}); end; end
om_write_info(OM);
OM.ready = true;
fprintf('Output_Manager: %d map variables, %d table outputs (%s, %s)\n', numel(OM.map), numel(OM.tab), ...
    strjoin(OM.set.map_format,'+'), strjoin(OM.set.point_format,'+'));
end

%--------------------------------------------------------------------------
function [OM, day_end] = om_step(OM, vals, Ccrown, d)
% Accumulates this step and writes the periods that end with this step
ends.hour  = true;
ends.day   = hour(d) == 23;
ends.month = ends.day && day(d) == eomday(year(d), month(d));
ends.year  = ends.month && month(d) == 12;
day_end = ends.day;

X = cell(OM.nvar,1);
for i = 1:OM.nvar
    if isfield(OM.var, 'routed') && OM.var(i).routed
        % Routing output (one value per routed cell) -> value at each routed cell
        % of the grid, so that the point and series tables can pick their cells
        Xg = nan(OM.num_cell, 1); Xg(OM.idx_routed) = double(vals{i}(:));
        X{i} = Xg;
    else
        X{i} = om_prep(OM.var(i), vals{i}, Ccrown);
    end
end

for k = 1:numel(OM.map)
    m = OM.map(k);
    if m.n == 0; m.t0 = dateshift(d, 'start', m.period); end
    m.A = om_acc(m.A, X{m.vi}, m.func); m.n = m.n + 1;
    if ends.(m.period)
        om_write_map(OM, m.vi, m.period, m.t0, m.n, om_final(m.A, m.n, m.func));
        m.A = []; m.n = 0;
    end
    OM.map(k) = m;
end

for k = 1:numel(OM.tab)
    g = OM.tab(k);
    if g.n == 0; g.t0 = dateshift(d, 'start', g.period); end
    for b = 1:numel(g.vi)
        Xb = X{g.vi(b)}(g.rows,:);
        if strcmp(g.type,'series') && size(Xb,1) > 1; Xb = Xb'; end % routed series: one row, one column per cell
        g.A{b} = om_acc(g.A{b}, Xb, g.funcs{b});
    end
    g.n = g.n + 1;
    if ends.(g.period); g = om_tab_close(g); end
    if day_end; om_write_tab(OM, g); g.buf = g.buf([]); end
    OM.tab(k) = g;
end
end

%--------------------------------------------------------------------------
function OM = om_close_all(OM)
% End of the run: write the incomplete periods
for k = 1:numel(OM.map)
    m = OM.map(k);
    if m.n > 0
        om_write_map(OM, m.vi, m.period, m.t0, m.n, om_final(m.A, m.n, m.func));
        m.A = []; m.n = 0; OM.map(k) = m;
    end
end
for k = 1:numel(OM.tab)
    g = OM.tab(k);
    if g.n > 0; g = om_tab_close(g); end
    om_write_tab(OM, g); g.buf = g.buf([]);
    OM.tab(k) = g;
end
end

%--------------------------------------------------------------------------
function g = om_tab_close(g)
% Moves the aggregated values of a table group to its write buffer
V = cell(1, numel(g.vi));
for b = 1:numel(g.vi); V{b} = om_final(g.A{b}, g.n, g.funcs{b}); end
g.buf(end+1) = struct('t0', g.t0, 'n', g.n, 'V', [V{:}]);
g.A = cell(1, numel(g.vi)); g.n = 0;
end

%--------------------------------------------------------------------------
function X = om_prep(v, X, Ccrown)
% Values of one variable as (cells x columns), crowns weighted if asked
X = double(X);
if strcmp(v.kind, 'series'); X = X(:)'; return; end
if isrow(X) && v.ncol_in == 1; X = X'; end
if strcmp(v.crown, 'weighted')
    X(Ccrown == 0) = 0;
    X = sum(X .* Ccrown, 2);
end
end

%--------------------------------------------------------------------------
function A = om_acc(A, X, func)
if isempty(A)
    A = X;
    return
end
switch func
    case {'mean','sum'}; A = A + X;
    case 'min';  A = min(A, X);
    case 'max';  A = max(A, X);
    case 'last'; A = X;
end
end

%--------------------------------------------------------------------------
function Y = om_final(A, n, func)
if strcmp(func, 'mean'); Y = A / n; else; Y = A; end
end

%--------------------------------------------------------------------------
function [nrec, k] = om_rec_index(t0, period)
% Number of records in the yearly file and index of the record of t0
yr = year(t0);
switch period
    case 'hour';  nrec = 24 * (365 + om_isleap(yr)); k = (day(t0,'dayofyear') - 1) * 24 + hour(t0) + 1;
    case 'day';   nrec = 365 + om_isleap(yr);        k = day(t0, 'dayofyear');
    case 'month'; nrec = 12;                          k = month(t0);
    case 'year';  nrec = 1;                           k = 1;
end
end

function b = om_isleap(yr)
b = (mod(yr,4) == 0 & mod(yr,100) ~= 0) | mod(yr,400) == 0;
end

function te = om_period_end(t0, period)
switch period
    case 'hour';  te = t0 + hours(1);
    case 'day';   te = t0 + days(1);
    case 'month'; te = t0 + calmonths(1);
    case 'year';  te = t0 + calyears(1);
end
end

%--------------------------------------------------------------------------
function u = om_units_out(u, func, period)
if strcmp(func, 'sum')
    if endsWith(u, '/h'); u = [u(1:end-2) '/' period];
    elseif ~strcmp(period, 'hour'); u = [u ' (summed over hours)']; end
end
end

function cm = om_cell_method(func)
switch func
    case 'mean'; cm = 'time: mean';
    case 'sum';  cm = 'time: sum';
    case 'min';  cm = 'time: minimum';
    case 'max';  cm = 'time: maximum';
    case 'last'; cm = 'time: point (value at the end of the period)';
end
end

%--------------------------------------------------------------------------
function om_write_map(OM, vi, period, t0, n, Y)
% Writes one map record (all formats)
v = OM.var(vi); g = OM.grid; nc = size(Y, 2);
G = NaN(g.nrow, g.ncol, nc, OM.class);
for c = 1:nc
    tmp = reshape(Y(:,c), g.nrow, g.ncol);
    tmp(~g.mask) = NaN;
    G(:,:,c) = tmp;
end
[nrec, k] = om_rec_index(t0, period);
base = sprintf('%s%s/%s/%s_%s_%d', OM.dir.maps, period, v.name, v.name, period, year(t0));
units = om_units_out(v.units, v.map_func, period);

if ismember('mat', OM.set.map_format)
    f = [base '.mat'];
    if ~isfile(f)
        S.(v.name) = NaN(g.nrow, g.ncol, nc, nrec, OM.class);
        if nc == 1; S.(v.name) = reshape(S.(v.name), g.nrow, g.ncol, nrec); end
        S.time = NaN(nrec, 1);   % start of each period (datenum), NaN = not written
        S.n_steps = zeros(nrec, 1);
        S.info = struct('variable',v.name, 'units',units, 'description',v.desc, 'period',period, ...
            'func',v.map_func, 'columns',{v.cols}, 'x',g.x, 'y',g.y, 'xy_units',g.xy_units, 'site',OM.site, ...
            'time_note','time = start of the period (datenum): datetime(time,''ConvertFrom'',''datenum'')', ...
            'dims', 'rows x cols x [columns] x time');
        save(f, '-struct', 'S', '-v7.3');
    end
    m = matfile(f, 'Writable', true);
    if nrec == 1; m.(v.name) = G;   % yearly file: a single record
    elseif nc == 1; m.(v.name)(:,:,k) = G;
    else; m.(v.name)(:,:,:,k) = G; end
    m.time(k,1) = datenum(t0);
    m.n_steps(k,1) = n;
end

if ismember('netcdf', OM.set.map_format)
    f = [base '.nc'];
    if ~isfile(f)
        if nc == 1; dims = {'x', g.ncol, 'y', g.nrow, 'time', nrec};
        else; dims = {'x', g.ncol, 'y', g.nrow, 'column', nc, 'time', nrec}; end
        nccreate(f, v.name, 'Dimensions', dims, 'Datatype', OM.class, 'FillValue', NaN, ...
                 'Format', 'netcdf4', 'DeflateLevel', 4);
        nccreate(f, 'time', 'Dimensions', {'time', nrec}, 'Datatype', 'double', 'FillValue', NaN);
        nccreate(f, 'n_steps', 'Dimensions', {'time', nrec}, 'Datatype', 'int32', 'FillValue', int32(0));
        nccreate(f, 'x', 'Dimensions', {'x', g.ncol}, 'Datatype', 'double');
        nccreate(f, 'y', 'Dimensions', {'y', g.nrow}, 'Datatype', 'double');
        ncwrite(f, 'x', double(g.x(:))); ncwrite(f, 'y', double(g.y(:)));
        ncwriteatt(f, 'x', 'units', g.xy_units); ncwriteatt(f, 'y', 'units', g.xy_units);
        ncwriteatt(f, 'time', 'units', 'days since 1970-01-01 00:00:00');
        ncwriteatt(f, 'time', 'long_name', 'start of the aggregation period');
        ncwriteatt(f, 'n_steps', 'long_name', 'number of model hours in the period');
        ncwriteatt(f, v.name, 'units', units);
        ncwriteatt(f, v.name, 'long_name', v.desc);
        ncwriteatt(f, v.name, 'cell_methods', om_cell_method(v.map_func));
        if nc > 1; ncwriteatt(f, v.name, 'columns', strjoin(v.cols, ',')); end
        ncwriteatt(f, '/', 'title', sprintf('T&C %s, %s %s', OM.site, period, v.name));
        ncwriteatt(f, '/', 'aggregation_period', period);
    end
    if nc == 1
        ncwrite(f, v.name, permute(G, [2 1]), [1 1 k]);
    else
        ncwrite(f, v.name, permute(G, [2 1 3]), [1 1 1 k]);
    end
    ncwrite(f, 'time', days(t0 - datetime(1970,1,1)), k);
    ncwrite(f, 'n_steps', int32(n), k);
end
end

%--------------------------------------------------------------------------
function om_write_tab(OM, g)
% Writes the buffered records of a table group (points or series)
if isempty(g.buf); return; end
if strcmp(g.type, 'points'); root = OM.dir.points; else; root = OM.dir.series; end
t0 = [g.buf.t0]'; nst = [g.buf.n]';
for p = 1:numel(g.files)
    D = cell2mat(arrayfun(@(r) r.V(p,:), g.buf(:), 'UniformOutput', false)); % records x columns
    base = sprintf('%s%s/%s_%s', root, g.period, g.files{p}, g.period);
    if ismember('csv', OM.set.point_format)
        f = [base '.csv'];
        new = ~isfile(f);
        fid = fopen(f, 'a');
        if fid < 0; error('Output_Manager: cannot write %s', f); end
        if new
            if strcmp(g.period, 'hour'); fprintf(fid, 'Date,%s\n', strjoin(g.cols, ','));
            else; fprintf(fid, 'Date,n_steps,%s\n', strjoin(g.cols, ',')); end
        end
        ds = cellstr(char(datetime(t0, 'Format', 'yyyy-MM-dd HH:mm')));
        fmt = [repmat(',%.7g', 1, size(D,2)) '\n'];
        for r = 1:size(D,1)
            if strcmp(g.period, 'hour'); fprintf(fid, ['%s' fmt], ds{r}, D(r,:));
            else; fprintf(fid, ['%s,%d' fmt], ds{r}, nst(r), D(r,:)); end
        end
        fclose(fid);
    end
    if ismember('mat', OM.set.point_format)
        yrs = year(t0);
        for yr = unique(yrs)'
            ir = find(yrs == yr);
            f = sprintf('%s_%d.mat', base, yr);
            [nrec, ~] = om_rec_index(t0(ir(1)), g.period);
            if ~isfile(f)
                S.data = NaN(nrec, numel(g.cols), OM.class);
                S.time = NaN(nrec, 1);
                S.n_steps = zeros(nrec, 1);
                S.columns = g.cols;
                S.info = om_tab_info(OM, g);
                save(f, '-struct', 'S', '-v7.3');
                clear S
            end
            m = matfile(f, 'Writable', true);
            kk = zeros(numel(ir),1);
            for r = 1:numel(ir); [~, kk(r)] = om_rec_index(t0(ir(r)), g.period); end
            if all(diff(kk) == 1)
                m.data(kk(1):kk(end),:) = cast(D(ir,:), OM.class);
                m.time(kk(1):kk(end),1) = datenum(t0(ir));
                m.n_steps(kk(1):kk(end),1) = nst(ir);
            else
                for r = 1:numel(ir)
                    m.data(kk(r),:) = cast(D(ir(r),:), OM.class);
                    m.time(kk(r),1) = datenum(t0(ir(r)));
                    m.n_steps(kk(r),1) = nst(ir(r));
                end
            end
        end
    end
end
end

%--------------------------------------------------------------------------
function info = om_tab_info(OM, g)
units = {}; desc = {}; func = {};
for b = 1:numel(g.vi)
    v = OM.var(g.vi(b)); nc = numel(v.cols);
    units = [units, repmat({om_units_out(v.units, g.funcs{b}, g.period)}, 1, nc)];
    desc  = [desc,  repmat({v.desc}, 1, nc)];
    func  = [func,  repmat(g.funcs(b), 1, nc)];
end
info = struct('columns',{g.cols}, 'units',{units}, 'description',{desc}, 'func',{func}, ...
    'period',g.period, 'site',OM.site, ...
    'time_note','time = start of the period (datenum): datetime(time,''ConvertFrom'',''datenum'')');
end

%--------------------------------------------------------------------------
function om_write_info(OM)
% Outputs_info.csv: one line per exported column
fid = fopen([OM.dir.root 'Outputs_info.csv'], 'w');
fprintf(fid, 'Output,Period,Func,Variable,Column,Units,Description\n');
for k = 1:numel(OM.map)
    v = OM.var(OM.map(k).vi);
    for c = 1:numel(v.cols)
        fprintf(fid, 'Map,%s,%s,%s,%s,%s,"%s"\n', OM.map(k).period, OM.map(k).func, v.name, v.cols{c}, ...
            om_units_out(v.units, OM.map(k).func, OM.map(k).period), v.desc);
    end
end
for k = 1:numel(OM.tab)
    g = OM.tab(k);
    for b = 1:numel(g.vi)
        v = OM.var(g.vi(b));
        if strcmp(g.type, 'series'); cols = g.cols; else; cols = v.cols; end % series: one column per outlet/cell
        for c = 1:numel(cols)
            fprintf(fid, '%s,%s,%s,%s,%s,%s,"%s"\n', g.type, g.period, g.funcs{b}, v.name, cols{c}, ...
                om_units_out(v.units, g.funcs{b}, g.period), v.desc);
        end
    end
end
fclose(fid);
end

%--------------------------------------------------------------------------
function om_write_status(OM, d, t, finished)
el = toc(OM.job_clock);
nstep = t - OM.job_t_start + 1;
sp = el / max(nstep, 1);
left = (OM.N_time_step - t) * sp;
fid = fopen(OM.dir.status, 'w');
if fid < 0; return; end
fprintf(fid, 'T&C run status (%s)\n', char(datetime('now','Format','yyyy-MM-dd HH:mm:ss')));
fprintf(fid, 'Output folder      : %s\n', OM.dir.root);
fprintf(fid, 'Model date         : %s\n', char(datetime(d,'Format','yyyy-MM-dd HH:mm')));
fprintf(fid, 'Step               : %d / %d (%.1f %%)\n', t, OM.N_time_step, 100 * t / OM.N_time_step);
fprintf(fid, 'This job           : %d steps in %.2f h (%.1f s per step)\n', nstep, el/3600, sp);
if finished
    fprintf(fid, 'State              : FINISHED\n');
else
    fprintf(fid, 'Estimated time left: %.1f h\n', left/3600);
end
fprintf(fid, 'Outputs written up to the end of: %s\n', char(datetime(d,'Format','yyyy-MM-dd')));
fprintf(fid, 'Last checkpoint    : %s\n', OM.last_ckpt);
fclose(fid);
end

%--------------------------------------------------------------------------
function f = om_find_checkpoint(store, rdate)
L = dir([store 'Checkpoint_*.mat']);
if isempty(L); error('Output_Manager: Restart_ID = 1 but no checkpoint in %s', store); end
tk = NaT(numel(L),1);
for i = 1:numel(L)
    tk(i) = datetime(L(i).name(12:26), 'InputFormat', 'yyyy-MM-dd_HHmm');
end
want = NaT;
try
    if ~isempty(rdate) && ~(isnumeric(rdate) && all(isnan(rdate))); want = datetime(rdate); end
catch
    want = NaT;
end
if isnat(want)
    [~, i] = max(tk);
else
    i = find(tk == want | dateshift(tk - hours(1), 'start', 'day') == dateshift(want, 'start', 'day'), 1, 'last');
    if isempty(i)
        error('Output_Manager: no checkpoint for Restart_Date %s in %s (available: %s)', char(want), store, strjoin({L.name}, ', '));
    end
end
f = [store L(i).name];
end

%--------------------------------------------------------------------------
function om_clean_checkpoints(store, keep)
L = dir([store 'Checkpoint_*.mat']);
if numel(L) <= keep; return; end
[~, o] = sort({L.name});  % names sort in time order
L = L(o);
for i = 1:numel(L) - keep
    delete([store L(i).name]);
end
end

%--------------------------------------------------------------------------
function om_trim_outputs(OM)
% Restart: removes the records that end after the checkpoint time OM.tc
tc = OM.tc;
if ~OM.ready; return; end
% Maps
for k = 1:numel(OM.map)
    v = OM.var(OM.map(k).vi); per = OM.map(k).period;
    L = dir(sprintf('%s%s/%s/%s_%s_*.*', OM.dir.maps, per, v.name, v.name, per));
    for i = 1:numel(L)
        f = [L(i).folder '/' L(i).name];
        [~, ~, ext] = fileparts(f);
        if strcmp(ext, '.mat')
            m = matfile(f, 'Writable', true);
            tt = m.time; bad = find(~isnan(tt) & om_period_end(datetime(tt,'ConvertFrom','datenum'), per) > tc);
            if isempty(bad); continue; end
            if all(isnan(tt) | ismember((1:numel(tt))', bad)); delete(f); continue; end
            sz = size(m, v.name);
            for r = bad'
                if numel(sz) == 3; m.(v.name)(:,:,r) = NaN(sz(1), sz(2), OM.class);
                else; m.(v.name)(:,:,:,r) = NaN(sz(1), sz(2), sz(3), OM.class); end
                m.time(r,1) = NaN; m.n_steps(r,1) = 0;
            end
        elseif strcmp(ext, '.nc')
            tt = ncread(f, 'time'); tdt = datetime(1970,1,1) + days(tt);
            bad = find(~isnan(tt) & om_period_end(tdt, per) > tc);
            if isempty(bad); continue; end
            if all(isnan(tt) | ismember((1:numel(tt))', bad)); delete(f); continue; end
            info = ncinfo(f, v.name); sz = info.Size;
            for r = bad'
                if numel(sz) == 3; ncwrite(f, v.name, NaN(sz(1), sz(2), OM.class), [1 1 r]);
                else; ncwrite(f, v.name, NaN(sz(1), sz(2), sz(3), OM.class), [1 1 1 r]); end
                ncwrite(f, 'time', NaN, r); ncwrite(f, 'n_steps', int32(0), r);
            end
        end
    end
end
% Points and series
for k = 1:numel(OM.tab)
    g = OM.tab(k);
    if strcmp(g.type, 'points'); root = OM.dir.points; else; root = OM.dir.series; end
    % Records still in the write buffer of the checkpoint are written again
    % after the restart: remove them from the csv files (appended files)
    tbuf = NaT;
    if ~isempty(g.buf); tbuf = min([g.buf.t0]); end
    for p = 1:numel(g.files)
        base = sprintf('%s%s/%s_%s', root, g.period, g.files{p}, g.period);
        f = [base '.csv'];
        if isfile(f)
            txt = splitlines(fileread(f));
            txt = txt(~cellfun(@isempty, txt));
            keep = true(numel(txt), 1);
            for r = 2:numel(txt)
                dt = datetime(txt{r}(1:16), 'InputFormat', 'yyyy-MM-dd HH:mm');
                keep(r) = om_period_end(dt, g.period) <= tc && ~(dt >= tbuf);
            end
            if any(~keep)
                fid = fopen(f, 'w'); fprintf(fid, '%s\n', txt{keep}); fclose(fid);
            end
        end
        L = dir([base '_*.mat']);
        for i = 1:numel(L)
            f = [L(i).folder '/' L(i).name];
            m = matfile(f, 'Writable', true);
            tt = m.time; bad = find(~isnan(tt) & om_period_end(datetime(tt,'ConvertFrom','datenum'), g.period) > tc);
            if isempty(bad); continue; end
            if all(isnan(tt) | ismember((1:numel(tt))', bad)); delete(f); continue; end
            nc = size(m, 'data', 2);
            m.data(bad, :) = NaN(numel(bad), nc, OM.class);
            m.time(bad, 1) = NaN; m.n_steps(bad, 1) = 0;
        end
    end
end
end

%--------------------------------------------------------------------------
function cat = om_catalogue()
% Units and descriptions of the usual output variables
% (variables not listed here can still be exported; their units show as '?')
cat = {
 'ANPP_H','gC/m2/day','Above-ground net primary production, high vegetation'
 'ANPP_L','gC/m2/day','Above-ground net primary production, low vegetation'
 'LAI_H','m2/m2','Leaf area index, high vegetation'
 'LAI_L','m2/m2','Leaf area index, low vegetation'
 'NDVI','-','Normalised difference vegetation index'
 'NPP_H','gC/m2/day','Net primary production, high vegetation'
 'NPP_L','gC/m2/day','Net primary production, low vegetation'
 'hc_H','m','Vegetation height, high vegetation'
 'hc_L','m','Vegetation height, low vegetation'
 'An_H','umolCO2/m2/s','Net assimilation, high vegetation'
 'An_L','umolCO2/m2/s','Net assimilation, low vegetation'
 'RA_H','gC/m2/day','Autotrophic respiration, high vegetation'
 'RA_L','gC/m2/day','Autotrophic respiration, low vegetation'
 'Rg_H','gC/m2/day','Growth respiration, high vegetation'
 'Rg_L','gC/m2/day','Growth respiration, low vegetation'
 'SAI_H','m2/m2','Stem area index, high vegetation'
 'SAI_L','m2/m2','Stem area index, low vegetation'
 'PHE_S_H','-','Phenology state, high vegetation'
 'PHE_S_L','-','Phenology state, low vegetation'
 'EG','mm/h','Evaporation from bare soil'
 'EICE','mm/h','Evaporation/sublimation from ice'
 'EIn_H','mm/h','Evaporation of intercepted water, high vegetation'
 'EIn_L','mm/h','Evaporation of intercepted water, low vegetation'
 'EIn_rock','mm/h','Evaporation of water intercepted by rock'
 'EIn_urb','mm/h','Evaporation of water intercepted by urban surfaces'
 'EWAT','mm/h','Evaporation from water and ponds'
 'T_H','mm/h','Transpiration, high vegetation'
 'T_L','mm/h','Transpiration, low vegetation'
 'ESN','mm/h','Evaporation/sublimation from the snowpack'
 'ELitter','mm/h','Evaporation from litter'
 'ESN_In','mm/h','Evaporation/sublimation of intercepted snow'
 'SND','m','Snow depth'
 'SWE','mm','Snow water equivalent'
 'ros','kg/m3','Snow density'
 'snow_albedo','-','Snow albedo'
 'Smelt','mm/h','Snow melt'
 'Imelt','mm/h','Ice melt'
 'SP_wc','mm','Snowpack liquid water content'
 'SWE_avalanched','mm','SWE change due to avalanches'
 'Cicew','-','Presence of frozen water (0/1)'
 'Cice','-','Presence of ice (0/1)'
 'IP_wc','mm','Ice pack liquid water content'
 'ICE_D','m','Ice thickness'
 'ICE','mm','Ice water equivalent'
 'U_SWE','mm','Unloaded snow from intercepted snow'
 'FROCK','mm','Water storage in fractured rock'
 'Gfin','W/m2','Ground heat flux (heat diffusion)'
 'Lk_rock','mm/h','Leakage from rock surface to bedrock'
 'Lk','mm/h','Leakage at the bottom of the soil (recharge)'
 'Lk_wat','mm/h','Leakage from water/ponds to bedrock'
 'ZWT','mm','Water table depth'
 'O','m3/m3','Soil water content (per layer)'
 'OF','m3/m3','Soil water content, first layer'
 'OH','m3/m3','Soil water content available to roots, high vegetation'
 'OL','m3/m3','Soil water content available to roots, low vegetation'
 'OS','m3/m3','Soil water content of the bare-soil evaporation layers'
 'Tdp','degC','Soil temperature (per layer)'
 'Tdp_L','degC','Soil temperature of the root zone, low vegetation'
 'Tdp_H','degC','Soil temperature of the root zone, high vegetation'
 'Q_channel','mm','Water in the channel, as a depth over the cell (x cellsize/channel width for the water depth in the channel)'
 'QpointC','mm/h','Channel outflow at the outlet and the points, as a depth over the cell (x cellsize^2/3.6e6 for m3/s)'
 'QpointH','mm/h','Hillslope surface outflow at the outlet and the points, as a depth over the cell (x cellsize^2/3.6e6 for m3/s)'
 'UpointC','m/s','Channel flow velocity at the outlet and the points'
 'UpointH','m/s','Hillslope surface flow velocity at the outlet and the points'
 'Q_exit','mm/h','Surface water leaving the domain, as a depth over the catchment'
 'Qsub_exit','mm/h','Subsurface water leaving the domain, as a depth over the catchment'
 'q_runon','mm/h','Run-on'
 'Rd','mm/h','Saturation excess runoff'
 'Rh','mm/h','Infiltration excess runoff'
 'f','mm/h','Infiltration'
 'In_rock','mm','Water intercepted by rock'
 'In_SWE','mm','Intercepted snow water equivalent'
 'In_L','mm','Intercepted water, low vegetation'
 'In_H','mm','Intercepted water, high vegetation'
 'Pr_liq','mm/h','Liquid precipitation'
 'Pr_sno','mm/h','Solid precipitation'
 'Pr_St','mm/h','Precipitation (forcing)'
 'Ta_St','degC','Air temperature (forcing)'
 'Tdew_St','degC','Dew point temperature (forcing)'
 'Ws_St','m/s','Wind speed (forcing)'
 'U_St','-','Relative humidity, fraction (forcing)'
 'N_St','W/m2','Incoming longwave radiation (forcing)'
 'ea_St','Pa','Vapour pressure (forcing)'
 'Pre_St','hPa','Atmospheric pressure (forcing)'
 'SAD1_St','W/m2','Diffuse shortwave radiation, band 1 (forcing)'
 'SAD2_St','W/m2','Diffuse shortwave radiation, band 2 (forcing)'
 'SAB1_St','W/m2','Direct shortwave radiation, band 1 (forcing)'
 'SAB2_St','W/m2','Direct shortwave radiation, band 2 (forcing)'
 'Ts','degC','Surface temperature (energy balance)'
 'Rn','W/m2','Net radiation'
 'QE','W/m2','Latent heat flux'
 'H','W/m2','Sensible heat flux'
 'HV','W/m2','Sensible heat flux from vegetation above snow'
 'G','W/m2','Ground heat flux'
 'Qfm','W/m2','Heat for freezing/melting'
 'Qv','W/m2','Heat advected by precipitation'
 'dQ_S','W/m2','Residual of the surface energy balance'
 'DQ_S','W/m2','Residual of the surface energy balance'
 'dQVEG','W/m2','Residual of the energy balance of snow-free vegetation'
 'Csno','-','Presence of snow (0/1)'
};
end
