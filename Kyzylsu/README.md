# Kyzylsu – set-up for `Launcher_Dev_Distributed.m`

Conversion of the former set-up `Case_study/Kyzylsu_distributed` (`Launcher_Kyzylsu_Distributed.m`,
`Inputs/PARAMETERS_SOIL.m`, `Inputs/INIT_COND_v2.m`) to the structure used for Shallap.

```
Kyzylsu/
  Parameters/      Site_Config_TC.xlsx, Options_and_NoVegParams.xlsx, Parameters_TC.xlsx, VPAR_T.xlsx
  Forcing/         Pluvio.mat, AWS_3900m.mat, AWS_OnGlacier.mat (original station data), Ca_Data.mat,
                   Kyzylsu_Pluvio_meteo.mat, Kyzylsu_hm_Ta_lapse.mat, Kyzylsu_Pr_ratio_v_el_const.mat
  Preprocessing/   Prepare_Kyzylsu_inputs.m, INPUTS/ (former dtm, albedo profile, initial snow maps),
                   OUTPUTS/Kyzylsu_100m_from_Kyzylsu_distributed/dtm_Kyzylsu_100m.mat
  Multipoint/      OUTPUTS/Kyzylsu_MultiPoints.txt
```

`Preprocessing/Prepare_Kyzylsu_inputs.m` regenerates everything in `Forcing/` (except the original
station files), `Preprocessing/OUTPUTS/` and `Multipoint/OUTPUTS/` (run it from `Preprocessing/`,
with the repository `Functions` folder on the path for `utm2deg`).

## Land-cover classes (VEG_CODE in the new dtm file)

The classes are the same as in the former `PARAMETERS_SOIL.m` (`ANSWER` 1–7). Ice cells (former code 8) are set to
7 (Rock), as the former launcher did (`ksv(ksv==8)=7`).

| Class | Veg_type | Ccrown | Crock | Cbare | Krock | Parameters |
|---|---|---|---|---|---|---|
| 1 | Fir | 1 | 0 | 0 | NaN | Former column 1 |
| 2 | Larch | 1 | 0 | 0 | NaN | Former column 2 |
| 3 | Grass_C3 | 1 | 0 | 0 | NaN | Former column 3 |
| 4 | Shrub | 0.9 | 0 | 0.1 | NaN | Former column 4 |
| 5 | BL_ever | 1 | 0 | 0 | NaN | Former column 5 |
| 6 | BL_dec | 1 | 0 | 0 | NaN | Former column 6 |
| 7 | Rock (and Ice) | 0 | 1 | 0 | 0.15 | Former column 4 (Shrub) without roots |

Cwat and Curb are 0 and Kbot is 0 for all classes. Classes present in the catchment: 2 (4 cells), 3 (8178),
6 (13) and 7 (8033, including 3372 ice cells).

`INIT_COND_v6` sets `WATtm1 = 50000` mm for class 3 (lake water, as for Shallap), which here is Grass. With
Cwat = 0, `HYDROLOGIC_UNIT` sets WAT to 0 at the first step and no lake evaporation or freezing is computed,
so the grass cells are not affected. The only effect is a jump in the catchment water storage (and its balance check)
at the first time step.

`OPT_Veg_Param = 2`: the Ccrown/Cwat/Curb/Crock/Cbare fractions are written as grids (CCROWN, CWATER, CURB,
CROCK, CBARE) in the dtm file with the same values as `VPAR_T.xlsx`. With `OPT_Veg_Param = 1` the launcher fails
in the output stage, because `VPAR.Ccrownr` is a cell array there.

## Differences from the former set-up (from the new launcher, which was not modified)

- **Air temperature**: the former set-up used an hourly lapse rate from Pluvio (3371 m) and AWS_3900m. Now the launcher uses a
  climatological lapse rate per hour and month (mean of those hourly rates, `Kyzylsu_hm_Ta_lapse.mat`).
  Dew point is lapsed with the same rate, so RH now varies with elevation (it was uniform before).
- **Precipitation**: the ratio is uniform (p1 = 0, p2 = 1), the same as the former Pmod = 0.
- **Shortwave**: the station files hold already-partitioned SW. `SWin` = SAB1+SAB2+SAD1+SAD2 is
  re-partitioned by `Automatic_Radiation_Partition_I`, and the launcher applies its fixed SAB/PARB elevation gradients.
- **Longwave and pressure**: the launcher applies an LWin elevation gradient (–0.031 W m-2 m-1) and computes pressure from
  elevation. The former set-up used the Pluvio values everywhere.
- **Bare-ice albedo**: the former piecewise profile (`Kyzylsu_Albedo_vs_elev.mat`, capped at 0.4) is replaced by a linear fit
  (Alb_v_el_MM/CC), bounded by Aice = 0.09 and Amax = 0.4.
- **Soil layers**: the former set-up used 10 layers, Zs = [0 10 20 50 100 150 200 300 400 700 1000] mm, with the deepest layer
  impermeable (SOIL_TH ≥ 100 % everywhere). The launcher builds Zs from `Zs_all` up to max(SOIL_TH), and its point outputs
  (`BC_O`) assume exactly 10 layers. A 100 cm soil (11 layers) stops the run, so SOIL_TH is capped at 99 cm
  (`soil_th_max_cm`). This gives Zs = [0 10 20 50 100 150 200 250 300 350 500] mm, with the 350–500 mm layer impermeable.
  The soil column is shallower than before. Use 100 cm once the launcher's point outputs handle `ms_max` ≠ 10.
  The initial soil moisture (`vi`) is the launcher default.
- **Vegetation initial conditions**: the per-class values from `INIT_COND_v2` are in `Parameters_TC.xlsx` (Ci = Ca on
  01-Oct-2021 = 413.19 ppm). The carbon pools (B_H, B_L) are set globally inside `INIT_COND_v6`, so they cannot be set per class.
- **Grass cutting**: the former `Mpar_L.jDay_cut = 180:243; LAI_cut = -0.1` for grass is not available in `PARAMETERS_SOIL_DEV`.
- **Initial snow**: SNOWD/SNOWALB come from `Kyzylsu_Init_*_LSS2_1999-10-16.mat` (former `diff_IniSND = 1`). Set
  `use_ini_snow` or the file names in the preprocessing script to change this.
