function R = forcing_terrain_radiation(R, h_S, zeta_S, DTM, Slo_top, Aspect, SvFn, Ct, HZ, Zasp, rho_g)
%FORCING_TERRAIN_RADIATION Shortwave forcing of one time step corrected for
% the terrain: slope and aspect, cast shadows, sky view and reflection from
% the surrounding terrain.
%
% R: struct with the shortwave components of one time step on a flat surface,
%    each a row vector with one value per cell (1 x num_cell) [W m-2]:
%    SAB1, SAB2, PARB (direct) and SAD1, SAD2, PARD (diffuse)
% h_S, zeta_S: solar height and azimuth [rad] (SetSunVariables)
% Grids (m_cell x n_cell): DTM [m], Slo_top [fraction], Aspect [deg from N],
%    Ct (terrain configuration factor), HZ, Zasp (horizon angles and
%    their azimuths, for Shadow_Effect)
% SvFn: sky view factor, row vector (1 x num_cell)
% rho_g: albedo of the surrounding terrain
%
% Returns R with the same fields, corrected (negative and NaN values set to
% 0; all components 0 when the sun is below about 5.7 deg).

num_cell = numel(DTM);
ShF = Shadow_Effect(DTM, h_S, zeta_S, HZ, Zasp); % cast shadow factor (grid)

% Cosine of the angle between the sun and the normal of the slope
cos_fst = cos(atan(Slo_top))*sin(h_S) + sin(atan(Slo_top)).*cos(h_S).*cos(zeta_S-Aspect*pi/180);
cos_fst(cos_fst<0) = 0;

Ctn = reshape(Ct, num_cell, 1)';
ShFn = reshape(ShF, num_cell, 1)';
cos_fstn = reshape(cos_fst, num_cell, 1)';

% Diffuse: sky view + terrain reflection
R.SAD1 = R.SAD1.*SvFn + Ctn.*rho_g.*(R.SAD1./sin(h_S).*cos_fstn + (1-SvFn).*R.SAD1);
R.SAD2 = R.SAD2.*SvFn + Ctn.*rho_g.*(R.SAD2./sin(h_S).*cos_fstn + (1-SvFn).*R.SAD2);
R.PARD = R.PARD.*SvFn + Ctn.*rho_g.*(R.PARD./sin(h_S).*cos_fstn + (1-SvFn).*R.PARD);
% Direct: slope and aspect, cast shadows
R.SAB1 = (R.SAB1./sin(h_S)).*cos_fstn.*ShFn;
R.SAB2 = (R.SAB2./sin(h_S)).*cos_fstn.*ShFn;
R.PARB = (R.PARB./sin(h_S)).*cos_fstn.*ShFn;

R.SAB1(R.SAB1<0) = 0;
R.SAB2(R.SAB2<0) = 0;
R.PARB(R.PARB<0) = 0;
R.PARD(R.PARD<0) = 0;
R.SAD1(R.SAD1<0) = 0;
R.SAD2(R.SAD2<0) = 0;
R.SAB1(isnan(R.SAB1)) = 0;
R.SAB2(isnan(R.SAB2)) = 0;
R.SAD1(isnan(R.SAD1)) = 0;
R.SAD2(isnan(R.SAD2)) = 0;
R.PARB(isnan(R.PARB)) = 0;
R.PARD(isnan(R.PARD)) = 0;

% Sun below ~5.7 deg: the slope correction (division by sin(h_S)) is unstable;
% set all shortwave to 0, as in the former launcher
if sin(h_S) <= 0.10
    R.SAB1(:) = 0; R.SAB2(:) = 0; R.SAD1(:) = 0; R.SAD2(:) = 0; R.PARB(:) = 0; R.PARD(:) = 0;
end
end
