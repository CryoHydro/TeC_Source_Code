%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute Soil Heat Flux     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[G0,T,Gn]=Soil_Heat_Profile_Normal(Tup,dt,Ttm1,ms,Zs,lanS,cv,Tdown,G0,Gn,OPZ)
%%%INPUTS
%%%%%% Boundary Mixed Condition
%Tdown = NaN; %% [K]
%G0 = NaN; %% [W/m^2]
%%%%
% OPZ case 1
% %%%% G0 Gn
% %%%%%%%%%%%%%
% OPZ case 2
% %%% Tup Tdown
% %%%%%%%%%%%%%%
% OPZ  case 3
% %%% Tup Gn
% OPZ case 4 
%%%% G0 Tdown 
% %%%%%%%%%%%%%%
dz= diff(Zs); %%%% [mm]  Thickness of the Layers
%%%%
if (OPZ == 1) || (OPZ == 3)
    Tdown= Ttm1(ms) - Gn*(0.001*dz(ms)*0.5)/lanS(ms); %%[K]
end
if (OPZ == 1) || (OPZ == 4)
    Tup = Ttm1(1) + G0*(0.001*dz(1)*0.5)/lanS(1);%% [K]
end
%%%%%%%%%%
Ttm1= [ Tup Ttm1 Tdown];
%%%%%%%
nit=10; %%% Number of Internal time step
RF=0;
%G0p = zeros(1,nit); Gnp = zeros(1,nit);
%%%


%%%% this is moved out of the Heat_CN_normal function - constants computed outside of iteration now
Ttm1=flip(Ttm1); lanS=flip(lanS);  cv=flip(cv);
lanS=reshape(lanS,length(lanS),1);
cv=reshape(cv,length(cv),1);
Ttm1=reshape(Ttm1,length(Ttm1),1);
Zs=sort(-Zs');
nz=length(Zs)-1; %% number of layers
nn = nz + 2; %%% number of nodes for temperature
dz_m= diff(Zs)*0.001;  %% [m]  Thickness of the Layers

DzC = [dz_m(1)*0.5 ; 0.5*(dz_m(1:nz-1)+dz_m(2:nz)) ; dz_m(end)*0.5];

%%% Harmonic Mean of lanS
%lanS_half(2:nz)= (dz(1:nz-1) + dz(2:nz))./(dz(1:nz-1)./(lanS(1:nz-1)) + dz(2:nz)./lanS(2:nz));
lanS_half = zeros(nz+1,1);
lanS_half(2:nz)= 0.5*(lanS(1:nz-1) + lanS(2:nz));
lanS_half(1)=lanS(1); 
lanS_half(nz+1)=lanS(nz);

%%% i ...  2:nn-1
%a = % subdiagonal a: coefficients of T(i-1)  %[1/s]
%c = a; % superdiagonal c: coefficients of T(i+1) %[1/s]

alp=0.5;   %% Tridiagonal system of Equation
dts=dt/nit;  % timestep

a = -(1-alp)*lanS_half(1:nz)./(cv(1:nz).*dz_m(1:nz).*DzC(1:nz));
c= -(1-alp)*lanS_half(2:nz+1)./(cv(1:nz).*dz_m(1:nz).*DzC(2:nz+1));

a_glob=[ NaN ; a ; NaN] ; 
c_glob=[ NaN ; c ; NaN] ; 
b_glob=[ NaN ; (1/dts)*ones(nz,1) - (a+c) ; NaN] ;

%preallocate d
d=zeros(nn,1);
%pre-compute slices
interior = 2:nn-1;
left = 1:nn-2;
right = 3:nn;


% Apply boundary conditions to coefficient matrices (before inner loop)
switch OPZ
    case 1
        %%%%%%%% Boundary Condition - Neumann 
        b_glob(nn)= 1;  a_glob(nn)= -1;
        b_glob(1)=-1;   c_glob(1)=1;

        % d boundary values use these constants, but updated inside loop
        d_nn_const = G0*DzC(nz+1)/lanS_half(nz+1);
        d_1_const = Gn*DzC(1)/lanS_half(1);
        
    case 2
        %%%%%%%% Boundary Condition - Dirichlet 
        b_glob(nn)= 1;  a_glob(nn)= 0;
        b_glob(1)=1;    c_glob(1)=0;
        % d boundary values use these constants
        d_nn_const = Tup;
        d_1_const = Tdown;
        
    case 3
        %%%%%%%%% Mixed Condition 
        b_glob(nn)= 1;  a_glob(nn)= 0;
        b_glob(1)=-1;   c_glob(1)=1;

        % d boundary values use these constants
        d_nn_const = Tup;
        d_1_const = Gn*DzC(1)/lanS_half(1);
        
    case 4
        %%%%%%%%% Mixed Condition 
        b_glob(nn)= 1;  a_glob(nn)= -1;
        b_glob(1)=1;    c_glob(1)=0;

        % d boundary values use these constants
        d_nn_const = G0*DzC(nz+1)/lanS_half(nz+1);
        d_1_const = Tdown;
end


a_int = a_glob(interior);  
c_int = c_glob(interior);

% compute coefficent - used in loop
ac_int=a_int + c_int;



% Build sparse matrix (for tridiagonal system)
A= sparse(1:nn,1:nn,b_glob,nn,nn)...
  +sparse(1:(nn-1),2:nn,c_glob(1:nn-1),nn,nn)...
  +sparse(2:nn,1:(nn-1),a_glob(2:nn),nn,nn);



%%% INNER LOOP- reduce unnecessary computation inside 
for kk=1:nit 
    % Compute full RHS
    % DON'T FACTOR OUT Ttm1(interior).*(ac_int + 1/dts) - this causes numerical differences
    d = Ttm1/dts - [NaN; a_int.*Ttm1(left); NaN] ...
        + [NaN; (ac_int).*Ttm1(interior); NaN] ...
        - [NaN; c_int.*Ttm1(right); NaN]+ [NaN ; RF./(cv.*dz_m) ; NaN]; 


    

    % Then OVERRIDE with boundary conditions (matching old code)
    d(1)=d_1_const;
    d(nn)=d_nn_const;
    
    % Solve
    Ttm1 = A\d;         
end

Tout = flip(Ttm1);   % flip back to match shapes in the previous version of this function
lanS = flip(lanS);
T=Tout(interior);% Temperature Layer %% [�C ]
T=T';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
if isnan(sum(T))
    %%%%%%%%%%%%%%%
    disp('NaN values in Soil Temperature')
    return
end
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%
if  (OPZ == 2) || (OPZ == 3)
    G0=-lanS(1).*((T(1)-Tup)./(0.001*dz(1)*0.5));  % [W/m^2]
    %G0=mean(G0p);
end
if  (OPZ == 2) || (OPZ == 4)
    Gn=-lanS(ms)*((Tdown-T(ms))./(0.001*dz(end)*0.5));
    %Gn=mean(Gnp);
end
%Tdown= T(nz) - Gn*(0.001*dz(nz)*0.5)/lanS(nz); %%[�C][K]
return
end

