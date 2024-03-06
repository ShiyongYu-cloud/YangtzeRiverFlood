function [LON,LAT,FX,FY,CF] = moisture_transport(lon,lat,lev,U,V,Q)
%function for calculating vertically integrated moisture flux and
%convergence
%INPUT
%lon: vector containing longitudinal coordinates
%lat: vector containing latitudinal coordinates
%lev: vector containing pressure levels from bottom to top (e.g., [1000,...,0]) in Pa 
%U:   matix containing zonal winds 
%V:   matix containing meridional winds 
%Q:   matix containing specific humidity
%OUTPUT
%LON: matrix containing longitudinal grids 
%LAT: matrix containing latitudinal grids 
%FX:  matrix containing longitudinal component of vertically integrated moisture flux 
%FY:  matrix containing latitudinal component of vertically integrated moisture flux 
%CF:  matrix containing vertically integrated moisture flux convergence 
%%
%% make grids
[LON,LAT] = meshgrid(lon,lat);
%% check dimension
M = numel(lat);
N = numel(lon);
L = numel(lev);
if size(U,1) ~= M || size(U,2) ~= N
    U = permute(U,[2 1 3]);        %transpose to lat-lon-lev
end
if size(V,1) ~= M || size(V,2) ~= N
    V = permute(V,[2 1 3]);        %transpose to lat-lon-lev
end
if size(Q,1) ~= M || size(Q,2) ~= N
    Q = permute(Q,[2 1 3]);        %transpose to lat-lon-lev
end
%% standard acceleration due to gravity 
g0 = 9.80665; % m/s2
%% convert millibar (hPa) to Pa (kg/m2.s) for pressure
lev = 100*lev; 
%% calculate moisture flux for each pressure level (kg.m/kg.s)
Q_U = Q.*U;
Q_V = Q.*V;
%% integrate moisture flux vertically using Simpson's rule(kg/m.s)
Q_U = reshape(Q_U,M*N,L); %convert to lat*lon,lev
Q_V = reshape(Q_V,M*N,L); %convert to lat*lon,lev
if any(diff(lev)<0)
    FX = -simps(lev,Q_U,2)/g0; 
    FY = -simps(lev,Q_V,2)/g0;
else
    FX = simps(lev,Q_U,2)/g0; 
    FY = simps(lev,Q_V,2)/g0;
end
FX = reshape(FX,M,N); %reshape back to lat-lon
FY = reshape(FY,M,N); %reshape back to lat-lon
%% calculate moisture flux divergence for each pressure level (kg/kg.s)
q_divv = zeros(M,N,L);
v_graq = zeros(M,N,L);
for i = 1:L
    [~,~,divv] = div_spherical(lon,lat,U(:,:,i),V(:,:,i));  %divergence of V
    [~,~,gqx,gqy] = grad_spherical(lon,lat,Q(:,:,i));       %gradient of q 
    q_divv(:,:,i) = Q(:,:,i).*divv;
    v_graq(:,:,i) = U(:,:,i).*gqx + V(:,:,i).*gqy;
end
%% integrate moisture flux divergence vertically using Simpson's rule (kg/m2.s)
q_divv = reshape(q_divv,M*N,L); %convert to lat*lon,lev
v_graq = reshape(v_graq,M*N,L); %convert to lat*lon,lev
if any(diff(lev)<0)
    DF1 = -simps(lev,q_divv,2)/g0; 
    DF2 = -simps(lev,v_graq,2)/g0;
else
    DF1 = simps(lev,q_divv,2)/g0; 
    DF2 = simps(lev,v_graq,2)/g0;
end
DF = DF1 + DF2;
DF = reshape(DF,M,N); %reshape back to lat-lon
%% convert divergence to convergence
CF = -DF; 
end