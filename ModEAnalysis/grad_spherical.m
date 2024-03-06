function [LON,LAT,GFX,GFY] = grad_spherical(lon,lat,F)
%function for calculating gradient of a scalar field on a spherical surface 
%note that grids can be unevenly spaced
%INPUT
%lon: scalar containing longitudinal coordinates
%lat: scalar containing latitudinal coordinates
%F:   matrix containing gridded scalar field
%OUTPUT
%LON: matrix containing longitudinal grid points
%LAT: matrix containing latitudinal grid points
%GFX: matrix containing the zonal component of gradient of F 
%GFY: matrix containing the meridional component of gradient of F
%% make data grids
[LON,LAT] = meshgrid(lon,lat);
%% check dimension
M = numel(lat);
N = numel(lon);
if size(lat,1) ~= M
    lat = lat';
end
if size(lon,2) ~= N
    lon = lon';
end
if size(F,1) ~= M || size(F,2) ~= N
    F = F';
end
%% make computing grids
ra = earth_radius('m');% earth radus in meters
d_lon = abs(diff(lon));
d_lat = abs(diff(lat));
dx = ra*deg2rad(d_lon);
dy = ra*deg2rad(d_lat); 
WT = repmat(cosd(lat),1,N-1);% latitude weight
DX = repmat(dx,M,1);
DX = DX.*WT; %weight longitudinal grids by latitude
DY = repmat(dy,1,N);
%% calculate the zonal and meridional components of gradient
GFX = zeros(M,N);
GFY = zeros(M,N);
for i = 1:M
    for j = 1:N
        if (i == 1) && (j == 1) %upperleft corner 
            DFDX = (F(i,j+1)-F(i,j))/DX(i,j);
            DFDY = (F(i+1,j)-F(i,j))/DY(i,j);
            GFX(i,j) = DFDX;
            GFY(i,j) = DFDY;
        elseif (i == 1) && (j == N) %upperright corner
            DFDX = (F(i,j)-F(i,j-1))/DX(i,j-1);
            DFDY = (F(i+1,j)-F(i,j))/DY(i,j);
            GFX(i,j) = DFDX;
            GFY(i,j) = DFDY;
        elseif (i == M) && (j == 1) %lowerleft corner
            DFDX = (F(i,j+1)-F(i,j))/DX(i,j);
            DFDY = (F(i,j)-F(i-1,j))/DY(i-1,j);
            GFX(i,j) = DFDX;
            GFY(i,j) = DFDY;
        elseif (i == M) && (j == N) %lowerright corner
            DFDX = (F(i,j)-F(i,j-1))/DX(i,j-1);
            DFDY = (F(i,j)-F(i-1,j))/DY(i-1,j);
            GFX(i,j) = DFDX;
            GFY(i,j) = DFDY;
        elseif (i == 1) && (j > 1 && j < N) %upper edge 
            DFDX = 0.5*(F(i,j+1)-F(i,j))/DX(i,j)+0.5*(F(i,j)-F(i,j-1))/DX(i,j-1);
            DFDY = (F(i+1,j)-F(i,j))/DY(i,j);
            GFX(i,j) = DFDX;
            GFY(i,j) = DFDY;
        elseif (i == M) && (j > 1 && j < N) %lower edge
            DFDX = 0.5*(F(i,j+1)-F(i,j))/DX(i,j)+0.5*(F(i,j)-F(i,j-1))/DX(i,j-1);
            DFDY = (F(i,j)-F(i-1,j))/DY(i-1,j);
            GFX(i,j) = DFDX;
            GFY(i,j) = DFDY;
        elseif (j == 1) && (i > 1 && i < M) %left edge
            DFDX = (F(i,j+1)-F(i,j))/DX(i,j);
            DFDY = 0.5*(F(i+1,j)-F(i,j))/DY(i,j)+0.5*(F(i,j)-F(i-1,j))/DY(i-1,j);
            GFX(i,j) = DFDX;
            GFY(i,j) = DFDY;
        elseif (j == N) && (i > 1 && i < M) %right edge
            DFDX = (F(i,j)-F(i,j-1))/DX(i,j-1);
            DFDY = 0.5*(F(i+1,j)-F(i,j))/DY(i,j)+0.5*(F(i,j)-F(i-1,j))/DY(i-1,j);
            GFX(i,j) = DFDX;
            GFY(i,j) = DFDY;
        else  %inner grids            
            DFDX = 0.5*(F(i,j+1)-F(i,j))/DX(i,j)+0.5*(F(i,j)-F(i,j-1))/DX(i,j-1);
            DFDY = 0.5*(F(i+1,j)-F(i,j))/DY(i,j)+0.5*(F(i,j)-F(i-1,j))/DY(i-1,j);
            GFX(i,j) = DFDX;
            GFY(i,j) = DFDY;
        end 
    end
end
end