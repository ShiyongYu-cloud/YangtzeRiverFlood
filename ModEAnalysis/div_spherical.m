function [LON,LAT,DIV] = div_spherical(lon,lat,U,V)
%function for calculating divergence of a vector field on a spherical surface 
%note that grids can be unevenly spaced
%INPUT
%lon: scalar containing longitudinal coordinates
%lat: scalar containing latitudinal coordinates
%U: matrix containing the zonal component of the vector field
%V: matrix containing the meridional component of the vector field
%OUTPUT
%LON: matrix containing longitudinal grid points
%LAT: matrix containing latitudinal grids points
%DIV: matrix containing divergence
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
if size(U,1) ~= M || size(U,2) ~= N
    U = U';
end
if size(V,1) ~= M || size(V,2) ~= N
    V = V';
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
%% calculate divergence
DIV = zeros(M,N);
for i = 1:M
    for j = 1:N
        if (i == 1) && (j == 1) %upperleft corner 
            DUDX = (U(i,j+1)-U(i,j))/DX(i,j);
            DVDY = (V(i+1,j)-V(i,j))/DY(i,j);
            DIV(i,j) = DUDX + DVDY;
        elseif (i == 1) && (j == N) %upperright corner
            DUDX = (U(i,j)-U(i,j-1))/DX(i,j-1);
            DVDY = (V(i+1,j)-V(i,j))/DY(i,j);
            DIV(i,j) = DUDX + DVDY;
        elseif (i == M) && (j == 1) %lowerleft corner
            DUDX = (U(i,j+1)-U(i,j))/DX(i,j);
            DVDY = (V(i,j)-V(i-1,j))/DY(i-1,j);
            DIV(i,j) = DUDX + DVDY;
        elseif (i == M) && (j == N) %lowerright corner
            DUDX = (U(i,j)-U(i,j-1))/DX(i,j-1);
            DVDY = (V(i,j)-V(i-1,j))/DY(i-1,j);
            DIV(i,j) = DUDX + DVDY;
        elseif (i == 1) && (j > 1 && j < N) %upper edge 
            DUDX = 0.5*(U(i,j+1)-U(i,j))/DX(i,j)+0.5*(U(i,j)-U(i,j-1))/DX(i,j-1);
            DVDY = (V(i+1,j)-V(i,j))/DY(i,j);
            DIV(i,j) = DUDX + DVDY;
        elseif (i == M) && (j > 1 && j < N) %lower edge
            DUDX = 0.5*(U(i,j+1)-U(i,j))/DX(i,j)+0.5*(U(i,j)-U(i,j-1))/DX(i,j-1);
            DVDY = (V(i,j)-V(i-1,j))/DY(i-1,j);
            DIV(i,j) = DUDX + DVDY;
        elseif (j == 1) && (i > 1 && i < M) %left edge
            DUDX = (U(i,j+1)-U(i,j))/DX(i,j);
            DVDY = 0.5*(V(i+1,j)-V(i,j))/DY(i,j)+0.5*(V(i,j)-V(i-1,j))/DY(i-1,j);
            DIV(i,j) = DUDX + DVDY;
        elseif (j == N) && (i > 1 && i < M) %right edge
            DUDX = (U(i,j)-U(i,j-1))/DX(i,j-1);
            DVDY = 0.5*(V(i+1,j)-V(i,j))/DY(i,j)+0.5*(V(i,j)-V(i-1,j))/DY(i-1,j);
            DIV(i,j) = DUDX + DVDY;            
        else  %inner grids            
            DUDX = 0.5*(U(i,j+1)-U(i,j))/DX(i,j)+0.5*(U(i,j)-U(i,j-1))/DX(i,j-1);
            DVDY = 0.5*(V(i+1,j)-V(i,j))/DY(i,j)+0.5*(V(i,j)-V(i-1,j))/DY(i-1,j);
            DIV(i,j) = DUDX + DVDY;            
        end 
    end
end
end