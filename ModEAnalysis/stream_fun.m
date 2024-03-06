function [LON,LAT,Psi] = stream_fun(lon,lat,U,V)
%% function for calculating stream function from a velocity field
%% on the spherical surface
%% INPUT
%lon: array containing evenly or unevently spaced longitudinal grid points
%lat: array containing evenly or unevently spaced latitudinal grid points
%U: matrix containing the zonal component of the velocity field
%V: matrix containing the meridional component of the velocity field
%% OUTPUT
%LON: matrix containing longitudinal grid points
%LAT: matrix containing latitudinal grids points
%Psi: matrix containing the stream function 
%%
%% make data grids
[LON,LAT] = meshgrid(lon,lat);
M = numel(lat);
%% integrating U over latitudinal grid points
Psi_u = cumsimps_lat(lat,U);
%% integrating V over longitudinal grid points
Psi_v = cumsimps_lon(lon,lat(1),V(1,:));
Psi_v = repmat(Psi_v,M,1);
%% combine U and V 
Psi = Psi_u - Psi_v;
return;
%%
function F = cumsimps_lat(lat,Z)
%% function for numerical approximation of a continuous function f(y) 
%% over latitudes on a spherical surface such that z(y) = df/dy using 
%% Simpson-rule column-wise cumulative summation 
%% INPUT
%lat: array containing evenly or unevently spaced latitudinal grid points
%Z:   array or matrix each COLUMN of which represents the value of the 
%     integrand at the latitudinal grid points y 
%% OUTPUT
%F: array or matrix containing the approximation of the integral F(y) 
%% determine the size of y and check the number of rows for Z 
M = numel(lat);
if size(lat,1) ~= M
    lat = lat';
end
if size(Z,1) ~= M
    error('Error: length of first non-singleton dim of Z must equal length(lat)!');
end
if M < 3
    error('Error: length(lat) must be at least three');
end
%% convert latitude to colatitude
lat = 90-lat; 
%% calculate latitudinal grid distances
N = size(Z,2);
ra = earth_radius('m');% earth radius in meters
y = ra*deg2rad(lat);
dy = diff(y,1,1);
DY = repmat(dy,1,N);
%%
F = zeros(size(Z),class(Z));
dy1 = DY(1:end-1,:);
dy2 = DY(2:end,:);
alpha = (dy1+dy2)./dy1/6;
a0 = alpha.*(2*dy1-dy2);
a1 = alpha.*(dy1+dy2).^2./dy2;
a2 = alpha.*dy1./dy2.*(2*dy2-dy1);
%% first cumulative value
state0 = warning('query','MATLAB:nearlySingularMatrix');
state0 = state0.state;
warning('off','MATLAB:nearlySingularMatrix')
C = vander(y(1:3))\Z(1:3,:);
F(2,:) = C(1,:).*(y(2,:).^3-y(1,:).^3)/3 +...
    C(2,:).*(y(2,:).^2-y(1,:).^2)/2 +...
    C(3,:).*dy(1,:);
warning(state0,'MATLAB:nearlySingularMatrix')
%% other cumulative values
F(3:2:end,:) = cumsum(a0(1:2:end,:).*Z(1:2:M-2,:) +...
    a1(1:2:end,:).*Z(2:2:M-1,:) +...
    a2(1:2:end,:).*Z(3:2:M,:),1);
F(4:2:end,:) = cumsum(a0(2:2:end,:).*Z(2:2:M-2,:) +...
    a1(2:2:end,:).*Z(3:2:M-1,:) +...
    a2(2:2:end,:).*Z(4:2:M,:),1) +...
    repmat(F(2,:),ceil((M-3)/2),1);
return;
%%
function F = cumsimps_lon(lon,lat,Z)
%% function for numerical approximation of a continuous function f(x) 
%% over longitudes on a spherical surface such that z(x) = df/dx using 
%% Simpson-rule row-wise cumulative summation 
%% INPUT
%lon: array containing evenly or unevently spaced longitudinal grid points
%lat: array containing evenly or unevently spaced latitudinal grid points
%Z:   array or matrix each row of which represents the value of the 
%     integrand at the longitudinal grid points x 
%% OUTPUT
%F: array or matrix containing the approximation of the integral F(x) 
%% determine the size of y and check the number of rows for Z 
N = numel(lon);
M = numel(lat);
if size(lon,2) ~= N
    lon = lon';
end    
if size(lat,1) ~= M
    lat = lat';
end   
if size(Z,1) ~= M 
    error('Error: length of first non-singleton dim of Z must equal length(lat)');
end
if size(Z,2) ~= N
    error('Error: length of second non-singleton dim of Z must equal length(lon)'); 
end
if N < 3
    error('Error: length(lon) must be at least three');
end
%% convert latitude to colatitude
lat = 90-lat; 
%% calculate longitudinal grid distances
ra = earth_radius('m');% earth radius in meters
x = ra*deg2rad(lon);
X = repmat(x,M,1); % matrix of longitudinal grid points
WT = repmat(sind(lat),1,N); % matrix of latitudinal weight 
X = X.*WT; % weight by latitude
DX = diff(X,1,2);
%%
F = zeros(size(Z),class(Z));
dx1 = DX(:,1:end-1);
dx2 = DX(:,2:end);
alpha = (dx1+dx2)./dx1/6;
a0 = alpha.*(2*dx1-dx2);
a1 = alpha.*(dx1+dx2).^2./dx2;
a2 = alpha.*dx1./dx2.*(2*dx2-dx1);
%% first cumulative value
for i = 1:M
    state0 = warning('query','MATLAB:nearlySingularMatrix');
    state0 = state0.state;
    warning('off','MATLAB:nearlySingularMatrix')
    C = vander(X(i,1:3)')\Z(i,1:3)';
    F(i,2) = C(1)*(X(i,2)^3-X(i,1)^3)/3 +...
    C(2)*(X(i,2)^2-X(i,1)^2)/2 +...
    C(3)*DX(i,1);
    warning(state0,'MATLAB:nearlySingularMatrix')
end
%% other cumulative values
F(:,3:2:end) = cumsum(a0(:,1:2:end).*Z(:,1:2:N-2) +...
    a1(:,1:2:end).*Z(:,2:2:N-1) +...
    a2(:,1:2:end).*Z(:,3:2:N),2);
F(:,4:2:end) = cumsum(a0(:,2:2:end).*Z(:,2:2:N-2) +...
    a1(:,2:2:end).*Z(:,3:2:N-1) +...
    a2(:,2:2:end).*Z(:,4:2:N),2) +...
    repmat(F(:,2),1,ceil((N-3)/2));
return;