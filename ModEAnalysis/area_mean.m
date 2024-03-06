function mean_D = area_mean(lon,lat,D)
%function for calculating area-weighted mean of data on a spherical surface
%note that grids can be unevenly spaced
%INPUT
%lon: vector containing longitudinal coordinates
%lat: vector containing latitudinal coordinates
%D:   2D matrix containing gridded data  
%OUTPUT
%mean_D: area-weighted value of D
%% check dimension
M = numel(lat);
N = numel(lon);
if size(D,1) == N && size(D,2) == M
    D = D';
end
%% calculate grid resolution
d_lon = diff(lon);
d_lat = diff(lat);
d_lon = abs(d_lon);
d_lat = abs(d_lat);
dy = zeros(M,1);
dx = zeros(1,N);
%% loop over for latitudinal grid distances
if M == 1
    dy(:,1) = mean(d_lon);
elseif M == 2
    dy(:,1) = d_lat;
else    
    for i = 1:M
         if i == 1
             dy(i) = d_lat(i); %reflective boundary condition
         elseif i > 1 && i < M
             dy(i) = d_lat(i-1)/2 + d_lat(i)/2;
         elseif i == M
             dy(i) = d_lat(i-1);%reflective boundary condition
         end
    end
end
%% loop over for longitudinal grid distances
if N == 1
    dx(1,:) = mean(d_lat);
elseif N == 2
    dx(1,:) = d_lon;
else    
    for j = 1:N
         if j == 1
             dx(j) = d_lon(j); %reflective boundary condition
         elseif j > 1 && j < N
             dx(j) = d_lon(j-1)/2 + d_lon(j)/2;
         elseif j == N
             dx(j) = d_lon(j-1); %reflective boundary condition
         end
    end  
end    
%% calculate latitude weight
wt = cosd(lat-dy/2);
WT = repmat(wt,1,N); 
%% calculate grid length 
ra = earth_radius('km');
dx = ra*deg2rad(dx);
dy = ra*deg2rad(dy);
%% build longitudinal and latitudinal grid matrix
DX = repmat(dx,M,1);
DX = DX.*WT; %weight by cosine of latitude
DY = repmat(dy,1,N);
%% calculate grid cell area
A = DX.*DY;
%% weight data by grid cell area
F = A/sum(A,'all'); %area weight of each grid
mean_D = sum(D.*F,'all','omitnan');
end