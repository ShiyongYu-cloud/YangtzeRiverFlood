function [lon_small,lat_small,data_small] = data_mask2D(data,lon,lat,mask)
%function for subsetting a large 2D data matrix
%INPUT 
%data_large: larger 2D data matrix
%lon: vector containing longitudinal grid points of the larger data matrix
%lat: vector containing latitudinal grid points of the larger data matrix
%mask: 4x2 matrix containing corner coordinates [lons,lats] of the smaller region
%OUTPUT
%lon_small: vector containing lontitudinal grid points of the masked data
%lat_small: veector containing latitudinal grid points of the masked data
%data_small: extracted 3D data matrix 
%%
%% read corner coordinates of the masked region
xmin = min(mask(:,1));
xmax = max(mask(:,1));
ymin = min(mask(:,2));
ymax = max(mask(:,2));
%% define mask id
lon_mask = lon>=xmin & lon<=xmax;
lat_mask = lat>=ymin & lat<=ymax;
%% extract longitudinal and latitudinal grid points of the mask
lon_small = lon(lon_mask);
lat_small = lat(lat_mask);
%% check dimension and extract data
M = numel(lat);
N = numel(lon);
if size(data,1) == M && size(data,2) == N       %lat-lon-time
    data_small = data(lat_mask,lon_mask);
elseif size(data,1) == N && size(data,2) == M   %lon-lat-time
    data_small = data(lon_mask,lat_mask);
end
end