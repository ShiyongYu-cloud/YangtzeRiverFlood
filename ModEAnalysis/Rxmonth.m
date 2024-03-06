%% Calculating annual maximum consecutive x-month total 
%% precipitation (Rxmonth) and annual mean temperture 
%% from the production of ModE 
%% https://www.wdc-climate.de/ui/project?acronym=ModE
clear;
clc;
%% read ModE-RA data
lat = ncread('ModE-RA_ensmean_totprec_anom_wrt_1901-2000_1421-2008_mon.nc','latitude');
lon = ncread('ModE-RA_ensmean_totprec_anom_wrt_1901-2000_1421-2008_mon.nc','longitude');
time_ra = ncread('ModE-RA_ensmean_totprec_anom_wrt_1901-2000_1421-2008_mon.nc','time');
ptot_mean_ra = ncread('ModE-RA_ensmean_totprec_anom_wrt_1901-2000_1421-2008_mon.nc','totprec');
ptot_stdd_ra = ncread('ModE-RA_ensstd_totprec_anom_wrt_1901-2000_1421-2008_mon.nc','totprec');
temp_mean_ra = ncread('ModE-RA_ensmean_temp2_anom_wrt_1901-2000_1421-2008_mon.nc','temp2');
temp_stdd_ra = ncread('ModE-RA_ensstd_temp2_anom_wrt_1901-2000_1421-2008_mon.nc','temp2');
p_clim = ncread('ModE-RA_ensmean_totprec_clim_1901-2000_mon.nc','totprec');
t_clim = ncread('ModE-RA_ensmean_temp2_clim_1901-2000_mon.nc','temp2');
%% read ModE-RAclim data
ptot_mean_cl = ncread('ModE-RAclim_ensmean_totprec_anom_1421-2008_mon.nc','totprec');
%% read ModE-Sim data
time_sm1 = ncread('ModE-Sim_set_1420-3_ensmean_totprec_abs_1420-1849_mon.nc','time');
time_sm2 = ncread('ModE-Sim_set_1850-1_ensmean_totprec_abs_1850-2009_mon.nc','time');
ptot_mean_sm1 = ncread('ModE-Sim_set_1420-3_ensmean_totprec_abs_1420-1849_mon.nc','totprec');
ptot_mean_sm2 = ncread('ModE-Sim_set_1850-1_ensmean_totprec_abs_1850-2009_mon.nc','totprec');
%% transpose to lat-lon-time
ptot_mean_ra = permute(ptot_mean_ra,[2 1 3]);
ptot_stdd_ra = permute(ptot_stdd_ra,[2 1 3]);
temp_mean_ra = permute(temp_mean_ra,[2 1 3]);
temp_stdd_ra = permute(temp_stdd_ra,[2 1 3]);
ptot_mean_cl = permute(ptot_mean_cl,[2 1 3]);
ptot_mean_sm1 = permute(ptot_mean_sm1,[2 1 3]);
ptot_mean_sm2 = permute(ptot_mean_sm2,[2 1 3]);
p_clim = permute(p_clim,[2 1 3]);
t_clim = permute(t_clim,[2 1 3]);
%% convert lon from -180:180 to 0:360 
lon = [lon(97:end);180+(180+lon(1:96))]; 
[LON,LAT] = meshgrid(lon,lat);
days = [31 28 31 30 31 30 31 31 30 31 30 31];
L_ra = numel(time_ra)/12; %number of years
L_sm1 = numel(time_sm1)/12; %number of years
L_sm2 = numel(time_sm2)/12; %number of years
%L_sm = L_sm1 + L_sm2;
%% reshape data to lat-lon-month-year order
ptot_mean_ra = reshape(ptot_mean_ra,numel(lat),numel(lon),12,L_ra);
ptot_stdd_ra = reshape(ptot_stdd_ra,numel(lat),numel(lon),12,L_ra);
ptot_mean_cl = reshape(ptot_mean_cl,numel(lat),numel(lon),12,L_ra);
ptot_mean_sm1 = reshape(ptot_mean_sm1,numel(lat),numel(lon),12,L_sm1);
ptot_mean_sm2 = reshape(ptot_mean_sm2,numel(lat),numel(lon),12,L_sm2);
temp_mean_ra = reshape(temp_mean_ra,numel(lat),numel(lon),12,L_ra);
temp_stdd_ra = reshape(temp_stdd_ra,numel(lat),numel(lon),12,L_ra);
ptot_mean_sm = cat(4,ptot_mean_sm1(:,:,:,2:end),ptot_mean_sm2(:,:,:,1:end-1)); %combine ModE-Sim sub-datasets
%% re-arrange data along 0-360 longitude
p_clim = [p_clim(:,97:end,:) p_clim(:,1:96,:)];
t_clim = [t_clim(:,97:end,:) t_clim(:,1:96,:)];
ptot_mean_ra = [ptot_mean_ra(:,97:end,:,:) ptot_mean_ra(:,1:96,:,:)]; 
ptot_stdd_ra = [ptot_stdd_ra(:,97:end,:,:) ptot_stdd_ra(:,1:96,:,:)]; 
ptot_mean_cl = [ptot_mean_cl(:,97:end,:,:) ptot_mean_cl(:,1:96,:,:)]; 
ptot_mean_sm = [ptot_mean_sm(:,97:end,:,:) ptot_mean_sm(:,1:96,:,:)]; 
temp_mean_ra = [temp_mean_ra(:,97:end,:,:) temp_mean_ra(:,1:96,:,:)]; 
temp_stdd_ra = [temp_stdd_ra(:,97:end,:,:) temp_stdd_ra(:,1:96,:,:)]; 
%% add climatology back to ModE-RA and ModE-RAClim 
ptot_mean_ra = ptot_mean_ra + p_clim; 
ptot_mean_cl = ptot_mean_cl + p_clim; 
temp_mean_ra = temp_mean_ra + t_clim; 
%% calculate annual max precip and mean temp
%x = 1; %annual maximum monthly rainfall
x = 3; %consecutive 3-month total precipitaion
MAXRx_mean_ra = zeros(numel(lat),numel(lon),L_ra);
MAXRx_mean_cl = zeros(numel(lat),numel(lon),L_ra);
MAXRx_mean_sm = zeros(numel(lat),numel(lon),L_ra);
MAXRx_stdd_ra = zeros(numel(lat),numel(lon),L_ra);
t_mean_ra = zeros(numel(lat),numel(lon),L_ra);
t_stdd_ra = zeros(numel(lat),numel(lon),L_ra);
for i = 1:L_ra
    % convert to monthly total precipitation
    for j = 1:12
        ptot_mean_ra(:,:,j,i) = ptot_mean_ra(:,:,j,i)*days(j)*24*60*60;
        ptot_mean_cl(:,:,j,i) = ptot_mean_cl(:,:,j,i)*days(j)*24*60*60;
        ptot_mean_sm(:,:,j,i) = ptot_mean_sm(:,:,j,i)*days(j)*24*60*60;
        ptot_stdd_ra(:,:,j,i) = ptot_stdd_ra(:,:,j,i)*days(j)*24*60*60;
    end
    % calculate consecutive x-month total precipitaion
    rxmonth_mean_ra = movsum(ptot_mean_ra(:,:,:,i),x,3,'omitnan'); 
    rxmonth_mean_cl = movsum(ptot_mean_cl(:,:,:,i),x,3,'omitnan'); 
    rxmonth_mean_sm = movsum(ptot_mean_sm(:,:,:,i),x,3,'omitnan'); 
    rxmonth_stdd_ra = movsum(ptot_stdd_ra(:,:,:,i),x,3,'omitnan'); 
    % find annual max consecutive x-month total precip
    [max_rxmonth_mean_ra,id] = max(rxmonth_mean_ra,[],3); %find max along months
    max_rxmonth_mean_cl = max(rxmonth_mean_cl,[],3); %find max along months
    max_rxmonth_mean_sm = max(rxmonth_mean_sm,[],3); %find max along months
    max_rxmonth_stdd_ra = rxmonth_stdd_ra(id);
    MAXRx_mean_ra(:,:,i) = max_rxmonth_mean_ra;
    MAXRx_mean_cl(:,:,i) = max_rxmonth_mean_cl;
    MAXRx_mean_sm(:,:,i) = max_rxmonth_mean_sm;
    MAXRx_stdd_ra(:,:,i) = max_rxmonth_stdd_ra;
    % calculate annual mean temp
    t_mean_ra(:,:,i) = mean(temp_mean_ra(:,:,:,i),3);
    t_stdd_ra(:,:,i) = mean(temp_stdd_ra(:,:,:,i),3); 
end
years_ra = 1421:1:2008;
years_ra = years_ra';
%% extract and calculate area-weighted mean
T_mean_ra = zeros(L_ra,1);
T_stdd_ra = zeros(L_ra,1);
P_mean_ra = zeros(L_ra,1);
P_mean_cl = zeros(L_ra,1);
P_mean_sm = zeros(L_ra,1);
P_stdd_ra = zeros(L_ra,1);
mask = [110 28; 115 28; 115 32; 110 32; 110 28];
[lon_small,lat_small,t_mean_small_ra] = data_mask3D(t_mean_ra,lon,lat,mask);
[~,~,t_stdd_small_ra] = data_mask3D(t_stdd_ra,lon,lat,mask);
[~,~,MAXRx_mean_small_ra] = data_mask3D(MAXRx_mean_ra,lon,lat,mask);
[~,~,MAXRx_mean_small_cl] = data_mask3D(MAXRx_mean_cl,lon,lat,mask);
[~,~,MAXRx_mean_small_sm] = data_mask3D(MAXRx_mean_sm,lon,lat,mask);
[~,~,MAXRx_stdd_small_ra] = data_mask3D(MAXRx_stdd_ra,lon,lat,mask);
for i = 1:L_ra
    T_mean_ra(i) = area_mean(lon_small,lat_small,t_mean_small_ra(:,:,i));
    T_stdd_ra(i) = area_mean(lon_small,lat_small,t_stdd_small_ra(:,:,i));
    P_mean_ra(i) = area_mean(lon_small,lat_small,MAXRx_mean_small_ra(:,:,i));
    P_mean_cl(i) = area_mean(lon_small,lat_small,MAXRx_mean_small_cl(:,:,i));
    P_mean_sm(i) = area_mean(lon_small,lat_small,MAXRx_mean_small_sm(:,:,i));
    P_stdd_ra(i) = area_mean(lon_small,lat_small,MAXRx_stdd_small_ra(:,:,i));
end
T_ra = [(T_mean_ra+T_stdd_ra/2)-273.15;flipud((T_mean_ra-T_stdd_ra/2)-273.15);(T_mean_ra(1)+T_stdd_ra(1)/2)-273.15];
P_ra = [(P_mean_ra+P_stdd_ra);flipud(P_mean_ra-P_stdd_ra);(P_mean_ra(1)+P_stdd_ra(1))];
YR_ra = [years_ra;flipud(years_ra);years_ra(1)];
%standardize precipitation
id1 = years_ra<1850;
id2 = years_ra>=1850;
P1_mean_norm_ra = (P_mean_ra(id1) - mean(P_mean_ra(id1)))/std(P_mean_ra(id1));
P2_mean_norm_ra = (P_mean_ra(id2) - mean(P_mean_ra(id2)))/std(P_mean_ra(id2));
P_mean_norm_ra = [P1_mean_norm_ra;P2_mean_norm_ra];
P1_mean_norm_cl = (P_mean_cl(id1) - mean(P_mean_cl(id1)))/std(P_mean_cl(id1));
P2_mean_norm_cl = (P_mean_cl(id2) - mean(P_mean_cl(id2)))/std(P_mean_cl(id2));
P_mean_norm_cl = [P1_mean_norm_cl;P2_mean_norm_cl];
P1_mean_norm_sm = (P_mean_sm(id1) - mean(P_mean_sm(id1)))/std(P_mean_sm(id1));
P2_mean_norm_sm = (P_mean_sm(id2) - mean(P_mean_sm(id2)))/std(P_mean_sm(id2));
P_mean_norm_sm = [P1_mean_norm_sm;P2_mean_norm_sm];
%% plot results of sensitivity analysis
figure(1)
subplot(2,2,1) %plot ModE-RA
plot(years_ra,P_mean_ra);
hold on
plot(years_ra(id1),ones(numel(years_ra(id1)),1).*mean(P_mean_ra(id1)))
hold on
plot(years_ra(id2),ones(numel(years_ra(id2)),1).*mean(P_mean_ra(id2)))
xlim([1490 2010]);
ylim([600 1000]);
xlabel('Year (CE)');
ylabel('Rx3month (mm)');
grid on
title('ModE-RA');
subplot(2,2,2) %plot ModE-Sim
plot(years_ra,P_mean_sm);
hold on
plot(years_ra(id1),ones(numel(years_ra(id1)),1).*mean(P_mean_sm(id1)))
hold on
plot(years_ra(id2),ones(numel(years_ra(id2)),1).*mean(P_mean_sm(id2)))
xlim([1490 2010]);
ylim([700 1000]);
xlabel('Year (CE)');
ylabel('Rx3month (mm)');
grid on
title('ModE-Sim');
subplot(2,2,3) %plot ModE-RAClim
plot(years_ra,P_mean_cl);
hold on
plot(years_ra(id1),ones(numel(years_ra(id1)),1).*mean(P_mean_cl(id1)))
hold on
plot(years_ra(id2),ones(numel(years_ra(id2)),1).*mean(P_mean_cl(id2)))
xlim([1490 2010]);
ylim([500 1000]);
xlabel('Year (CE)');
ylabel('Rx3month (mm)');
grid on
title('ModE-RAClim');
subplot(2,2,4);
plot(years_ra,P_mean_norm_ra);
xlim([1490 2010]);
ylim([-3 3]);
xlabel('Year (CE)');
ylabel('Normalized Rx3month');
grid on
title('ModE-RA z-score');
%% identify precipitation-rich years
events = years_ra(P_mean_norm_ra>0.5);
%calculate frequency
[yrs,freq] = flood_freq(events,1480,2020,31);
%bootstrap
sim = 1;
L = 10;
nbsamples = 5000;
Xb = bootstrap_series(freq,sim,L,nbsamples); %bootstrap
CI = prctile(Xb',[2.5 25 50 75 97.5]); %calcualte confidence interval
CI = CI';
CI(CI<0) = 0;
CI95 = [CI(:,1);flipud(CI(:,5));CI(1,1)]; 
CI50 = [CI(:,2);flipud(CI(:,4));CI(1,2)]; 
YRS = [yrs;flipud(yrs);yrs(1)];
% plot frequency
figure(2)
fill(YRS,CI95,'r');
hold on
fill(YRS,CI50,'b');
hold on
plot(yrs,movmean(freq,31),'g');
xlim([1480 2020]);
xlabel('Year (CE)');
ylabel('Frequency');
% plotannual mean temperature
figure(3)
fill(YR_ra,T_ra,'blue');
hold on
plot(years_ra,T_mean_ra-273.1,'r');
xlim([1500 2010]);
xlabel('Year (CE)');
ylabel('Temperature ({\circ}C)');