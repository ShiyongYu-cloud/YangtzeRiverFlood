clc;
clear;
close all
%% process the instrumental data
instrumental_record = 'Yichang_QZmax.txt';
flood_stage = 52;  %local flood stage
critical_h = 2;    %critical height above the flood stage
core_year = 2019;  %year sediment core was collected
%process instrumental data
ins = ins_process(instrumental_record,flood_stage,critical_h,core_year);
%lablel flood years in the instrumental record
ins.labels = cell(length(ins.years),1);
for i = 1:length(ins.years)
    if ismember(ins.years(i),ins.flood_years) 
        ins.labels{i} = 'flood';
    else
        ins.labels{i} = 'non_flood';    
    end
end
%plot results
figure(1)
plot(ins.years,ins.Zmax);
hold on
plot(ins.flood_years,ins.flood_Zmax,'o');
hold on
plot([ins.years(1) ins.years(end)],[flood_stage+critical_h flood_stage+critical_h],'r');
xlabel('Year (CE)');
ylabel('Peak stage (m)');
%%
%% process the sedimentary data 
sed_data = 'DGH_EM3.txt';
age_model = 'DGH_ages.txt';
%process sedimentry data
low_pass = 81; % 41_point moving min filter
high_pass = 21; % 11-point moving mean filter
sed = sed_process(sed_data,age_model,low_pass,high_pass);
%plot result
figure(2)
subplot(3,1,1); plot(sed.age_mean,sed.EM,sed.age_mean,sed.EM_trend);
xlim([min(sed.age_mean) max(sed.age_mean)]);
ylabel('EM fraction');
set(gca,'Xdir','reverse');
subplot(3,1,2); plot(sed.age_mean,sed.EM_detrend,sed.age_mean,sed.EM_mean);
xlim([min(sed.age_mean) max(sed.age_mean)]);
xlabel('Calendar age (BP)');
ylabel('detrended EM');
set(gca,'Xdir','reverse');
subplot(3,1,3); bar(sed.age_mean,sed.EM_score);
xlim([min(sed.age_mean) max(sed.age_mean)]);
xlabel('Calendar age (BP)');
ylabel('EM score');
set(gca,'Xdir','reverse');
%% validate the sedimentary record against instrumental data
[TPR,FPR,AUC] = my_roc(ins,sed);
figure(3)
plot(FPR,TPR);
hold on
plot([0 1],[0 1],'r-');
grid on;
axis square
xticks([0 0.2 0.4 0.6 0.8 1.0]);
yticks([0 0.2 0.4 0.6 0.8 1.0]);
xlabel('False positive rate'); 
ylabel('True positive rate');
%% find the optimal threshold on EM score
h_opt = h_threshold(ins,sed);
%% identify flood events from the sedimentary record
sed_flood = find_flood(sed,h_opt);
figure(4)
bar(1950-sed.age_mean,sed.EM_score);
hold on
plot(sed_flood.age_mean,sed_flood.EM_score,'r.');
hold on
line([1950-sed.age_mean(1) 1950-sed.age_mean(end)],[h_opt h_opt],'Color','green');
ylim([0 1]);
xlabel('Calendar age (BP)');
ylabel('EM score');
%% regress EM score on Qmax
%extract data during the instrumental period
[x_em,y_qm] = EM_Qmax(sed_flood,ins);
%bootstrap the linear regression model
[PAR,par] = reg_bootstrap(x_em,y_qm);
figure(5)
tiledlayout(2,3);
nexttile([2 2])
plot(x_em,y_qm,'o')
xlabel('EM score');
ylabel('Peak Discharge (m^3/s)');
axis square;
nexttile
histogram(PAR(:,1))
xlabel('Intercept');
nexttile
histogram(PAR(:,2))
xlabel('Slope');
%% predict Qmax from the EM score of sedimentary flood events
X = [ones(size(sed_flood.EM_score)) sed_flood.EM_score];
QMAX = X*PAR';
sed_flood.Qmax_mean = mean(QMAX,2);
sed_flood.Qmax_mean = round(sed_flood.Qmax_mean/100)*100;
sed_flood.Qmax_prct = prctile(QMAX',[2.5 50 97.5])';
sed_flood.Qmax_prct = round(sed_flood.Qmax_prct/100)*100;
figure(6)
x = sed_flood.age_mean;
y = sed_flood.Qmax_mean;
xmins = sed_flood.age_mean - sed_flood.age_min;
xplus = sed_flood.age_max - sed_flood.age_mean;
ymins = sed_flood.Qmax_mean - sed_flood.Qmax_prct(:,1);
yplus = sed_flood.Qmax_prct(:,3) - sed_flood.Qmax_mean;
errorbar(x,y,ymins,yplus,xmins,xplus,'sq','MarkerFaceColor','red');
grid on;
xlabel('Calendar age (CE)');
ylabel('Peak discharge (m^3/s)');