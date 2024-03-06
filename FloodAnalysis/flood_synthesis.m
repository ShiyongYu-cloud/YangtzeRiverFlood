%Synthesizing individual sedimetnary record and bootstrapping
clc;
clear;
%% load sedimentary and instrumental flood data
data1 = 'sedimentary_flood_DGH.txt';
data2 = 'sedimentary_flood_DJH.txt';
data3 = 'sedimentary_flood_CTH.txt';
data4 = 'instrumental_flood_YC.txt';
h = fopen(data1);
DGH_data = textscan(h,'%f %f','headerlines',1,'delimiter','\t');
fclose(h);
h = fopen(data2);
DJH_data = textscan(h,'%f %f','headerlines',1,'delimiter','\t');
fclose(h);
h = fopen(data3);
CTH_data = textscan(h,'%f %f','headerlines',1,'delimiter','\t');
fclose(h);
INS_data = load(data4);
%
DGH.age = DGH_data{1};
DGH.err = DGH_data{2};
DJH.age = DJH_data{1};
DJH.err = DJH_data{2};
CTH.age = CTH_data{1};
CTH.err = CTH_data{2};
%% calculate flood frequency during the instrumental era
ins_min = min(INS_data);
DGH_ins = DGH.age(DGH.age >= ins_min);
DJH_ins = DJH.age(DJH.age >= ins_min);
CTH_ins = CTH.age(CTH.age >= ins_min);
a1 = 1320;
a2 = 1730;
a3 = 1760;
a4 = 1890;
b0 = 2020;
W = 31; % window size
[~,DGH_freq_ins] = flood_freq(DGH_ins,a4,b0,W); % Donggang
[~,DJH_freq_ins] = flood_freq(DJH_ins,a4,b0,W); % Dajing  
[~,CTH_freq_ins] = flood_freq(CTH_ins,a4,b0,W); % Chenta
[INS_yrs,INS_freq] = flood_freq(INS_data,a4,b0,W);% Instrumental 
%% build up empirical cdf
[X1_DGH,F1_DGH] = ecdf(DGH_freq_ins,0,0.4,100); % Donggang
[X1_DJH,F1_DJH] = ecdf(DJH_freq_ins,0,0.4,100); % Dajing
[X1_CTH,F1_CTH] = ecdf(CTH_freq_ins,0,0.4,100); % Chenta  
[X1_INS,F1_INS] = ecdf(INS_freq,0,0.4,100);     % Instrumental
% plot results
figure(1)
subplot(1,3,1); %Donggang
plot(X1_DGH,F1_DGH,'r',X1_INS,F1_INS,'b');
xlabel('Flood frequency');
ylabel('Cumulative probability');
xlim([0 0.4]);
xticks((0:0.1:0.4));
xtickformat('%.1f');
yticks((0:0.2:1));
ytickformat('%.1f');
legend
grid on
%
subplot(1,3,2); % Dajing
plot(X1_DJH,F1_DJH,'r',X1_INS,F1_INS,'b');
xlabel('Flood frequency');
ylabel('Cumulative probability');
xlim([0 0.4]);
xticks((0:0.1:0.4));
xtickformat('%.1f');
yticks((0:0.2:1));
ytickformat('%.1f');
legend
grid on
%
subplot(1,3,3); % Chenta
plot(X1_CTH,F1_CTH,'r',X1_INS,F1_INS,'b');
xlabel('Flood frequency');
ylabel('Cumulative probability');
xlim([0 0.4]);
xticks((0:0.1:0.4));
xtickformat('%.1f');
yticks((0:0.2:1));
ytickformat('%.1f');
legend
grid on
%% calibrate flood frequency against instrumental data 
sim = 1;
L = 10;
nbsamples = 5000;
% instrumental
Xb_INS = bootstrap_series(INS_freq,sim,L,nbsamples); %bootstrap
INS_CI = prctile(Xb_INS',[2.5 25 50 75 97.5]); %calcualte confidence interval
INS_CI = INS_CI';
INS_CI(INS_CI<0) = 0;
INS_CI95 = [INS_CI(:,1);flipud(INS_CI(:,5));INS_CI(1,1)]; 
INS_CI50 = [INS_CI(:,2);flipud(INS_CI(:,4));INS_CI(1,2)]; 
INS_YRS = [INS_yrs;flipud(INS_yrs);INS_yrs(1)];
% Donggang Lake
[DGH_yrs,DGH_freq] = flood_freq(DGH.age,a1,b0,W);
id = find(DGH_yrs >= min(INS_yrs));
F1_star = interp1(X1_DGH,F1_DGH,DGH_freq(id));
DGH_star = interp1(F1_INS,X1_INS,F1_star);
DGH_freq(id) = DGH_star;
Xb_DGH = bootstrap_series(DGH_freq,sim,L,nbsamples); %bootstrap
DGH_CI = prctile(Xb_DGH',[2.5 25 50 75 97.5]); %calcualte confidence interval
DGH_CI = DGH_CI';
DGH_CI(DGH_CI<0) = 0;
DGH_CI95 = [DGH_CI(:,1);flipud(DGH_CI(:,5));DGH_CI(1,1)]; 
DGH_CI50 = [DGH_CI(:,2);flipud(DGH_CI(:,4));DGH_CI(1,2)]; 
DGH_YRS = [DGH_yrs;flipud(DGH_yrs);DGH_yrs(1)];
% Dajing Lake
[DJH_yrs,DJH_freq] = flood_freq(DJH.age,a2,b0,W);
id = find(DJH_yrs >= min(INS_yrs));
F1_star = interp1(X1_DJH,F1_DJH,DJH_freq(id));
DJH_star = interp1(F1_INS,X1_INS,F1_star);
DJH_freq(id) = DJH_star;
Xb_DJH = bootstrap_series(DJH_freq,sim,L,nbsamples); %bootstrap
DJH_CI = prctile(Xb_DJH',[2.5 25 50 75 97.5]); %calcualte confidence interval
DJH_CI = DJH_CI';
DJH_CI(DJH_CI<0) = 0;
DJH_CI95 = [DJH_CI(:,1);flipud(DJH_CI(:,5));DJH_CI(1,1)]; 
DJH_CI50 = [DJH_CI(:,2);flipud(DJH_CI(:,4));DJH_CI(1,2)]; 
DJH_YRS = [DJH_yrs;flipud(DJH_yrs);DJH_yrs(1)]; 
% Chenta Lake
[CTH_yrs,CTH_freq] = flood_freq(CTH.age,a3,b0,W);
id = find(CTH_yrs >= min(INS_yrs));
F1_star = interp1(X1_CTH,F1_CTH,CTH_freq(id));
CTH_star = interp1(F1_INS,X1_INS,F1_star);
CTH_freq(id) = CTH_star;
Xb_CTH = bootstrap_series(CTH_freq,sim,L,nbsamples); %bootstrap
CTH_CI = prctile(Xb_CTH',[2.5 25 50 75 97.5]); %calcualte confidence interval
CTH_CI = CTH_CI';
CTH_CI(CTH_CI<0) = 0;
CTH_CI95 = [CTH_CI(:,1);flipud(CTH_CI(:,5));CTH_CI(1,1)]; 
CTH_CI50 = [CTH_CI(:,2);flipud(CTH_CI(:,4));CTH_CI(1,2)]; 
CTH_YRS = [CTH_yrs;flipud(CTH_yrs);CTH_yrs(1)]; 
% plot results
figure(2)
subplot(4,1,1); %Donggang Lake
fill(DGH_YRS,DGH_CI95,'b');
hold on
fill(DGH_YRS,DGH_CI50,'r');
hold on
plot(DGH_yrs,DGH_CI(:,3),'g');
ylabel('Flood frequency');
xlabel('Year (CE)');
xlim([1280 2020]);
%
subplot(4,1,2); % Dajing Lake
fill(DJH_YRS,DJH_CI95,'b');
hold on
fill(DJH_YRS,DJH_CI50,'r');
hold on
plot(DJH_yrs,DJH_CI(:,3),'g');
ylabel('Flood frequency');
xlabel('Year (CE)');
xlim([1720 2020]);
%
subplot(4,1,3); % Chenta Lake
fill(CTH_YRS,CTH_CI95,'b');
hold on
fill(CTH_YRS,CTH_CI50,'r');
hold on
plot(CTH_yrs,CTH_CI(:,3),'g');
ylabel('Flood frequency');
xlabel('Year (CE)');
xlim([1740 2020]);
subplot(4,1,4); % Instrumental
fill(INS_YRS,INS_CI95,'b');
hold on
fill(INS_YRS,INS_CI50,'r');
hold on
plot(INS_yrs,INS_CI(:,3),'g');
ylabel('Flood frequency');
xlabel('Year (CE)');
xlim([1880 2020]);
%% synthesize individual flood records
id1 = find(DGH_yrs<1780);
id2 = find(DGH_yrs>=1780 & DGH_yrs<1840);
id3 = find(DGH_yrs>=1840);
id4 = find(DJH_yrs>=1780 & DJH_yrs<1840); 
id5 = find(DJH_yrs>=1840); 
id6 = find(CTH_yrs>=1840); 
freq = zeros(length(DGH_freq),1);
freq(id1) = DGH_freq(id1)*1.0;
freq(id2) = 0.8*DGH_freq(id2)+0.2*DJH_freq(id4);
freq(id3) = 0.7*DGH_freq(id3)+0.2*DJH_freq(id5)+0.1*CTH_freq(id6);
Xb = bootstrap_series(freq,sim,L,nbsamples); %bootstrap
Xb_CI = prctile(Xb',[2.5 25 50 75 97.5]); %calcualte confidence interval
Xb_CI = Xb_CI';
Xb_CI(Xb_CI<0) = 0;
Xb_CI95 = [Xb_CI(:,1);flipud(Xb_CI(:,5));Xb_CI(1,1)]; 
Xb_CI50 = [Xb_CI(:,2);flipud(Xb_CI(:,4));Xb_CI(1,2)]; 
% plot results
figure(3)
fill(DGH_YRS,Xb_CI95,'b');
hold on
fill(DGH_YRS,Xb_CI50,'r');
hold on
plot(DGH_yrs,Xb_CI(:,3),'g');
ylabel('Flood frequency');
xlabel('Year (CE)');
xlim([1490 2010]);