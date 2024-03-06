function flood = find_flood(sed,h)
%function for finding flood events given a critical EM score
%sed: structure containing sedimentary data
%h: critical EM score
%flood: structure containing flood events
%%
%initial search by finding peaks
[pks,locs] = findpeaks(sed.EM_score,sed.age_mean,'MinPeakDistance',1);
id1 = ismember(sed.age_mean,locs);
depth = sed.depth(id1);
EM_score = sed.EM_score(id1);
age_min = sed.age_min(id1);
age_max = sed.age_max(id1);
age_medn = sed.age_medn(id1);
age_mean = sed.age_mean(id1);
%refine the results by settitng peaks above a small positive number
id2 = find(pks>h);
depth = depth(id2);
EM_score = EM_score(id2);
age_min = age_min(id2);
age_max = age_max(id2);
age_mean = age_mean(id2);
age_medn = age_medn(id2);
%refine the result again by choosing unique ages
[~,ia,~] = unique(round(age_mean),'stable');
flood.depth = depth(ia);
flood.EM_score = EM_score(ia);
flood.age_min = age_min(ia);
flood.age_max = age_max(ia);
flood.age_mean = age_mean(ia);
flood.age_medn = age_medn(ia);
%round and convert ages to AD/BC scale
flood.age_min = 1950-round(flood.age_min);
flood.age_max = 1950-round(flood.age_max);
flood.age_mean = 1950-round(flood.age_mean);
flood.age_medn = 1950-round(flood.age_medn);
%swap min and max ages 
temp = flood.age_min;
flood.age_min = flood.age_max;
flood.age_max = temp;
end