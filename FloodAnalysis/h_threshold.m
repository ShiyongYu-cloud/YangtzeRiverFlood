function h_opt = h_threshold(ins,sed)
%function for finding the optimal discrimination threshold
%INPUT
%sed_flood: structure containing sedimetnary flood events
%ins: structure containing instrumental flood events
%OUTPUT
%h_opt: optimal estimate of the discrimination threshold
%% calculate instrumental flood frequency
years = ins.years;    
ins_event = zeros(length(years),1);
id = strcmpi(ins.labels,'flood');
ins_event(id) = 1;
L = 31;
win = gausswin(L);
ins_freq = conv(ins_event,win,'same')/L;
%% calcualte sedimentary flood frequency
N = 500;
h = linspace(0,1,N);
COR_F = zeros(N,1);
for j = 1:N
    %find flood events
    sed_flood = find_flood(sed,h(j));
    sed_event = zeros(length(years),1);
       for i = 1:length(years)
            if ismember(years(i),sed_flood.age_mean)
                sed_event(i) = 1;
            else
                sed_event(i) = 0;
            end
       end    
    sed_freq = conv(sed_event,win,'same')/L;
    COR_F(j) = corr(sed_freq,ins_freq,'rows','complete');
end
id = find(h>0.05);
h = h(id);
COR_F = COR_F(id);
id = find(COR_F == max(COR_F(2:end,1)));
h_opt = h(id(end)); %pick the largest one
end