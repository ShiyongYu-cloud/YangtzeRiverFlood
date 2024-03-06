function [years,freq] = flood_freq(data,A,B,W)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
data = data(:);
years = A:1:B;
years = years';
win = gausswin(W);
%win = ones(W,1);
flood_id = ismember(years,data);
f1ood_event = double(flood_id);
freq = conv(f1ood_event,win,'same');
freq = freq/W;
end