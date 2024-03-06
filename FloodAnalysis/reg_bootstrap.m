function [PAR,par] = reg_bootstrap(x,y)
%function for bootstrapping a linear regression model
%INPUT
%x: independent variable
%y: dependent variable
%OUTPUT
%PAR: all bootstrapped parameters
%par: structure containing bootstrapped mean and 95% ci of parameters
%%
%% Perform a linear regression, and compute the residuals.
y = y(:);
if size(x,1) ~= size(y,1)
    x = x';
end    
X = [ones(size(y)) x];
b = regress(y,X);
yfit = X*b;
residu = y - yfit;
%% bootstrap the residuals
PAR = bootstrp(1000,@(bootrap)regress(yfit+bootrap,X),residu);
%% calculate mean and 95% CI
par.mean = mean(PAR);
par.prct = prctile(PAR,[2.5 50 97.5]);
end