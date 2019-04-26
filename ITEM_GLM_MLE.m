function [beta_est, sig2_est] = ITEM_GLM_MLE(Y, X, V, msg)
% _
% Maximum Likelihood Estimation of General Linear Model
% FORMAT [beta_est, sig2_est] = ITEM_GLM_MLE(Y, X, V, msg)
% 
%     Y   - an n x v data matrix of v time series with n data points
%     X   - an n x p design matrix of p regressors with n data points
%     V   - an n x n covariance matrix embodying covariance assumptions
%     msg - a string used as a message on the SPM progress bar
% 
%     beta_est - a p x v matrix (estimated regression coefficients)
%     sig2_est - a 1 x v vector (estimated residual variance)
% 
% FORMAT [beta_est, sig2_est] = ITEM_GLM_MLE(Y, X, V) returns the maximum
% likelihood parameter estimates for a general linear model with data Y,
% design matrix X, covariance matrix V.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 29/10/2014, 14:05 (V0.0)
%  Last edit: 26/04/2019, 17:35 (V0.2)


% Get model dimensions
%-------------------------------------------------------------------------%
v = size(Y,2);                  % number of time series
n = size(Y,1);                  % number of data points
p = size(X,2);                  % number of parameters
d = floor(v/100);

% Init progress bar
%-------------------------------------------------------------------------%
if nargin < 4, msg = 'Estimate GLM ...'; end;
Finter = spm('FigName','ITEM_GLM_MLE: estimate');
spm_progress_bar('Init', 100, msg, '');

% Prepare parameter estimation
%-------------------------------------------------------------------------%
P = inv(V);                     % inverse of covariance matrix

% Calculate parameter estimates
%-------------------------------------------------------------------------%
beta_est = (X'*P*X)^-1 * X'*P*Y;% estimated regression coefficients
resi_est = Y - X*beta_est;      % estimated residuals/errors/noise
sig2_est = zeros(1,v);          % estimated residual variance
for j = 1:v
    sig2_est(j) = 1/n * (resi_est(:,j))' * P * (resi_est(:,j));
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;
% An implementation without a loop would look like this:
%   sig2_est = 1/n * sum(resi_est.^2), if V = I
% However, it would be impossible to make a progress bar then.

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');