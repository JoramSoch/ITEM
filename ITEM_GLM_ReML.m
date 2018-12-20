function [V, s2] = ITEM_GLM_ReML(Y, X, Q1, Q2, msg)
% _
% Restricted Maximum Likelihood for General Linear Model
% FORMAT [V, s2] = ITEM_GLM_ReML(Y, X, Q1, Q2, msg)
% 
%     Y   - an n x v data matrix of v time series with n data points
%     X   - an n x p design matrix of p regressors with n data points
%     Q1  - an n x n covariance matrix, the first variance component
%     Q2  - an n x n covariance matrix, the second variance component
%     msg - a string used as a message in the MATLAB command window
% 
%     V   - an n x n covariance matrix embodying covariance assumptions
%     s2  - a  1 x 2 vector of variance factors for the two components
% 
% FORMAT [V, s2] = ITEM_GLM_ReML(Y, X, Q1, Q2, msg) performs a restricted
% maximum likelihood (ReML) analysis to determine the mixing of covariance
% components Q1 and Q2 into the full covariance matrix V, given a general
% linear model with data Y and design matrix X.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 03/12/2018, 17:00 (V0.1)
%  Last edit: 04/12/2018, 15:00 (V0.1)


% Get model dimensions
%-------------------------------------------------------------------------%
v = size(Y,2);                  % number of time series
n = size(Y,1);                  % number of data points
p = size(X,2);                  % number of parameters

% Calculate sample covariance matrix
%-------------------------------------------------------------------------%
YY   = (1/v) * (Y*Y');
Q{1} = Q1;
Q{2} = Q2;

% Perform restricted maximum likelihood
%-------------------------------------------------------------------------%
if nargin < 5, msg = 'ReML estimation'; end;
fprintf('\n> %s:\n', msg);
[V, s2] = spm_reml(YY, X, Q);
fprintf('> ReML estimates: V = %3.3f x Q1 + %3.3f x Q2.\n', s2(1), s2(2));
if size(s2,1) > 1, s2 = s2'; end;