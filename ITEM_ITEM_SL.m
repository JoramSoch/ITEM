function Yp = ITEM_ITEM_SL(Y1, X1, V1, X2, V2, SLs, msg)
% _
% Searchlight-Based Estimation of Inverse Transformed Encoding Model
% FORMAT Yp = ITEM_ITEM_SL(Y1, X1, V1, X2, V2, SLs, msg)
% 
%     Y1  - an n1 x q  training data matrix of q variables with n data points
%     X1  - an n1 x v  training design matrix of v instances with n data points
%     V1  - an n1 x n1 training covariance matrix embodying correlation assumptions
%     X2  - an n2 x v  test design matrix of v instances with n data points
%     V2  - an n2 x n2 test covariance matrix embodying correlation assumptions
%     SLs - a   v x 1  cell array containing searchlight indices
%     msg - a string used as a message on the SPM progress bar
% 
%     Yp  - an n2 x q x v array of predictions for test data in all instances
% 
% FORMAT Yp = ITEM_ITEM_SL(Y1, X1, V1, X2, V2, SLs, msg) performs an
% ITEM-based searchlight analysis (or searchlight-based ITEM analysis)
% to decode out-of-sample data (*2) by learning from in-sample data (*1)
% and returning predictions for out-of-sample data points (Yp).
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 10/05/2019, 12:55 (V0.2)
%  Last edit: 30/11/2021, 14:48 (V0.3)


% Get model dimensions
%-------------------------------------------------------------------------%
n1 = size(X1,1);                % number of data points (training)
n2 = size(X2,1);                % number of data points (test)
q  = size(Y1,2);                % number of variables
v  = size(X1,2);                % number of instances
d  = floor(v/100);

% Init progress bar
%-------------------------------------------------------------------------%
if nargin < 7, msg = 'Searchlight-based ITEM analysis...'; end;
Finter = spm('FigName','ITEM_ITEM_SL: estimate');
spm_progress_bar('Init', 100, msg, '');

% Prepare parameter estimation
%-------------------------------------------------------------------------%
P1 = inv(V1);                   % precision matrix                  [n x n]
W1 = sqrtm(P1);                 % whitening matrix (training)       [n x n]
W2 = sqrtm(inv(V2));            % whitening matrix (test)           [n x n]

% Calculate auxiliary variables
%-------------------------------------------------------------------------%
X1  =[X1, ones(n1,1)];          % augmented design matrix           [n x v]
X2  =[X2, ones(n2,1)];          % augmented design matrix           [n x v]
PX  = P1 * X1;                  % precision times design            [n x v]
XPY = PX'* Y1; % = X1'* P1 * Y1 % design times data (whitened)      [v x q]
WX  = W2 * X2;                  % whitened design matrix (test)     [n x v]

% Perform decoding analysis
%-------------------------------------------------------------------------%
Yp = zeros(n2,q,v);
for j = 1:v
    vj = [SLs{j}, (v+1)];       % indices of current voxel's searchlight
    B1 = inv(X1(:,vj)'*PX(:,vj)) * XPY(vj,:);   % plus implicit baseline
    Yp(:,:,j) = WX(:,vj) * B1;  % generate predictions in test data set
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;
clear PX XPY WX
clear vj B1

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');