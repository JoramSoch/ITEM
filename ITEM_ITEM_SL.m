function DA = ITEM_ITEM_SL(Y1, X1, V1, Y2, X2, V2, SLs, type, msg)
% _
% Searchlight-Based Estimation of Inverse Transformed Encoding Model
% FORMAT DA = ITEM_ITEM_SL(Y1, X1, V1, Y2, X2, V2, v2v, vXv, type, msg)
% 
%     Y1   - an n1 x q  training data matrix of q variables with n data points
%     X1   - an n1 x v  training design matrix of v instances with n data points
%     V1   - an n1 x n1 training covariance matrix embodying correlation assumptions
%     Y2   - an n2 x q  test data matrix of q variables with n data points
%     X2   - an n2 x v  test design matrix of v instances with n data points
%     V2   - an n2 x n2 test covariance matrix embodying correlation assumptions
%     SLs  - a   v x 1  cell array containing searchlight indices
%     type - a string naming the type of analysis ('class' or 'recon')
%     msg  - a string used as a message on the SPM progress bar
% 
%     DA   - a 1 x v matrix of decoding accuracies (classification) or
%            a q x v matrix of correlation coefficients (reconstruction)
% 
% FORMAT DA = ITEM_ITEM_SL(Y1, X1, V1, Y2, X2, V2, SLs, type, msg) performs
% an ITEM-based searchlight analysis (or searchlight-based ITEM analysis)
% to decode out-of-sample data (*2) by learning from in-sample data (*1)
% with searchlights (SLs) using classification or reconstruction (type).
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 10/05/2019, 12:55 (V0.2)
%  Last edit: 24/11/2021, 17:36 (V0.3)


% Get model dimensions
%-------------------------------------------------------------------------%
q = size(Y1,2);                 % number of variables
v = size(X1,2);                 % number of instances
d = floor(v/100);

% Init progress bar
%-------------------------------------------------------------------------%
if nargin < 10, msg = 'Searchlight-based ITEM analysis...'; end;
Finter = spm('FigName','ITEM_ITEM_SL: estimate');
spm_progress_bar('Init', 100, msg, '');

% Prepare parameter estimation
%-------------------------------------------------------------------------%
P1 = inv(V1);                   % precision matrix                  [n x n]
W1 = sqrtm(P1);                 % whitening matrix (training)       [n x n]
W2 = sqrtm(inv(V2));            % whitening matrix (test)           [n x n]

% Calculate auxiliary variables
%-------------------------------------------------------------------------%
PX  = P1 * X1;                  % precision times design            [n x v]
XPY = PX'* Y1; % = X1'* P1 * Y1 % design times data (whitened)      [v x q]
WX  = W2 * X2;                  % whitened design matrix (test)     [n x v]

% Prepare decoding analysis
%-------------------------------------------------------------------------%
if strcmp(type,'class')         % classification between classes
    i  = find(sum(Y2,2)>0)';    % trials of interest, as determined by contrast
    Yt = Y2(i,:);               % trial-wise true conditions to be decoded
    DA = zeros(1,v);
end;
if strcmp(type,'recon')         % reconstruction of modulators
    i  = cell(1,q);             % trials of interest, as determined by contrast
    Yt = cell(1,q);             % trial-wise true variables to be decoded
    for k = 1:q
        i{k}  = find(Y2(:,k)~=0);
        Yt{k} = Y2(i{k},k);
    end;
    DA = zeros(q,v);
end;

% Perform decoding analysis
%-------------------------------------------------------------------------%
for j = 1:v
    vj = SLs{j};                % indices of current voxel's searchlight
    B1 = inv(X1(:,vj)'*PX(:,vj)) * XPY(vj,:);
    Yp = WX(:,vj) * B1;
    if strcmp(type,'class')     % classification: proportion correct
        Yc    = Yp(i,:)==repmat(max(Yp(i,:),[],2),[1 q]);
        DA(j) = (1/numel(i)) * sum(sum(Yt.*Yc,2),1);
    end;
    if strcmp(type,'recon')     % reconstruction: correlation coefficient
        for k = 1:q
            DA(k,j) = corr(Yt{k},Yp(i{k},k));
        end;
    end;
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;
clear PX XPY WX
clear Yt Yp Yc

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');