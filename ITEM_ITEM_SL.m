function DA = ITEM_ITEM_SL(Y_in, X_in, V_in, Y_out, X_out, V_out, v2v, vXv, type, msg)
% _
% Searchlight-Based Estimation of Inverse Transformed Encoding Model
% FORMAT DA = ITEM_ITEM_SL(Y_in, X_in, V_in, Y_out, X_out, V_out, v2v, vXv, type, msg)
% 
%     Y_in  - an n_in x q training data matrix of q variables with n data points
%     X_in  - an n_in x v training design matrix of v instances with n data points
%     V_in  - an n_in x n_in training covariance matrix embodying correlation assumptions
%     Y_out - an n_out x q test data matrix of q variables with n data points
%     X_out - an n_out x v test design matrix of v instances with n data points
%     V_out - an n_out x n_out test covariance matrix embodying correlation assumptions
%     v2v   - a  v x v logical matrix indicating searchlight assignments
%     vXv   - a  v x v logical matrix prescribing covariance analysis
%     type  - a string naming the type of analysis ('class' or 'recon')
%     msg   - a string used as a message on the SPM progress bar
% 
%     DA    - a  1 x v matrix of decoding accuracies (classification) or
%             a  q x v matrix of correlation coefficients (reconstruction)
% 
% FORMAT DA = ITEM_ITEM_SL(Y_in, X_in, V_in, Y_out, X_out, V_out, v2v, vXv, type, msg)
% performs an ITEM-based searchlight analysis (or searchlight-based ITEM analysis)
% to decode out-of-sample data (*_out) by learning from in-sample data (*_in)
% with searchlights (v2v/vXv) using classification or reconstruction (type).
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 10/05/2019, 12:55 (V0.2)
%  Last edit: 10/05/2019, 15:10 (V0.2)


% Get model dimensions
%-------------------------------------------------------------------------%
q = size(Y_in,2);               % number of variables
v = size(X_in,2);               % number of instances
d = floor(v/100);

% Init progress bar
%-------------------------------------------------------------------------%
if nargin < 10, msg = 'Searchlight-based ITEM analysis...'; end;
spm_progress_bar('Init', 100, msg, '');

% Prepare parameter estimation
%-------------------------------------------------------------------------%
P_in  = inv(V_in);              % precision matrix                  [n x n]
W_in  = sqrt(P_in);             % whitening matrix (training)       [n x n]
W_out = sqrtm(inv(V_out));      % whitening matrix (test)           [n x n]

% Calculate auxiliary variables
%-------------------------------------------------------------------------%
PX  = P_in * X_in;              % precision times design            [n x v]
XPY = X_in'* P_in * Y_in;       % design times data (whitened)      [v x q]
WX  = W_out* X_out;             % whitened design matrix (test)     [n x v]

% Estimate design covariance
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_ITEM_SL: estimate (1)');
qq = sum(sum(vXv));             % number of all voxel-to-voxel combinations
pi = zeros(1,qq);               % row indices
pj = zeros(1,qq);               % col indices
pn = zeros(1,qq);               % XPX values
pc = 0;                         % qq counter
for j = 1:v
    vj     = find(vXv(j,:));    % indices of current voxel's co-voxels
    qj     = numel(vj);         % number of voxels
    qi     = pc + [1:qj];       % vector indices
    pc     = pc + qj;           % number of entries
    pi(qi) = j;
    pj(qi) = vj;
    pn(qi) = X_in(:,j)'*PX(:,vj);
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;
XPX = sparse(pi,pj,pn);         % design covariance matrix          [v x v]
clear qq pi pj pn pc vj qj qi
clear PX

% Prepare decoding analysis
%-------------------------------------------------------------------------%
if strcmp(type,'class')         % classification between classes
    i = find(sum(Y_out,2)>0)';  % trials of interest, as determined by contrast
    Y_true = Y_out(i,:);        % trial-wise true conditions to be decoded
    DA = zeros(1,v);
end;
if strcmp(type,'recon')         % reconstruction of modulators
    i = cell(q,1);              % trials of interest, as determined by contrast
    Y_true = cell(1,q);         % trial-wise true variables to be decoded
    for k = 1:q
        i{k} = find(Y_out(:,k)~=0)';
        Y_true{k} = Y_out(i{k},k);
    end;
    DA = zeros(q,v);
end;

% Perform decoding analysis
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_ITEM_SL: estimate (2)');
for j = 1:v
    vj     = find(v2v(j,:));    % indices of current voxel's searchlight
    B_in   = inv(XPX(vj,vj)) * XPY(vj,:);
    Y_pred = WX(:,vj) * B_in;
    if strcmp(type,'class')     % classification: proportion correct
        Y_class = Y_pred(i,:)==repmat(max(Y_pred(i,:),[],2),[1 q]);
        DA(j)   = (1/numel(i)) * sum(sum(Y_true.*Y_class,2),1);
    end;
    if strcmp(type,'recon')     % reconstruction: correlation coefficient
        for k = 1:q
            DA(k,j) = corr(Y_true{k},Y_pred(i{k},k));
        end;
    end;
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;
clear XPX XPY
clear Y_true Y_pred Y_class

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');