function G = ITEM_load_gammas(SPM, m_ind)
% _
% Load Trial-Wise Parameter Estimates from First-Level Model
% FORMAT G = ITEM_load_gammas(SPM, m_ind)
% 
%     SPM   - a structure specifying an estimated GLM
%     m_ind - a 1 x v vector indexing in-mask voxels
% 
%     G     - a t x v parameter matrix (t: trials; v: voxels)
% 
% FORMAT G = ITEM_load_gammas(SPM, m_ind) loads only in-mask trial-wise
% estimates from a first-level (scan-wise) and returns a parameter matrix.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 08/08/2017, 15:10 (V0.0)
%  Last edit: 05/12/2018, 11:45 (V0.1)


% Load mask if necessary
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(m_ind)
    [M m_dim m_ind] = MA_load_mask(SPM);
end;
clear M m_dim

% Get data dimensions
%-------------------------------------------------------------------------%
load(strcat(SPM.swd,'/','ITEM_est_1st_lvl','/','GLM1.mat'));
t = sum(GLM1.tr);
v = numel(m_ind);
d = ceil(t/100);

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_load_gammas: load');
spm_progress_bar('Init',100,'Load trial-wise parameter estimates...','');

% Load gamma estimates
%-------------------------------------------------------------------------%
G = zeros(t,v);
for h = 1:numel(SPM.Sess)
    for k = 1:GLM1.tr(h)
        i = GLM1.Sess(h).t(k);
        g_str  = strcat(GLM1.swd,'/',GLM1.Vgamma(i).fname);
        g_hdr  = spm_vol(g_str);
        g_img  = spm_read_vols(g_hdr);
        G(i,:) = g_img(m_ind);
        if mod(i,d) == 0, spm_progress_bar('Set',(i/t)*100); end;
    end;
end;
clear g_str g_hdr g_img

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');