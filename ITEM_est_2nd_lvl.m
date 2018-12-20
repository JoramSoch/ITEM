function ITEM_est_2nd_lvl(SPM, mode)
% _
% Estimate Second-Level (Trial-Wise) Model
% FORMAT ITEM_est_2nd_lvl(SPM, mode)
%     SPM  - a structure specifying an estimated GLM
%     mode - a string indicating how to perform ReML estimation
% 
% FORMAT ITEM_est_2nd_lvl(SPM, mode) analyzes trial-wise parameter
% estimates in a second-level (trial-wise) GLM which uses a trial-wise
% design matrix with one matrix row per trial, plus some augmentation
% to handle nuisance conditions and regressors of no interest.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 22/11/2018, 12:45 (V0.1)
%  Last edit: 19/12/2018, 21:05 (V0.1)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get SPM.mat if necessary
%-------------------------------------------------------------------------%
if nargin == 0
    SPM_mat = spm_select(1,'^SPM\.mat$','Select SPM.mat!');
    SPM_dir = fileparts(SPM_mat); load(SPM_mat);
    SPM.swd = SPM_dir;
    ITEM_est_2nd_lvl(SPM);
    return
end;

% Set estimation mode if necessary
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(mode), mode = 'spm_reml'; end;

% Change to SPM.swd if specified
%-------------------------------------------------------------------------%
orig_dir = pwd;
try
    cd(SPM.swd);
catch
    SPM.swd = pwd;
end

% Get number of sessions
%-------------------------------------------------------------------------%
s = numel(SPM.Sess);

% Load GLM.mat in sub-directory
%-------------------------------------------------------------------------%
load(strcat(SPM.swd,'/','ITEM_est_1st_lvl','/','GLM1.mat'));

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_est_2nd_lvl: load');

% Load mask image
%-------------------------------------------------------------------------%
[M m_dim m_ind] = MA_load_mask(SPM);

% Load gamma estimates
%-------------------------------------------------------------------------%
G = ITEM_load_gammas(SPM, m_ind);
v = numel(m_ind);


%=========================================================================%
% E S T I M A T I O N                                                     %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_est_2nd_lvl: estimate');

% Preallocate beta estimates
%-------------------------------------------------------------------------%
B  = NaN(sum(GLM1.pr),prod(m_dim));
S2 = NaN(numel(GLM1.Sess),prod(m_dim));

% Estimate second-level model
%-------------------------------------------------------------------------%
for h = 1:numel(SPM.Sess)
    
    % Get data and design
    %---------------------------------------------------------------------%
    ih    = sum(GLM1.pr(1:h-1))+[1:GLM1.pr(h)];
    Yh    = G(GLM1.Sess(h).t,:);
    Xh    = GLM1.Sess(h).T;
    Qh{1} = eye(GLM1.tr(h));
    Qh{2} = GLM1.Sess(h).U;
    
    % Restricted maximum likelihood
    %---------------------------------------------------------------------%
    if strcmp(mode,'spm_reml')
       [Vh, s2] = ITEM_GLM_ReML(Yh, Xh, Qh{1}, Qh{2}, sprintf('ITEM_est_2nd_lvl: ReML estimation for session %d',h));
    end;
    if strcmp(mode,'spm_reml_sc')
        YY = (1/v)*(Yh*Yh');
       [Vh, s2] = spm_reml_sc(YY, Xh, Qh);
        Vh = full(Vh);
        s2 = full(s2)';
    end;
    
    % Maximum likelihood estimation
    %---------------------------------------------------------------------%
    [B(ih,m_ind), S2(h,m_ind)] = ITEM_GLM_MLE(Yh, Xh, Vh, sprintf('Analyze trial-wise response amplitudes from session %d',h));
    
    % Store (co)variance(s)
    %---------------------------------------------------------------------%
    GLM2.Sess(h).p = ih;
    GLM2.Sess(h).Q = Qh;
    GLM2.Sess(h).V = Vh;
    GLM2.Sess(h).s2= s2;
    clear Yh Xh Qh Vh s2
    
end;


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_est_2nd_lvl: save');
spm_progress_bar('Init',100,'Save condition-wise parameter estimates...','');

% Initialise image files
%-------------------------------------------------------------------------%
H = MA_init_header(SPM, false);
GLM2.swd = strcat(SPM.swd,'/','ITEM_est_2nd_lvl','/');
if ~exist(GLM2.swd,'dir'), mkdir(GLM2.swd); end;
cd(GLM2.swd);

% Save beta estimates
%-------------------------------------------------------------------------%
d = ceil(sum(GLM1.pr)/100);
for h = 1:numel(SPM.Sess)
    for k = 1:GLM1.pr(h)
        i = GLM2.Sess(h).p(k);
        H.fname   = strcat('beta_',MF_int2str0(i,4),'.nii');
        H.descrip = sprintf('ITEM_est_2nd_lvl: parameter estimate; session %d, condition %d', h, k);
        spm_write_vol(H,reshape(B(i,:),m_dim));
        GLM2.Vbeta(i) = H;
        if mod(i,d) == 0, spm_progress_bar('Set',(i/sum(GLM1.pr))*100); end;
    end;
end;

% Save sigma^2 estimates
%-------------------------------------------------------------------------%
for h = 1:numel(SPM.Sess)
    H.fname   = strcat('sigma_',MF_int2str0(h,4),'.nii');
    H.descrip = sprintf('ITEM_est_2nd_lvl: variance estimate; session %d', h);
    spm_write_vol(H,reshape(S2(h,:),m_dim));
    GLM2.Vsigma(h) = H;
    spm_progress_bar('Set',(h/s)*100);
end;

% Save GLM structure
%-------------------------------------------------------------------------%
save(strcat(GLM2.swd,'GLM2.mat'),'GLM2');

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);