function ITEM_est_1st_lvl(SPM, mode, nind)
% _
% Estimate First-Level (Scan-Wise) Model
% FORMAT ITEM_est_1st_lvl(SPM, mode, nind)
%     SPM  - a structure specifying an estimated GLM
%     mode - a string indicating how to handle filter regressors
%     nind - a vector indexing non-trial-wise nuisance conditions
% 
% FORMAT ITEM_est_1st_lvl(SPM, mode, nind) calculates trial-wise parameter
% estimates for a first-level (scan-wise) GLM which uses a trial-wise
% design matrix with one HRF regressor per trial, plus some augmentation
% to handle nuisance conditions and regressors of no interest.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 22/11/2018, 08:20 (V0.1)
%  Last edit: 20/12/2018, 06:00 (V0.1)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get SPM.mat if necessary
%-------------------------------------------------------------------------%
if nargin == 0
    SPM_mat = spm_select(1,'^SPM\.mat$','Select SPM.mat!');
    SPM_dir = fileparts(SPM_mat); load(SPM_mat);
    SPM.swd = SPM_dir;
    ITEM_est_1st_lvl(SPM);
    return
end;

% Estimate model if necessary
%-------------------------------------------------------------------------%
if ~isfield(SPM.xVi,'V')
    SPM_mat = strcat(SPM.swd,'/','SPM.mat');
    MA_GLM_AR_only(SPM_mat); load(SPM_mat);
    ITEM_est_1st_lvl(SPM);
    return
end;

% Set estimation mode if necessary
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(mode), mode = 'DCT'; end;

% Set nuisance conditions if necessary
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(nind), nind = []; end;

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

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_est_1st_lvl: load');

% Load mask image
%-------------------------------------------------------------------------%
[M m_dim m_ind] = MA_load_mask(SPM);

% Load time series
%-------------------------------------------------------------------------%
Y = MA_load_data(SPM, m_ind);
v = numel(m_ind);


%=========================================================================%
% E S T I M A T I O N   ( 1 ) :   T R I A L - W I S E   D E S I G N       %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_est_1st_lvl: estimate (1)');

% Preallocate GLM structure
%-------------------------------------------------------------------------%
GLM1.n = zeros(1,s);            % number of scans per session
GLM1.p = zeros(1,s);            % number of regressors of interest
GLM1.r = zeros(1,s);            % number of regressors of no interest
GLM1.t = zeros(1,s);            % number of trials per session

% Convert standard to first-level model
%-------------------------------------------------------------------------%
for h = 1:numel(SPM.Sess)
    
    % Preallocate onsets and durations
    %---------------------------------------------------------------------%
    ind = [];
    ons = [];
    dur = [];
    T   = [];
    p   = 0;
    c   = 0;
    n   = numel(SPM.Sess(h).row);
    if ~isempty(nind)
        names0     = cell(1,numel(nind));
        onsets0    = cell(1,numel(nind));
        durations0 = cell(1,numel(nind));
        pmod0      = struct([]);
    end;
    
    % Collect trials from conditions
    %---------------------------------------------------------------------%
    for i = 1:numel(SPM.Sess(h).U)
        if ~ismember(i,nind)
            ind = [ind; i*ones(numel(SPM.Sess(h).U(i).ons),1)];
            ons = [ons; SPM.Sess(h).U(i).ons];
            dur = [dur; SPM.Sess(h).U(i).dur];
            T_i = ones(numel(SPM.Sess(h).U(i).ons),1);
            if ~strcmp(SPM.Sess(h).U(i).P(1).name,'none')
                for j = 1:numel(SPM.Sess(h).U(i).P)
                    P_j = SPM.Sess(h).U(i).P(j).P;      % get modulator
                    P_j = P_j - mean(P_j);              % mean-centering
                    P_j = P_j./ max(abs(P_j));          % normalization
                    T_i = [T_i, P_j];
                end;
            end;
            T = blkdiag(T, T_i);
            p = p + size(T_i,2);
        else % non-trial-wise nuisance condition
            c = c + 1;
            names0{c}     = strcat('N',num2str(c));
            onsets0{c}    = SPM.Sess(h).U(i).ons;
            durations0{c} = SPM.Sess(h).U(i).dur;
            if ~strcmp(SPM.Sess(h).U(i).P(1).name,'none')
                for j = 1:numel(SPM.Sess(h).U(i).P)
                    pmod0(s).name{j}  = SPM.Sess(h).U(i).P(j).name;
                    pmod0(s).param{j} = SPM.Sess(h).U(i).P(j).P;
                end;
            end;
        end;
    end;
    clear T_i P_j
    
    % Convert from scans to seconds
    %---------------------------------------------------------------------%
    if strcmp(SPM.xBF.UNITS,'scans')
        ons = ons*SPM.xY.RT;
        dur = dur*SPM.xY.RT;
    end;
    
    % Sort trials by onset time
    %---------------------------------------------------------------------%
    [os, is] = sort(ons);
    ind = ind(is);
    ons = ons(is);
    dur = dur(is);
    T   = T(is,:);
    t   = size(T,1);
    names     = cell(1,t);
    onsets    = cell(1,t);
    durations = cell(1,t);
    
    % Get regressors of no interest
    %---------------------------------------------------------------------%
    if ~isempty(SPM.Sess(h).C) % movement parameters
        R = SPM.Sess(h).C.C;
    end;
    
    % Get trial-wise design matrix
    %---------------------------------------------------------------------%
    settings.n  = n;
    settings.TR = SPM.xY.RT;
    for k = 1:t
        names{k}     = strcat('trl-',MF_int2str0(k,2));
        onsets{k}    = ons(k);
        durations{k} = dur(k);
    end;
    if ~isempty(nind)
        [X,  L]  = ITEM_get_des_mat(names,  onsets,  durations,  [],    [], [], settings);
        [Xc, Lc] = ITEM_get_des_mat(names0, onsets0, durations0, pmod0, [],  R, settings);
         X = [X, Xc];
         L = [L, Lc];
         r = size(Xc,2);
    else
        [X, L] = ITEM_get_des_mat(names, onsets, durations, [], [], R, settings);
         r = size(R,2);
    end;
    
    % Add discrete cosine set
    %---------------------------------------------------------------------%
    if strcmp(mode,'DCT') % temporal filter
        k = size(SPM.xX.K(h).X0,2);
        for l = 1:k, L = [L, {strcat('K',num2str(l))}]; end;
        X = [X, SPM.xX.K(h).X0];
        r = r + k;
    end;
    
    % Add constant regressor
    %---------------------------------------------------------------------%
    L = [L, {'const.'}]; % implicit baseline
    X = [X, ones(n,1)];
    r = r + 1;
    
    % Normalize to unit range
    %---------------------------------------------------------------------%
    X = X ./ repmat(max(abs(X)), [n 1]);
    
    % Store session information
    %---------------------------------------------------------------------%
    GLM1.n(h) = n;
    GLM1.p(h) = p;
    GLM1.r(h) = r;
    GLM1.t(h) = t;
    GLM1.pr(h)= p+r;
    GLM1.tr(h)= t+r;
    % first-level (scan-wise) design matrix
    GLM1.Sess(h).n = SPM.Sess(h).row;
    GLM1.Sess(h).X = X;
    GLM1.Sess(h).L = L;
    GLM1.Sess(h).K = SPM.xX.K(h);
    GLM1.Sess(h).V = zeros(GLM1.n(h),GLM1.n(h));
    % second-level (trial-wise) design matrix
    GLM1.Sess(h).t = sum(GLM1.tr(1:h-1))+[1:GLM1.tr(h)];
    GLM1.Sess(h).T = blkdiag(T, eye(r)); % augmented
    GLM1.Sess(h).U = zeros(GLM1.tr(h),GLM1.tr(h));
    % By now, the following relation should hold:
    % X [n x (p+r)] = X_t [n x (t+r)] * T [(t+r) x (p+r)].
    GLM1.Sess(h).ind = ind;
    GLM1.Sess(h).ons = ons;
    GLM1.Sess(h).dur = dur;
    GLM1.Sess(h).is  = is;
    
end;

% Setup SPM for ReML
%-------------------------------------------------------------------------%
SPM_ReML.xY     = SPM.xY;
SPM_ReML.xM     = SPM.xM;
SPM_ReML.xX.X   = blkdiag(GLM1.Sess.X);
if strcmp(mode,'DCT')
    SPM_ReML.xX.K = 1;
end;
if strcmp(mode,'KXY')
    SPM_ReML.xX.K = SPM.xX.K;
end;
SPM_ReML.xX.iB  = cumsum(GLM1.tr);
SPM_ReML.xX.iG  = [];
SPM_ReML.xVi.Vi = SPM.xVi.Vi;

% Estimate non-sphericity
%-------------------------------------------------------------------------%
SPM_ReML.xVi = spm_est_non_sphericity(SPM_ReML);
for h = 1:numel(SPM.Sess)
    GLM1.Sess(h).V = SPM_ReML.xVi.V(SPM.Sess(h).row,SPM.Sess(h).row);
    GLM1.Sess(h).U = (GLM1.Sess(h).X'*inv(GLM1.Sess(h).V)*GLM1.Sess(h).X)^(-1);
end;
clear SPM_ReML


%=========================================================================%
% E S T I M A T I O N   ( 2 ) :   T R I A L - W I S E   E S T I M A T E S %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_est_1st_lvl: estimate (2)');

% Preallocate gamma estimates
%-------------------------------------------------------------------------%
G  = NaN(sum(GLM1.tr),prod(m_dim));
S2 = NaN(numel(SPM.Sess),prod(m_dim));

% Estimate parameters of first-level model
%-------------------------------------------------------------------------%
for h = 1:numel(SPM.Sess)
    
    % Get data and design
    %---------------------------------------------------------------------%
    Yh = Y(GLM1.Sess(h).n,:);
    Xh = GLM1.Sess(h).X;
    Vh = GLM1.Sess(h).V;
    
    % Filter data and design
    %---------------------------------------------------------------------%
    if strcmp(mode,'KXY')
        Xh = spm_filter(GLM1.Sess(h).K, Xh);
        Yh = spm_filter(GLM1.Sess(h).K, Yh);
    end;
    
    % Estimate parameters
    %---------------------------------------------------------------------%
    [G(GLM1.Sess(h).t,m_ind), S2(h,m_ind)] = ITEM_GLM_MLE(Yh, Xh, Vh, sprintf('Estimate trial-wise response amplitudes for session %d',h));
    
end;


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_est_1st_lvl: save');
spm_progress_bar('Init',100,'Save trial-wise parameter estimates...','');

% Initialise image files
%-------------------------------------------------------------------------%
H = MA_init_header(SPM, false);
GLM1.swd = strcat(SPM.swd,'/','ITEM_est_1st_lvl','/');
if ~exist(GLM1.swd,'dir'), mkdir(GLM1.swd); end;
cd(GLM1.swd);

% Save gamma estimates
%-------------------------------------------------------------------------%
d = ceil(sum(GLM1.tr)/100);
for h = 1:numel(SPM.Sess)
    for k = 1:GLM1.tr(h)
        i = GLM1.Sess(h).t(k);
        H.fname   = strcat('gamma_',MF_int2str0(i,4),'.nii');
        H.descrip = sprintf('ITEM_est_1st_lvl: parameter estimate; session %d, trial %d', h, k);
        spm_write_vol(H,reshape(G(i,:),m_dim));
        GLM1.Vgamma(i) = H;
        if mod(i,d) == 0, spm_progress_bar('Set',(i/sum(GLM1.tr))*100); end;
    end;
end;

% Save sigma^2 estimates
%-------------------------------------------------------------------------%
for h = 1:numel(SPM.Sess)
    H.fname   = strcat('sigma_',MF_int2str0(h,4),'.nii');
    H.descrip = sprintf('ITEM_est_1st_lvl: variance estimate; session %d', h);
    spm_write_vol(H,reshape(S2(h,:),m_dim));
    GLM1.Vsigma(h) = H;
    spm_progress_bar('Set',(h/s)*100);
end;

% Save GLM structure
%-------------------------------------------------------------------------%
save(strcat(GLM1.swd,'GLM1.mat'),'GLM1');

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);