function ITEM_dec_recon_SL(SPM, rad, c, con)
% _
% Decoding from Trials for Reconstruction (searchlight-based)
% FORMAT ITEM_dec_recon_SL(SPM, rad, c, con)
%     SPM - a structure specifying an estimated GLM
%     rad - a scalar specifying the searchlight radius in mm (e.g. 6)
%     c   - a 1 x p vector with +1s indicating variables to reconstruct or
%           a q x p matrix with one +1 in each row indexing each variable
%     con - a string without spaces describing the contrast (e.g. 'PE')
% 
% FORMAT ITEM_dec_recon_SL(SPM, rad, c, con) performs searchlight decoding
% by an inverse transformed encoding model for reconstruction of selected
% variables indicated by c using searchlights of size rad.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 10/05/2019, 14:50 (V0.2)
%  Last edit: 30/11/2021, 15:53 (V0.3)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get SPM.mat if necessary
%-------------------------------------------------------------------------%
if nargin == 0
    SPM_mat = spm_select(1,'^SPM\.mat$','Select SPM.mat!');
    SPM_dir = fileparts(SPM_mat); load(SPM_mat);
    SPM.swd = SPM_dir;
    ITEM_dec_recon_SL(SPM);
    return
end;

% Set searchlight radius if necessary
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(rad)
    rad = 2*abs(SPM.xY.VY(1).mat(1,1));
end;

% Set decoding contrast if necessary
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(c)
    c = zeros(1,0);
    k = 0;
    while isempty(c) && k < numel(SPM.Sess(1).U(1))
    	k = k + 1;
        if ~strcmp(SPM.Sess(1).U(k).P(1).name,'none')
            c(1,1) = SPM.Sess(1).U(k).P(1).i(2);
        end;
    end;
    if isempty(c)
        c(1,1) = 1;
    end;
end;

% Set contrast name if necessary
%-------------------------------------------------------------------------%
if nargin < 4 || isempty(con)
    if size(c,1) > 1, d = sum(c,1); else, d = c; end;
    con = '';
    ind = find(d);
    for l = 1:numel(ind)
        con = strcat(con,int2str(ind(l)));
        if l < numel(ind), con = strcat(con,','); end;
    end;
    clear d ind
end;

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
Finter = spm('FigName','ITEM_dec_recon_SL: load');

% Load mask image
%-------------------------------------------------------------------------%
[M, m_dim, m_ind] = MA_load_mask(SPM);
[m_img, m_xyz]    = spm_read_vols(SPM.VM);
clear m_img

% Load gamma estimates
%-------------------------------------------------------------------------%
G = ITEM_load_gammas(SPM, m_ind);
v = numel(m_ind);
d = floor(v/100);


%=========================================================================%
% E S T I M A T I O N   ( 1 ) :   P R E - C A L C U L A T I O N S         %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_dec_recon_SL: estimate (1)');
spm_progress_bar('Init', 100, 'Determine searchlight voxels...', '');

% Get voxels per searchlight
%-------------------------------------------------------------------------%
SLs  = cell(v,1);
VpSL = NaN(size(M));
XYZ  = m_xyz(:,m_ind);
for j = 1:v
    xyz_cent = m_xyz(:,m_ind(j));
    vox_ind  = find(sqrt(sum((XYZ - repmat(xyz_cent,[1 v])).^2)) <= rad);
    SLs{j}   = vox_ind;
    VpSL(m_ind(j)) = numel(vox_ind);
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;
clear xyz_cent vox_ind

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');


%=========================================================================%
% E S T I M A T I O N   ( 2 ) :   I N V E R S E   M O D E L               %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_dec_recon_SL: estimate (2)');

% Augment decoding contrast if necessary
%-------------------------------------------------------------------------%
if size(c,1) > 1, c = sum(c,1); end;
c = [c, zeros(1, GLM1.p(1)-numel(c))];
q = numel(find(c));

% Cycle through recording sessions
%-------------------------------------------------------------------------%
for h = 1:s
    
    % "data" - the T matrix
    %---------------------------------------------------------------------%
    Th = GLM1.Sess(h).T(1:GLM1.t(h),1:GLM1.p(h));
    Yh = Th(:,c==1);
    ITEM.Sess(h).Y = Yh;
    clear Yh Th
    
    % "design" - gamma estimates
    %---------------------------------------------------------------------%
    Xh = G(GLM1.Sess(h).t(1:GLM1.t(h)),:);
    ITEM.Sess(h).X = Xh;
    clear Xh
    
    % "covariance" - the U matrix
    %---------------------------------------------------------------------%
    % Yh  = G(GLM1.Sess(h).t(1:GLM1.t(h)),:);
    % Xh  = GLM1.Sess(h).T(1:GLM1.t(h),1:GLM1.p(h));
    Qh{1} = eye(GLM1.t(h));
    Qh{2} = GLM1.Sess(h).U(1:GLM1.t(h),1:GLM1.t(h));
    
    % Restricted maximum likelihood
    %---------------------------------------------------------------------%
    % [Vh, s2] = ITEM_GLM_ReML(Yh, Xh, Qh{1}, Qh{2}, sprintf('ITEM_dec_class: ReML estimation for session %d',h));
    ITEM.Sess(h).V  = Qh{2}; % Vh;
    ITEM.Sess(h).s2 = [0 1]; % s2;
    % clear Yh Xh Qh s2
    clear Qh
    
end;


%=========================================================================%
% E S T I M A T I O N   ( 3 ) :   C R O S S - V A L I D A T I O N         %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_dec_recon_SL: estimate (3)');

% Cycle through cross-validation folds
%-------------------------------------------------------------------------%
for g = 1:s
    
    % List sessions for this fold
    %---------------------------------------------------------------------%
    fold = 1:s;
    fold = fold(fold~=g);
    
    % Establish (in-sample) training data set
    %---------------------------------------------------------------------%
    Y_in = vertcat(ITEM.Sess(fold).Y);
    X_in = vertcat(ITEM.Sess(fold).X);
    V_in = blkdiag(ITEM.Sess(fold).V);
    
    % Establish (out-of-sample) test data set
    %---------------------------------------------------------------------%
    X_out = ITEM.Sess(g).X;
    V_out = ITEM.Sess(g).V;
    
    % Perform searchlight-based ITEM analysis
    %---------------------------------------------------------------------%
    ITEM.Sess(g).Yp = ITEM_ITEM_SL(Y_in, X_in, V_in, X_out, V_out, SLs, sprintf('Searchlight-based reconstruction of session %d',g));
    clear Y_in X_in V_in X_out V_out
    
end;

% Remove data matrix from ITEM structure
%-------------------------------------------------------------------------%
ITEM.Sess = rmfield(ITEM.Sess,'X');
ITEM.Sess = rmfield(ITEM.Sess,'V');


%=========================================================================%
% E S T I M A T I O N   ( 4 ) :   D E C O D I N G   A C C U R A C Y       %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_dec_recon_SL: estimate (4)');

% Preallocate oos & cv correlation coefficients
%-------------------------------------------------------------------------%
oosCC = NaN(q,numel(M),s);
cvCC  = NaN(q,numel(M));

% Calculate out-of-sample correlation coefficients
%-------------------------------------------------------------------------%
for g = 1:s
    Y_true  = ITEM.Sess(g).Y;
    Y_recon = ITEM.Sess(g).Yp;
    for k = 1:q
        spm_progress_bar('Init', 100, sprintf('Calculate correlation coefficient for session %d, variable %d',g,k), '');
        i_eff = find(Y_true(:,k)~=0)';
        for j = 1:v
            oosCC(k,m_ind(j),g) = corr(Y_true(i_eff,k),Y_recon(i_eff,k,j));
            if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
        end;
    end;
end;
avgCC = mean(oosCC,3);
clear Y_true Y_recon i_eff

% Calculate cross-validated correlation coefficients
%-------------------------------------------------------------------------%
Y_true  = vertcat(ITEM.Sess(1:s).Y);
Y_recon = vertcat(ITEM.Sess(1:s).Yp);
for k = 1:q
    spm_progress_bar('Init', 100, sprintf('Calculate correlation coefficient across all sessions, variable %d',k), '');
    i_eff = find(Y_true(:,k)~=0)';
    for j = 1:v
        cvCC(k,m_ind(j)) = corr(Y_true(i_eff,k),Y_recon(i_eff,k,j));
        if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
    end;
end;
clear Y_true Y_recon i_eff

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_dec_recon_SL: save');

% Initialize image files
%-------------------------------------------------------------------------%
H = MA_init_header(SPM, false);
ITEM.swd = strcat(SPM.swd,'/','ITEM_dec_recon','/','ITEM_',con,'_SL-',num2str(rad),'mm','/');
if ~exist(ITEM.swd,'dir'), mkdir(ITEM.swd); end;
cd(ITEM.swd);

% Save out-of-sample correlations
%-------------------------------------------------------------------------%
i = find(c);
d = floor(log10(s))+1;
for h = 1:s
    for k = 1:q
        H.fname   = strcat('oosCC_',MF_int2str0(i(k),4),'_S',MF_int2str0(h,d),'.nii');
        H.descrip = sprintf('ITEM_dec_recon_SL: out-of-sample correlation coefficient; session %d, regressor %d', h, i(k));
        spm_write_vol(H,reshape(oosCC(k,:,h),m_dim));
        ITEM.VoosCC(h,k) = H;
    end;
end;

% Save averaged correlations
%-------------------------------------------------------------------------%
for k = 1:q
    H.fname   = strcat('avgCC_',MF_int2str0(i(k),4),'.nii');
    H.descrip = sprintf('ITEM_dec_recon_SL: averaged correlation coefficient; %d sessions, regressor %d', s, i(k));
    spm_write_vol(H,reshape(avgCC(k,:),m_dim));
    ITEM.VavgCC(1,k) = H;
end;

% Save cross-validated correlations
%-------------------------------------------------------------------------%
for k = 1:q
    H.fname   = strcat('cvCC_',MF_int2str0(i(k),4),'.nii');
    H.descrip = sprintf('ITEM_dec_recon_SL: cross-validated correlation coefficient; %d sessions, regressor %d', s, i(k));
    spm_write_vol(H,reshape(cvCC(k,:),m_dim));
    ITEM.VcvCC(1,k) = H;
end;

% Save voxels per searchlight
%-------------------------------------------------------------------------%
H.fname   = strcat('VpSL.nii');
H.descrip = sprintf('ITEM_dec_recon_SL: voxels per searchlight; %s mm radius', num2str(rad));
spm_write_vol(H,reshape(VpSL,m_dim));
ITEM.VVpSL = H;

% Complete ITEM structure
%-------------------------------------------------------------------------%
ITEM.rad       = rad;
ITEM.SLs       = SLs;
ITEM.Recon.c   = c;
ITEM.Recon.con = con;

% Save ITEM structure
%-------------------------------------------------------------------------%
save(strcat(ITEM.swd,'ITEM.mat'),'ITEM','-v7.3');

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);