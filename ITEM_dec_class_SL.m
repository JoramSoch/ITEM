function ITEM_dec_class_SL(SPM, rad, c, con)
% _
% Decoding from Trials for Classification (searchlight-based)
% FORMAT ITEM_dec_class_SL(SPM, rad, c, con)
%     SPM - a structure specifying an estimated GLM
%     rad - a scalar specifying the searchlight radius in mm (e.g. 6)
%     c   - a 1 x p contrast vector with +1s & -1s for binary decoding or
%           a q x p contrast matrix with +1s only for n-ary decoding
%     con - a string without spaces describing the contrast (e.g. 'WM')
% 
% FORMAT ITEM_dec_class_SL(SPM, rad, c, con) performs searchlight decoding
% by an inverse transformed encoding model for classification between
% classes indicated by c using searchlights of size rad.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 10/05/2019, 11:15 (V0.2)
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
    ITEM_dec_class_SL(SPM);
    return
end;

% Set region of interest if necessary
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(rad)
    rad = 2*abs(SPM.xY.VY(1).mat(1,1));
end;

% Set decoding contrast if necessary
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(c)
    c = zeros(2,0);
    c(1,SPM.Sess(1).U(1).P.i(1)) = 1;
    c(2,SPM.Sess(1).U(2).P.i(1)) = 1;
end;

% Set contrast name if necessary
%-------------------------------------------------------------------------%
if nargin < 4 || isempty(con)
    if size(c,1) == 1, d = [(1/2)*(c+1); (-1/2)*(c-1)]; else, d = c; end;
    con = '';
    for k = 1:size(d,1)
        dk = find(d(k,:));
        for l = 1:numel(dk)
            con = strcat(con,int2str(dk(l)));
            if l < numel(dk), con = strcat(con,','); end;
        end;
        if k < size(d,1), con = strcat(con,'-vs-'); end;
    end;
    clear d dk
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
Finter = spm('FigName','ITEM_dec_class_SL: load');

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
Finter = spm('FigName','ITEM_dec_class_SL: estimate (1)');
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
Finter = spm('FigName','ITEM_dec_class_SL: estimate (2)');

% Augment decoding contrast if necessary
%-------------------------------------------------------------------------%
if size(c,1) == 1, c = [1*(c==1); 1*(c==-1)]; end;
c = [c, zeros(size(c,1), GLM1.p(1)-size(c,2))];
q = size(c,1);

% Cycle through recording sessions
%-------------------------------------------------------------------------%
for h = 1:s
    
    % "data" - the T matrix
    %---------------------------------------------------------------------%
    Th = GLM1.Sess(h).T(1:GLM1.t(h),1:GLM1.p(h));
    Yh = Th * c';
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
Finter = spm('FigName','ITEM_dec_class_SL: estimate (3)');

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
    Y_out = ITEM.Sess(g).Y;
    X_out = ITEM.Sess(g).X;
    V_out = ITEM.Sess(g).V;
    
    % Perform searchlight-based ITEM analysis
    %---------------------------------------------------------------------%
    ITEM.Sess(g).Yp = ITEM_ITEM_SL(Y_in, X_in, V_in, X_out, V_out, SLs, sprintf('Searchlight-based classification for session %d',g));
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
Finter = spm('FigName','ITEM_dec_class_SL: estimate (3)');

% Preallocate oos & cv decoding accuracies
%-------------------------------------------------------------------------%
oosDA = NaN(1,numel(M),s);
cvDA  = NaN(1,numel(M));

% Calculate out-of-sample decoding accuracies
%-------------------------------------------------------------------------%
for g = 1:s
    spm_progress_bar('Init', 100, sprintf('Calculate decoding accuracy for session %d',g), '');
    Y_true  = ITEM.Sess(g).Y;
    Y_class = ITEM.Sess(g).Yp;
    i_eff   = find(sum(Y_true,2)>0)';
    t_eff   = numel(i_eff);
    for j = 1:v
        for i = i_eff
            k = find(Y_class(i,1:q,j)==max(Y_class(i,1:q,j)));
            Y_class(i,k,j) = 1;
        end;
        oosDA(1,m_ind(j),g) = (1/t_eff) * sum(diag(Y_true'*Y_class(:,:,j)));
        if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
    end;
end;
avgDA = mean(oosDA,3);
clear Y_true Y_class i_eff t_eff

% Calculate cross-validated decoding accuracy
%-------------------------------------------------------------------------%
spm_progress_bar('Init', 100, 'Calculate decoding accuracy across all sessions', '');
Y_true  = vertcat(ITEM.Sess(1:s).Y);
Y_class = vertcat(ITEM.Sess(1:s).Yp);
i_eff   = find(sum(Y_true,2)>0)';
t_eff   = numel(i_eff);
for j = 1:v
    for i = i_eff
        k = find(Y_class(i,1:q,j)==max(Y_class(i,1:q,j)));
        Y_class(i,k,j) = 1;
    end;
    cvDA(m_ind(j)) = (1/t_eff) * sum(diag(Y_true'*Y_class(:,:,j)));
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;
clear Y_true Y_class i_eff t_eff

% Clear progress bar
%-------------------------------------------------------------------------%
spm_progress_bar('Clear');


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_dec_class_SL: save');

% Initialize image files
%-------------------------------------------------------------------------%
H = MA_init_header(SPM, false);
ITEM.swd = strcat(SPM.swd,'/','ITEM_dec_class','/','ITEM_',con,'_SL-',num2str(rad),'mm','/');
if ~exist(ITEM.swd,'dir'), mkdir(ITEM.swd); end;
cd(ITEM.swd);

% Save out-of-sample accuracies
%-------------------------------------------------------------------------%
d = floor(log10(s))+1;
for h = 1:s
    H.fname   = strcat('oosDA_S',MF_int2str0(h,d),'.nii');
    H.descrip = sprintf('ITEM_dec_class_SL: out-of-sample decoding accuracy; session %d', h);
    spm_write_vol(H,reshape(oosDA(1,:,h),m_dim));
    ITEM.VoosDA(h,1) = H;
end;

% Save average accuracy
%-------------------------------------------------------------------------%
H.fname   = strcat('avgDA.nii');
H.descrip = sprintf('ITEM_dec_class_SL: average decoding accuracy; %d sessions', s);
spm_write_vol(H,reshape(avgDA,m_dim));
ITEM.VavgDA = H;

% Save cross-validated accuracy
%-------------------------------------------------------------------------%
H.fname   = strcat('cvDA.nii');
H.descrip = sprintf('ITEM_dec_class_SL: cross-validated decoding accuracy; %d sessions', s);
spm_write_vol(H,reshape(cvDA,m_dim));
ITEM.VcvDA = H;

% Save voxels per searchlight
%-------------------------------------------------------------------------%
H.fname   = strcat('VpSL.nii');
H.descrip = sprintf('ITEM_dec_class_SL: voxels per searchlight; %s mm radius', num2str(rad));
spm_write_vol(H,reshape(VpSL,m_dim));
ITEM.VVpSL = H;

% Complete ITEM structure
%-------------------------------------------------------------------------%
ITEM.rad       = rad;
ITEM.SLs       = SLs;
ITEM.Class.c   = c;
ITEM.Class.con = con;

% Save ITEM structure
%-------------------------------------------------------------------------%
save(strcat(ITEM.swd,'ITEM.mat'),'ITEM','-v7.3');

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);