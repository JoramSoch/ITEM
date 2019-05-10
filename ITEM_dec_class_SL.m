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
%  Last edit: 10/05/2019, 14:20 (V0.2)


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
    if size(c,1) == 1, d = [(1/2)*(c+1); (-1/2)*(c-1)]; else d = c; end;
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
[M m_dim m_ind] = MA_load_mask(SPM);
[m_img m_xyz]   = spm_read_vols(SPM.VM);
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
Finter = spm('FigName','ITEM_dec_class_SL: estimate(1)');
spm_progress_bar('Init', 100, 'Determine searchlight voxels...', '');

% Get voxels per searchlight
%-------------------------------------------------------------------------%
v2v   = false(v,v);
vXv   = false(v,v);
% NOTE: The logical matrix "v2v" indicates whether a voxel (column) belongs
%       to the searchlight around a voxel (row). The logical matrix "vXv"
%       indicates whether there is a searchlight that two voxels (column/
%       row) both belong to.
VpSL  = NaN(size(M));
xyz   = m_xyz(:,m_ind);
for j = 1:v
    xyz_cent = m_xyz(:,m_ind(j));
    v2v_ind  = find(sqrt(sum((xyz - repmat(xyz_cent,[1 v])).^2)) <= rad);
    v2v(j,v2v_ind) = true;                  % voxels belonging to one SL
    VpSL(m_ind(j)) = numel(v2v_ind);        % number of voxels in this SL
    vXv(v2v_ind,v2v_ind) = true;            % all voxel-pairs in this SL
    if mod(j,d) == 0, spm_progress_bar('Set',(j/v)*100); end;
end;

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
c = [c, zeros(size(c,1), GLM1.pr(1)-size(c,2))];
q = size(c,1);

% Cycle through recording sessions
%-------------------------------------------------------------------------%
for h = 1:s
    
    % "data" - the T matrix
    %---------------------------------------------------------------------%
    Th = GLM1.Sess(h).T;
    Yh = Th * c';
    ITEM.Sess(h).Y = Yh;
    clear Yh Th
    
    % "design" - gamma estimates
    %---------------------------------------------------------------------%
    Xh = G(GLM1.Sess(h).t,:);
    ITEM.Sess(h).X = Xh;
    clear Xh
    
    % "covariance" - the U matrix
    %---------------------------------------------------------------------%
    Yh    = G(GLM1.Sess(h).t,:);
    Xh    = GLM1.Sess(h).T;
    Qh{1} = eye(GLM1.tr(h));
    Qh{2} = GLM1.Sess(h).U;
    
    % Restricted maximum likelihood
    %---------------------------------------------------------------------%
    % [Vh, s2] = ITEM_GLM_ReML(Yh, Xh, Qh{1}, Qh{2}, sprintf('ITEM_dec_class_SL: ReML estimation for session %d',h));
    ITEM.Sess(h).V  = Qh{2}; % Vh;
    ITEM.Sess(h).s2 = [0 1]; % s2;
    clear Yh Xh Qh Vh s2
    
end;


%=========================================================================%
% E S T I M A T I O N   ( 3 ) :   C R O S S - V A L I D A T I O N         %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_dec_class_SL: estimate (3)');

% Preallocate (out-of-sample) decoding accuracy
%-------------------------------------------------------------------------%
oosDA = NaN(s,numel(M));

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
    oosDA(g,m_ind) = ITEM_ITEM_SL(Y_in, X_in, V_in, Y_out, X_out, V_out, v2v, vXv, 'class', sprintf('Searchlight-based classification for session %d',g));
    
end;

% Calculate (cross-validated) decoding accuracy
%-------------------------------------------------------------------------%
cvDA = mean(oosDA,1);


%=========================================================================%
% S A V E   R E S U L T S                                                 %
%=========================================================================%

% Init progress bar
%-------------------------------------------------------------------------%
Finter = spm('FigName','ITEM_dec_class_SL: save');

% Initialize image files
%-------------------------------------------------------------------------%
H = MA_init_header(SPM, false);
ITEM.swd = strcat(SPM.swd,'/','ITEM_dec_class','/','ITEM_',con,'_',num2str(rad),'mm','/');
if ~exist(ITEM.swd,'dir'), mkdir(ITEM.swd); end;
cd(ITEM.swd);

% Save out-of-sample accuracies
%-------------------------------------------------------------------------%
d = floor(log10(s))+1;
for h = 1:s
    H.fname   = strcat('oosDA_S',MF_int2str0(h,d),'.nii');
    H.descrip = sprintf('ITEM_dec_class_SL: out-of-sample decoding accuracy; session %d', h);
    spm_write_vol(H,reshape(oosDA(h,:),m_dim));
    ITEM.VoosDA(h) = H;
end;

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
ITEM.Class.c   = c;
ITEM.Class.con = con;

% Save ITEM structure
%-------------------------------------------------------------------------%
save(strcat(ITEM.swd,'ITEM.mat'),'ITEM');

% Return to origin
%-------------------------------------------------------------------------%
cd(orig_dir);