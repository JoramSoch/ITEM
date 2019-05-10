function ITEM_review(SPM, step, what, sess)
% _
% Display/Review Function for ITEM Analyses
% FORMAT ITEM_review(SPM, step, what, sess)
%     SPM  - a structure specifying an estimated GLM
%     step - a string specifying the analysis step
%     what - a string specifying what exactly to show
%     sess - an integer specifying the session to show
% 
% FORMAT ITEM_review(SPM, step, what, sess) displays relevant quantities
% from an ITEM analysis step applied to a GLM indicated by SPM.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 04/12/2018, 13:10 (V0.1)
%  Last edit: 10/05/2019, 07:40 (V0.2)


%=========================================================================%
% P R E P A R A T I O N                                                   %
%=========================================================================%

% Get SPM.mat if necessary
%-------------------------------------------------------------------------%
if nargin == 0
    SPM_mat = spm_select(1,'^SPM\.mat$','Select SPM.mat!');
    SPM_dir = fileparts(SPM_mat); load(SPM_mat);
    SPM.swd = SPM_dir;
    ITEM_review(SPM);
    return
end;

% Set analysis step if necessary
%-------------------------------------------------------------------------%
if nargin < 2 || isempty(step)
    list = [];
    inds = [];
    if exist(strcat(SPM.swd,'/','ITEM_est_1st_lvl'),'dir')
        list = [list, {'estimate: first-level'}];
        inds = [inds, 1];
    end;
    if exist(strcat(SPM.swd,'/','ITEM_est_2nd_lvl'),'dir')
        list = [list, {'estimate: second-level'}];
        inds = [inds, 2];
    end;
    if exist(strcat(SPM.swd,'/','ITEM_dec_class'),'dir')
        list = [list, {'decode: classify'}];
        inds = [inds, 3];
    end;
    if exist(strcat(SPM.swd,'/','ITEM_dec_recon'),'dir')
        list = [list, {'decode: reconstruct'}];
        inds = [inds, 4];
    end;
    stind = spm_input('Select analysis step to review', 1, 'm', list, inds);
    steps = {'est-1st-lvl', 'est-2nd-lvl', 'dec-class', 'dec-recon'};
    step  = steps{stind};
end;

% Select decoding analysis if necessary
%-------------------------------------------------------------------------%
if strncmp(step, 'dec', 3)
    if strcmp(step, 'dec-class')
        files = dir(strcat(SPM.swd,'/','ITEM_dec_class','/','ITEM_*.mat'));
        list  = cell(1,numel(files));
        for i = 1:numel(files)
            load(strcat(SPM.swd,'/','ITEM_dec_class','/',files(i).name));
            list{i} = sprintf('classification of %s from %s', ITEM.Class.con, ITEM.Class.reg);
        end;
    end;
    if strcmp(step, 'dec-recon')
        files = dir(strcat(SPM.swd,'/','ITEM_dec_recon','/','ITEM_*.mat'));
        list  = cell(1,numel(files));
        for i = 1:numel(files)
            load(strcat(SPM.swd,'/','ITEM_dec_recon','/',files(i).name));
            list{i} = sprintf('reconstruction of %s from %s', ITEM.Recon.con, ITEM.Recon.reg);
        end;
    end;
    fidx = spm_input('Select decoding analysis to inspect', '+1', 'm', list, [1:numel(files)]);
    file = files(fidx).name;
end;

% Select what exactly to show if necessary
%-------------------------------------------------------------------------%
if nargin < 3 || isempty(what)
    if strcmp(step, 'est-1st-lvl')
        list  = {'design matrices', 'induced covariance'};
        whats = {'des-mat', 'ind-cov'};
    end;
    if strcmp(step, 'est-2nd-lvl')
        list  = {'covariance matrices', 'covariance components'};
        whats = {'cov-mat', 'cov-comp'};
    end;
    if strcmp(step, 'dec-class')
        list  = {'inverted model', 'model inversion', 'decoding accuracies'};
        whats = {'inv-mod', 'mod-inv', 'dec-acc'};
    end;
    if strcmp(step, 'dec-recon')
        list  = {'inverted model', 'model inversion', 'correlation coefficients'};
        whats = {'inv-mod', 'mod-inv', 'corr-coeff'};
    end;
    whind = spm_input('Select what exactly to show', '+1', 'm', list, [1:numel(list)]);
    what  = whats{whind};
end;

% Load ITEM.mat and get number of sessions
%-------------------------------------------------------------------------%
switch step
    case 'est-1st-lvl'
        ITEM_mat = strcat(SPM.swd,'/','ITEM_est_1st_lvl','/','GLM1.mat');
    case 'est-2nd-lvl'
        ITEM_mat = strcat(SPM.swd,'/','ITEM_est_2nd_lvl','/','GLM2.mat');
    case 'dec-class'
        ITEM_mat = strcat(SPM.swd,'/','ITEM_dec_class','/',file);
    case 'dec-recon'
        ITEM_mat = strcat(SPM.swd,'/','ITEM_dec_recon','/',file);
end;
load(ITEM_mat);
switch step
    case 'est-1st-lvl', s = numel(GLM1.Sess);
    case 'est-2nd-lvl', s = numel(GLM2.Sess);
    case 'dec-class',   s = numel(ITEM.Sess);
    case 'dec-recon',   s = numel(ITEM.Sess);
end;

% Select the session to dispay if necessary
%-------------------------------------------------------------------------%
if ~strcmp(what,'dec-acc') & ~strcmp(what,'corr-coeff')
    if nargin < 4 || isempty(sess)
        list = cellstr([repmat('Session ',[s 1]),num2str([1:s]')])';
        sess = spm_input('Select the session to display:', '+1', 'm', list, [1:s]);
    end;
    h = sess;
end;


%=========================================================================%
% D I S P L A Y                                                           %
%=========================================================================%

% Case: estimate first-level model
%-------------------------------------------------------------------------%
if strcmp(step, 'est-1st-lvl')
    
    % Display design matrices
    %---------------------------------------------------------------------%
    if strcmp(what,'des-mat')
        % open figure
        figure('Name', sprintf('ITEM_est_1st_lvl: design matrices (Session %d)', h), 'Color', [1 1 1], 'Position', [50 50 1600 900]);
        colormap gray;
        % extract matrices
        X  = GLM1.Sess(h).X*GLM1.Sess(h).T;
        Xt = GLM1.Sess(h).X;
        T  = GLM1.Sess(h).T;
        % standard design matrix
        subplot(1,3,1);
        imagesc(X);
        caxis([-max(max(abs(X))), +max(max(abs(X)))]);
        axis off;
        title(sprintf('X = X_t T [%d x %d]', GLM1.n(h), GLM1.pr(h)), 'FontSize', 12);
        % first-level (scan-wise) design matrix
        subplot(1,3,2);
        imagesc(Xt);
        caxis([-1, +1]);
        axis off;
        title(sprintf('X_t [%d x %d]', GLM1.n(h), GLM1.tr(h)), 'FontSize', 12);
        % second-level (trial-wise) design matrix
        subplot(1,3,3);
        imagesc(T);
        caxis([-1, +1]);
        axis off;
        title(sprintf('T_{ } [%d x %d]', GLM1.tr(h), GLM1.pr(h)), 'FontSize', 12);
        % delete matrices
        clear X Xt T
    end;
    
    % Display induced covariance
    %---------------------------------------------------------------------%
    if strcmp(what,'ind-cov')
        % open figure
        figure('Name', sprintf('ITEM_est_1st_lvl: induced covariance (Session %d)', h), 'Color', [1 1 1], 'Position', [50 50 1600 900]);
        colormap jet;
        % extract matrices
        Xt = GLM1.Sess(h).X;
        Vi = inv(GLM1.Sess(h).V);
        Ui = inv(GLM1.Sess(h).U);
        um = max(max(abs(Ui(1:GLM1.t(h),1:GLM1.t(h)))));
        % trial-wise design matrix
        subplot(1,5,1);
        imagesc(Xt);
        caxis([-1, +1]);
        axis off;
        title(sprintf('X_t [%d x %d]', GLM1.n(h), GLM1.tr(h)), 'FontSize', 12);
        % temporal covariance matrix
        subplot(1,5,[2,3]);
        imagesc(Vi);
        caxis([0, +1]);
        axis square off;
        title(sprintf('V^{-1}_{ } [%d x %d]', GLM1.n(h), GLM1.n(h)), 'FontSize', 12);
        % variance-covariance matrix
        subplot(1,5,[4,5]);
        imagesc(Ui);
        caxis([-um, +um]);
        axis square off;
        title(sprintf('U^{-1} = X_t^T V^{-1} X_t [%d x %d]', GLM1.tr(h), GLM1.tr(h)), 'FontSize', 12);
        % delete matrices
        clear Xt Vi Ui um
    end;
    
end;

% Case: estimate second-level model
%-------------------------------------------------------------------------%
if strcmp(step, 'est-2nd-lvl')
    
    % Display covariance matrices
    %---------------------------------------------------------------------%
    if strcmp(what,'cov-mat')
        % open figure
        figure('Name', sprintf('ITEM_est_2nd_lvl: covariance matrices (Session %d)', h), 'Color', [1 1 1], 'Position', [50 50 1600 900]);
        colormap jet;
        % extract matrices
        load(strcat(SPM.swd,'/','ITEM_est_1st_lvl','/','GLM1.mat'));
        t  = GLM1.t; clear GLM1
        V  = GLM2.Sess(h).V;
        Q1 = GLM2.Sess(h).Q{1};
        Q2 = GLM2.Sess(h).Q{2};
        vm = max(max(abs(V(1:t(h),1:t(h)))));
        um = max(max(abs(Q2(1:t(h),1:t(h)))));
        % covariance matrix
        subplot(1,3,1);
        imagesc(V);
        caxis([-vm, +vm]);
        axis square off;
        title(sprintf('V_{ } [%d x %d]', size(V)), 'FontSize', 12);
        % natural covariance
        subplot(1,3,2);
        imagesc(Q1);
        caxis([0, 1]);
        axis square off;
        title(sprintf('= %3.3f x I_t', GLM2.Sess(h).s2(1)), 'FontSize', 12);
        % induced covariance
        subplot(1,3,3);
        imagesc(Q2);
        caxis([-um, +um]);
        axis square off;
        title(sprintf('+ %3.3f x U_{ }', GLM2.Sess(h).s2(2)), 'FontSize', 12);
        % delete matrices
        clear t V Q1 Q2 Vi Ui vm um
    end;
    
    % Display covariance components
    %---------------------------------------------------------------------%
    if strcmp(what,'cov-comp')
        % open figure
        figure('Name', sprintf('ITEM_est_2nd_lvl: covariance components (Session %d)', h), 'Color', [1 1 1], 'Position', [50 50 1600 900]);
        colormap jet;
        % extract matrices
        load(strcat(SPM.swd,'/','ITEM_est_1st_lvl','/','GLM1.mat'));
        t  = GLM1.t; clear GLM1
        Q1 = GLM2.Sess(h).Q{1};
        Q2 = GLM2.Sess(h).Q{2};
        um = max(max(abs(Q2(1:t(h),1:t(h)))));
        % natural covariance
        subplot(1,2,1);
        imagesc(Q1);
        caxis([0, 1]);
        axis square off;
        title(sprintf('Q_1 = I_t [%d x %d]', size(Q1)), 'FontSize', 12);
        % induced covariance
        subplot(1,2,2);
        imagesc(Q2);
        caxis([-um, +um]);
        axis square off;
        title(sprintf('Q_2 = U_{ } [%d x %d]', size(Q2)), 'FontSize', 12);
        % delete matrices
        clear t Q1 Q2 um
    end;
    
end;

% Cases: decode by classification or reconstruction
%-------------------------------------------------------------------------%
if strcmp(step, 'dec-class') || strcmp(step, 'dec-recon')
    
    % Display inverted model
    %---------------------------------------------------------------------%
    if strcmp(what,'inv-mod')
        % open figure
        figure('Name', sprintf('ITEM_dec_%s: inverted model (Session %d)', step(5:end), h), 'Color', [1 1 1], 'Position', [50 50 1600 900]);
        colormap jet;
        % extract matrices
        load(strcat(SPM.swd,'/','ITEM_est_1st_lvl','/','GLM1.mat'));
        t  = GLM1.t; clear GLM1
        Y  = ITEM.Sess(h).Y;
        X  = ITEM.Sess(h).X;
        V  = ITEM.Sess(h).V;
        xm = max(max(abs(X(1:end-1,:))));
        vm = max(max(abs(V(1:t(h),1:t(h)))));
        % data matrix
        subplot(1,4,1);
        imagesc(Y);
        caxis([-1, +1]);
        axis off;
        title(sprintf('"data" T [%d x %d]', size(Y)), 'FontSize', 12);
        % design matrix
        subplot(1,4,2);
        imagesc(X);
        caxis([-xm, +xm]);
        axis off;
        title(['"design" \Gamma', sprintf(' [%d x %d]', size(X))], 'FontSize', 12);
        % covariance matrix
        subplot(1,4,[3,4]);
        imagesc(V);
        caxis([-vm, +vm]);
        axis square off;
        title(['"covariance" U', sprintf(' [%d x %d]', size(V))], 'FontSize', 12);
        % delete matrices
        clear t Y X V xm vm
    end;
    
    % Display model inversion
    %---------------------------------------------------------------------%
    if strcmp(what,'mod-inv')
        % open figure
        figure('Name', sprintf('ITEM_dec_%s: model inversion (Session %d)', step(5:end), h), 'Color', [1 1 1], 'Position', [50 50 1600 900]);
        colormap gray;
        % extract matrices
        Yp = ITEM.Sess(h).Yp;
        if strcmp(step, 'dec-class'), Yc = ITEM.Sess(h).Yc; end;
        if strcmp(step, 'dec-recon'), Yr = ITEM.Sess(h).Yr; end;
        Yt = ITEM.Sess(h).Yt;
        mp = max(max(abs(Yp(1:end-1,1:end-1))));
        % predicted design matrix
        subplot(1,4,[1,2]);
        imagesc(Yp);
        caxis([-mp, +mp]);
        set(gca,'Box','On');
        set(gca,'XTick',[],'YTick',[]);
        title('predicted design matrix', 'FontSize', 12);
        % classified conditions / reconstructed variables
        subplot(1,4,3);
        if strcmp(step, 'dec-class'), imagesc(Yc); end;
        if strcmp(step, 'dec-recon'), imagesc(Yr); end;
        caxis([-1, +1]);
        set(gca,'Box','On');
        set(gca,'XTick',[],'YTick',[]);
        if strcmp(step, 'dec-class'), title('classfied conditions', 'FontSize', 12);    end;
        if strcmp(step, 'dec-recon'), title('reconstructed variables', 'FontSize', 12); end;
        % true conditions / true variables
        subplot(1,4,4);
        imagesc(Yt);
        caxis([-1, +1]);
        set(gca,'Box','On');
        set(gca,'XTick',[],'YTick',[]);
        if strcmp(step, 'dec-class'), title('true conditions', 'FontSize', 12); end;
        if strcmp(step, 'dec-recon'), title('true variables', 'FontSize', 12);  end;
        % delete matrices
        clear Yp Yc Yr Yt mp
    end;
    
end;
    
% Case: decode by classification
%-------------------------------------------------------------------------%
if strcmp(step, 'dec-class')
    
    % Display decoding accuracies
    %---------------------------------------------------------------------%
    if strcmp(what,'dec-acc')
        figure('Name', 'ITEM_dec_class: decoding accuracies (all sessions)', 'Color', [1 1 1], 'Position', [50 50 1600 900]);
        q  = size(ITEM.Class.c,1);
        DA = horzcat(ITEM.Sess.DA);
        hold on;
        bar(0, mean(DA), 'r');
        bar([1:s], DA, 'b');
        plot([-100, 100], [1/q, 1/q], '-k', 'LineWidth', 2);
        axis([-0.5, max([s+0.5, 10.5]), 0, 1]);
        grid on;
        set(gca,'Box', 'On');
        set(gca,'XTick', [0, 1:s], 'XTickLabel', [{'avg'}, cellstr(num2str([1:s]'))']);
        legend('cvDA', 'oosDA', 'chance', 'Location', 'NorthEast');
        xlabel('recording session', 'FontSize', 12);
        ylabel('decoding accuracy', 'FontSize', 12);
        title(sprintf('Classification of %s from %s', ITEM.Class.con, ITEM.Class.reg), 'FontSize', 16);
    end;
    
end;

% Case: decode by reconstruction
%-------------------------------------------------------------------------%
if strcmp(step, 'dec-recon')
    
    % Display correlation coefficients
    %---------------------------------------------------------------------%
    if strcmp(what,'corr-coeff')
        figure('Name', 'ITEM_dec_recon: correlation coefficients (all sessions)', 'Color', [1 1 1], 'Position', [50 50 1600 900]);
        c  = ITEM.Recon.c;
        CC = vertcat(ITEM.Sess.CC);
        hold on;
        bar([-1,0], [mean(CC); mean(CC)], 'grouped');
        bar([1:s], CC, 'grouped');
        axis([-0.5, max([s+0.5, 10.5]), -1, +1]);
        grid on;
        set(gca,'Box', 'On');
        set(gca,'XTick', [0, 1:s], 'XTickLabel', [{'avg'}, cellstr(num2str([1:s]'))']);
        legend(cellstr(num2str(find(c)'))', 'Location', 'NorthEast');
        xlabel('recording session', 'FontSize', 12);
        ylabel('correlation coefficient', 'FontSize', 12);
        title(sprintf('Reconstruction of %s from %s', ITEM.Recon.con, ITEM.Recon.reg), 'FontSize', 16);
    end;
    
end;