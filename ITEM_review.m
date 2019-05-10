function ITEM_review(SPM, step, what)
% _
% Display/Review Function for ITEM Analyses
% FORMAT ITEM_review(SPM, step)
%     SPM  - a structure specifying an estimated GLM
%     step - a string specifying the analysis step
%     what - a string specifying what exactly to show
% 
% FORMAT ITEM_review(SPM, step) displays relevant quantities from an
% ITEM analysis step applied to a GLM indicated by SPM.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 04/12/2018, 13:10 (V0.1)
%  Last edit: 08/05/2019, 17:05 (V0.1)


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
    steps = {'estimate', 'decode'};
    step  = steps{spm_input('Analysis Step (1):', 1, 'b', steps, [1 2])};
    if strcmp(step, 'estimate')
        steps = {'first-level', 'second-level'};
        stind = spm_input('Analysis Step (2):', '+1', 'b', steps, [1 2]);
    else
        steps = {'classify', 'reconstruct'};
        stind = spm_input('Analysis Step (2):', '+1', 'b', steps, [3 4]);
    end;
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
    file = files(spm_input('Decoding Analysis:', '+1', 'm', list, [1:numel(files)])).name;
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
    what = whats{spm_input('What Exactly to Show:', '+1', 'm', list, [1:numel(list)])};
end;

% Load ITEM.mat file
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


%=========================================================================%
% D I S P L A Y                                                           %
%=========================================================================%

% Get number of sessions
%-------------------------------------------------------------------------%
switch step
    case 'est-1st-lvl', s = numel(GLM1.Sess);
    case 'est-2nd-lvl', s = numel(GLM2.Sess);
    case 'dec-class',   s = numel(ITEM.Sess);
    case 'dec-recon',   s = numel(ITEM.Sess);
end;
csd2 = ceil(s/2);
fsd4 = floor(csd2/2);

% Case: estimate first-level model
%-------------------------------------------------------------------------%
if strcmp(step, 'est-1st-lvl')
    
    % Display design matrices
    %---------------------------------------------------------------------%
    if strcmp(what,'des-mat')
        figure('Name', 'ITEM_est_1st_lvl: design matrices', 'Color', [1 1 1], 'Position', [50 50 1600 900]);
        colormap gray;
        for h = 1:s
            % extract matrices
            X  = GLM1.Sess(h).X*GLM1.Sess(h).T;
            Xt = GLM1.Sess(h).X;
            T  = GLM1.Sess(h).T;
            % standard design matrix
            subplot(csd2, 7, (mod(h,csd2)+csd2*(mod(h,csd2)==0)-1)*7+1+4*(h>csd2));
            imagesc(X);
            caxis([-max(max(abs(X))), +max(max(abs(X)))]);
            axis off;
            title(sprintf('S%d: X = X_t T', h), 'FontSize', 12);
            % first-level (scan-wise) design matrix
            subplot(csd2, 7, (mod(h,csd2)+csd2*(mod(h,csd2)==0)-1)*7+2+4*(h>csd2));
            imagesc(Xt);
            caxis([-1, +1]);
            axis off;
            title(sprintf('X_t [%d x %d]', GLM1.n(h), GLM1.tr(h)), 'FontSize', 12);
            % second-level (trial-wise) design matrix
            subplot(csd2, 7, (mod(h,csd2)+csd2*(mod(h,csd2)==0)-1)*7+3+4*(h>csd2));
            imagesc(T);
            caxis([-1, +1]);
            axis off;
            title(sprintf('T_{ } [%d x %d]', GLM1.tr(h), GLM1.pr(h)), 'FontSize', 12);
            % delete matrices
            clear X Xt T
        end;
    end;
    
    % Display induced covariance
    %---------------------------------------------------------------------%
    if strcmp(what,'ind-cov')
        figure('Name', 'ITEM_est_1st_lvl: induced covariance', 'Color', [1 1 1], 'Position', [50 50 1600 900]);
        colormap jet;
        for h = 1:s
            % extract matrices
            Xt = GLM1.Sess(h).X;
            V  = GLM1.Sess(h).V;
            Ui = inv(GLM1.Sess(h).U);
            um = max(max(abs(Ui(1:GLM1.t(h),1:GLM1.t(h)))));
            % trial-wise design matrix
            subplot(csd2, 7, (mod(h,csd2)+csd2*(mod(h,csd2)==0)-1)*7+1+4*(h>csd2));
            imagesc(Xt);
            caxis([-1, +1]);
            axis off;
            title(sprintf('S%d: X_t [%d x %d]', h, GLM1.n(h), GLM1.tr(h)), 'FontSize', 12);
            % temporal covariance matrix
            subplot(csd2, 7, (mod(h,csd2)+csd2*(mod(h,csd2)==0)-1)*7+2+4*(h>csd2));
            imagesc(V);
            caxis([0, +1]);
            axis square off;
            title(sprintf('V_{ } [%d x %d]', GLM1.n(h), GLM1.n(h)), 'FontSize', 12);
            % variance-covariance matrix
            subplot(csd2, 7, (mod(h,csd2)+csd2*(mod(h,csd2)==0)-1)*7+3+4*(h>csd2));
            imagesc(Ui);
            caxis([-um, +um]);
            axis square off;
            title('U^{-1} = X_t^T V^{-1} X_t', 'FontSize', 12);
            % delete matrices
            clear Xt V Ui um
        end;
    end;
    
end;

% Case: estimate second-level model
%-------------------------------------------------------------------------%
if strcmp(step, 'est-2nd-lvl')
    
    % Display covariance matrices
    %---------------------------------------------------------------------%
    if strcmp(what,'cov-mat')
        figure('Name', 'ITEM_est_2nd_lvl: covariance matrices', 'Color', [1 1 1], 'Position', [50 50 1600 900]);
        colormap jet;
        for h = 1:s
            % extract matrices
            load(strcat(SPM.swd,'/','ITEM_est_1st_lvl','/','GLM1.mat'));
            t  = GLM1.t; clear GLM1
            V  = GLM2.Sess(h).V;
            Q1 = GLM2.Sess(h).Q{1};
            Q2 = GLM2.Sess(h).Q{2};
            vm = max(max(abs(V(1:t(h),1:t(h)))));
            um = max(max(abs(Q2(1:t(h),1:t(h)))));
            % covariance matrix
            subplot(csd2, 7, (mod(h,csd2)+csd2*(mod(h,csd2)==0)-1)*7+1+4*(h>csd2));
            imagesc(V);
            caxis([-vm, +vm]);
            axis square off;
            title(sprintf('S%d: V_{ } [%d x %d]', h, size(V)), 'FontSize', 12);
            % natural covariance
            subplot(csd2, 7, (mod(h,csd2)+csd2*(mod(h,csd2)==0)-1)*7+2+4*(h>csd2));
            imagesc(Q1);
            caxis([0, 1]);
            axis square off;
            title(sprintf('= %3.3f x I_t', GLM2.Sess(h).s2(1)), 'FontSize', 12);
            % induced covariance
            subplot(csd2, 7, (mod(h,csd2)+csd2*(mod(h,csd2)==0)-1)*7+3+4*(h>csd2));
            imagesc(Q2);
            caxis([-um, +um]);
            axis square off;
            title(sprintf('+ %3.3f x U_{ }', GLM2.Sess(h).s2(2)), 'FontSize', 12);
            % delete matrices
            clear V Q1 Q2 Vi Ui vm um
        end;
    end;
    
    % Display covariance components
    %---------------------------------------------------------------------%
    if strcmp(what,'cov-comp')
        figure('Name', 'ITEM_est_2nd_lvl: covariance components', 'Color', [1 1 1], 'Position', [50 50 1600 900]);
        colormap jet;
        % extract matrices
        load(strcat(SPM.swd,'/','ITEM_est_1st_lvl','/','GLM1.mat'));
        t  = GLM1.t; clear GLM1
        Q1 = GLM2.Sess(1).Q{1};
        Q2 = GLM2.Sess(1).Q{2};
        um = max(max(abs(Q2(1:t(1),1:t(1)))));
        % natural covariance
        subplot(1,2,1);
        imagesc(Q1);
        caxis([0, 1]);
        axis square off;
        title(sprintf('Session 1: Q_1 = I_t [%d x %d]', size(Q1)), 'FontSize', 12);
        % induced covariance
        subplot(1,2,2);
        imagesc(Q2);
        caxis([-um, +um]);
        axis square off;
        title(sprintf('Session 1: Q_2 = U_{ } [%d x %d]', size(Q2)), 'FontSize', 12);
        % delete matrices
        clear Q1 Q2 um
    end;
    
end;

% Cases: decode by classification or reconstruction
%-------------------------------------------------------------------------%
if strcmp(step, 'dec-class') || strcmp(step, 'dec-recon')
    
    % Display inverted model
    %---------------------------------------------------------------------%
    if strcmp(what,'inv-mod')
        figure('Name', sprintf('ITEM_dec_%s: inverted model', step(5:end)), 'Color', [1 1 1], 'Position', [50 50 1600 900]);
        colormap jet;
        for h = 1:s
            % extract matrices
            load(strcat(SPM.swd,'/','ITEM_est_1st_lvl','/','GLM1.mat'));
            t  = GLM1.t; clear GLM1
            Y  = ITEM.Sess(h).Y;
            X  = ITEM.Sess(h).X;
            V  = ITEM.Sess(h).V;
            xm = max(max(abs(X(1:end-1,:))));
            vm = max(max(abs(V(1:t(h),1:t(h)))));
            % data matrix
            subplot(csd2, 7, (mod(h,csd2)+csd2*(mod(h,csd2)==0)-1)*7+1+4*(h>csd2));
            imagesc(Y);
            caxis([-1, +1]);
            axis off;
            title(sprintf('S%d: "data" T [%d x %d]', h, size(Y)), 'FontSize', 12);
            % design matrix
            subplot(csd2, 7, (mod(h,csd2)+csd2*(mod(h,csd2)==0)-1)*7+2+4*(h>csd2));
            imagesc(X);
            caxis([-xm, +xm]);
            axis off;
            title(['"design" \Gamma', sprintf(' [%d x %d]', size(X))], 'FontSize', 12);
            % covariance matrix
            subplot(csd2, 7, (mod(h,csd2)+csd2*(mod(h,csd2)==0)-1)*7+3+4*(h>csd2));
            imagesc(V);
            caxis([-vm, +vm]);
            axis square off;
            title(['"covariance" U', sprintf(' [%d x %d]', size(V))], 'FontSize', 12);
            % delete matrices
            clear Y X V xm vm
        end;
    end;
    
    % Display model inversion
    %---------------------------------------------------------------------%
    if strcmp(what,'mod-inv')
        figure('Name', sprintf('ITEM_dec_%s: model inversion', step(5:end)), 'Color', [1 1 1], 'Position', [50 50 1600 900]);
        colormap gray;
        for h = 1:s
            % extract matrices
            Yp = ITEM.Sess(h).Yp;
            if strcmp(step, 'dec-class'), Yc = ITEM.Sess(h).Yc; end;
            if strcmp(step, 'dec-recon'), Yr = ITEM.Sess(h).Yr; end;
            Yt = ITEM.Sess(h).Yt;
            mp = max(max(abs(Yp(1:end-1,1:end-1))));
            % predicted design matrix
            subplot(csd2, 7, (mod(h,csd2)+csd2*(mod(h,csd2)==0)-1)*7+1+4*(h>csd2));
            imagesc(Yp);
            caxis([-mp, +mp]);
            set(gca,'Box','On');
            set(gca,'XTick',[],'YTick',[]);
            title(sprintf('S%d: predicted DM', h), 'FontSize', 12);
            % classified conditions
            subplot(csd2, 7, (mod(h,csd2)+csd2*(mod(h,csd2)==0)-1)*7+2+4*(h>csd2));
            if strcmp(step, 'dec-class'), imagesc(Yc); end;
            if strcmp(step, 'dec-recon'), imagesc(Yr); end;
            caxis([-1, +1]);
            set(gca,'Box','On');
            set(gca,'XTick',[],'YTick',[]);
            if strcmp(step, 'dec-class'), title('classfied conditions', 'FontSize', 12);    end;
            if strcmp(step, 'dec-recon'), title('reconstructed variables', 'FontSize', 12); end;
            % true conditions
            subplot(csd2, 7, (mod(h,csd2)+csd2*(mod(h,csd2)==0)-1)*7+3+4*(h>csd2));
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
    
end;
    
% Case: decode by classification
%-------------------------------------------------------------------------%
if strcmp(step, 'dec-class')
    
    % Display decoding accuracies
    %---------------------------------------------------------------------%
    if strcmp(what,'dec-acc')
        figure('Name', 'ITEM_dec_class: decoding accuracies', 'Color', [1 1 1], 'Position', [50 50 1600 900]);
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
        figure('Name', 'ITEM_dec_recon: correlation coefficients', 'Color', [1 1 1], 'Position', [50 50 1600 900]);
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