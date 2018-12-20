function [X, L] = ITEM_get_des_mat(names, onsets, durations, pmod, orth, R, settings)
% _
% Create Design Matrix for First-Level fMRI Data Analysis
% FORMAT [X, L] = ITEM_get_des_mat(names, onsets, durations, pmod, orth, R, settings)
% 
%     names     - a  1 x c cell array of condition names
%     onsets    - a  1 x c cell array of condition onsets
%     durations - a  1 x c cell array of condition durations
%     pmod      - a  1 x c structure with the following fields:
%     o name    - a  1 x p cell array of modulator names
%     o param   - a  1 x p cell array of modulator values
%     orth      - a  1 x c logical vector indicating orthogonalization
%     R         - an n x r design matrix of multiple regressors (e.g. RPs)
%     settings  - a structure variable with the following fields:
%     o n       - an integer, the number of scans
%     o TR      - a  scalar, the fMRI repetition time
%     o dt      - a  scalar, the microtime resolution
%     o HRF     - a  string, the convolution function (e.g. 'spm_hrf')
%     o conv    - a  string, the convolution mode ('none'/'stand'/'persist')
%     o mc      - a  logical, indicating mean-centering (for R)
%     o RPs     - a  logical, indicating realignment parameters (for R)
% 
%     X         - an n x p first-level fMRI design matrix
%     L         - a  1 x p cell array of regressor labels
% 
% FORMAT [X, L] = ITEM_get_des_mat(...) takes names, onsets, durations,
% pmod, orth and R in well-known SPM format and creates a design matrix
% following the specified settings (see below).
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 14/11/2017, 17:35 (V0.0)
%  Last edit: 22/11/2018, 11:45 (V0.1)


% Set default values
%-------------------------------------------------------------------------%
if isempty(orth) || nargin < 5              % orthogonalization
    orth = num2cell(false(size(names)));
end;
if isempty(R) || nargin < 6                 % additional regressors
    R = [];
end;
if isempty(settings) || nargin < 7          % specification settings
    settings = struct([]);
end;

% Enact default settings
%-------------------------------------------------------------------------%
if ~isempty(settings)
    
    % Timing information
    %---------------------------------------------------------------------%
    if ~isfield(settings,'n'),    settings.n    = 1000;      end;
    if ~isfield(settings,'TR'),   settings.TR   = 2;         end;
    
    % Onset regressors & parametric modulators
    %---------------------------------------------------------------------%
    if ~isfield(settings,'dt'),   settings.dt   = 0.01;      end;
    if ~isfield(settings,'hrf'),  settings.hrf  = 'spm_hrf'; end;
    if ~isfield(settings,'conv'), settings.conv = 'persist'; end;
    
    % Additional regressors
    %---------------------------------------------------------------------%
    if ~isfield(settings,'mc'),   settings.mc   = true;      end;
    if ~isfield(settings,'RPs'),
        if size(R,2) == 6,        settings.RPs  = true;      end;
        if size(R,2) ~= 6,        settings.RPs  = false;     end;
    end;
    
end;

% Specify time vector
%-------------------------------------------------------------------------%
t = [0:settings.dt:(settings.n*settings.TR)];
z = zeros(1,2*numel(t)-1);
TRdt = round(settings.TR/settings.dt);

% Specify, normalize and partition HRF
%-------------------------------------------------------------------------%
if strcmp(settings.hrf,'spm_hrf')
    HRF = spm_hrf(settings.dt)';
else
    eval(strcat('HRF = ',settings.hrf,'(settings.dt);'));
end;
 HRF     = HRF./max(HRF);
[hrf, i] = max(HRF);
 HRF_1st = HRF(1:i);
 HRF_2nd = HRF((i+1):end);

% Preallocate design and labels
%-------------------------------------------------------------------------%
X = [];
L = [];

% Create design matrix
%-------------------------------------------------------------------------%
for i = 1:numel(names)
    
    % Onset regressor
    %---------------------------------------------------------------------%
    o = round(onsets{i}./settings.dt);
    d = round(durations{i}./settings.dt);
    % "standard" convolution
    if strcmp(settings.conv,'stand') || strcmp(settings.conv,'none')
        y = z;
        for k = 1:numel(o)
            y((o(k)+1):(o(k)+d(k))) = 1;
        end;
        x = conv(y,HRF);
        x = x(1:numel(z));
    end;
    % "persistent" convolution
    if strcmp(settings.conv,'persist')
        x = zeros(size(z));
        for k = 1:numel(o)
            y = z;
            HRF_ik = [HRF_1st, hrf*ones(1,d(k)), HRF_2nd];
            y((o(k)+1):(o(k)+numel(HRF_ik))) = HRF_ik;
            x = x + y;
        end;
    end;
    % no convolution
    if strcmp(settings.conv,'none')
        x = y;
    end;
    % add regressor
    X = [X, x(1:TRdt:(settings.n-1)*TRdt+1)'];
    L = [L, names(i)];
    
    % Check for PMs
    %---------------------------------------------------------------------%
    if i <= numel(pmod)
        if ~isempty(pmod(i).name)
            for j = 1:numel(pmod(i).name)

                % Parametric modulator
                %---------------------------------------------------------%
                p = pmod(i).param{j};
                p = p-mean(p);
                % "standard" convolution
                if strcmp(settings.conv,'stand') || strcmp(settings.conv,'none')
                    y = z;
                    for k = 1:numel(o)
                        y((o(k)+1):(o(k)+d(k))) = p(k);
                    end;
                    x = conv(y,HRF);
                    x = x(1:numel(z));
                end;
                % "persistent" convolution
                if strcmp(settings.conv,'persist')
                    x = zeros(size(z));
                    for k = 1:numel(p)
                        y = z;
                        HRF_ijk = p(k)*[HRF_1st, hrf*ones(1,d(k)), HRF_2nd];
                        y((o(k)+1):(o(k)+numel(HRF_ijk))) = HRF_ijk;
                        x = x + y;
                    end;
                end;
                % no convolution
                if strcmp(settings.conv,'none')
                    x = y;
                end;
                % add modulator
                X = [X, x(1:TRdt:(settings.n-1)*TRdt+1)'];
                L = [L, {sprintf('%s x %s',names{i}, pmod(i).name{j})}];

            end;
        end;
    end;
    
end;

% Add additional regressors
%-------------------------------------------------------------------------%
if ~isempty(R)
    if settings.mc
        R = R - repmat(mean(R),[settings.n 1]);
    end;
    if settings.RPs
        R(:,end-2:end) = R(:,end-2:end) * (180/pi);
    end;
    for i = 1:size(R,2)
        X = [X, R(:,i)];
        L = [L, {strcat('R',num2str(i))}];
    end;
end;