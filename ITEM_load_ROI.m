function [roi_ind] = ITEM_load_ROI(GLM, ROI)
% _
% Load Voxel Indices from Region of Interest Image
% FORMAT [roi_ind] = ITEM_load_ROI(GLM, ROI)
% 
%     GLM     - a structure specifying a trial-wise GLM
%     ROI     - a filepath to a region of interest image
% 
%     roi_ind - a 1 x v vector of ROI voxel indices
% 
% FORMAT [roi_ind] = ITEM_load_ROI(GLM, ROI) obtains the indices of
% voxels that are non-NaN and non-zero in a region of interest image.
% 
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 20/12/2018, 05:45 (V0.1)
%  Last edit: 20/12/2018, 05:45 (V0.1)


% Get number of sessions
%-------------------------------------------------------------------------%
s = numel(GLM.Sess);

% Load ROI image
%-------------------------------------------------------------------------%
roi_hdr = spm_vol(ROI);
roi_dim = roi_hdr.dim;
roi_img = spm_read_vols(roi_hdr);

% Select ROI voxels
%-------------------------------------------------------------------------%
roi_ind = find(~isnan(roi_img) & roi_img~=0);
if numel(roi_ind) > (s-1)*GLM.tr(1)-1
    roi_mat = [roi_ind, roi_img(roi_ind)];
    roi_mat = flipud(sortrows(roi_mat, 2));
    roi_ind = roi_mat(1:(s-1)*GLM.tr(1)-1,1)';
    % NOTE: When the number of voxels in a ROI exceeds the number of data
    % points in the training data of one cross-validation fold, then the
    % ROI is restricted to those voxels with the highest values.
else
    roi_ind = roi_ind';
end;