%% Roi-based MVPA for CAN Lab
%
% Load spmT images from each subject, apply ROI mask, calculate 
% correlations between all trials.

if isunix % if we are Hammer, a unix system
    addpath(genpath('/gpfs/group/n/nad12/RSA/Scripts/CoSMoMVPA-master'))
else % if not on unix, assume we are on Anvil
    addpath(genpath('S:\nad12\CoSMoMVPA-master'))
end

% add the functions subfolder to the MATLAB search path
path = fileparts(mfilename('fullpath'));
addpath([path filesep 'functions'])

%% Set analysis parameters
subjects   = {'18y404'  '18y566'  '20y297' '20y415'  '20y439'}; % '20y396' <-- this subjects doesn't have a model run
rois       = {'rLTG_left'}; %
study_path = '/gpfs/group/n/nad12/RSA/Analysis_ret/FAMEret8RSA_hrf'; % path to model

% initalizing z_all and rho_all cell arrays
z_all   = cell(1,length(rois));
rho_all = cell(1,length(rois));

for ss = 1:length(subjects)
 
    % Edit the SPM.mat file to use paths here on Hammer
    if isunix % only execute if we are on a Unix system (i.e., Hammer)
        spm_changepath(fullfile(study_path, subjects{ss}, 'SPM.mat'), 'S:\nad12\FAME8', '/gpfs/group/n/nad12/RSA')
        spm_changepath(fullfile(study_path, subjects{ss}, 'SPM.mat'), '\', '/')
    end
    
    % This subjects data_path and spm_path
    data_path  = fullfile(study_path, subjects{ss});
    spm_path   = fullfile(data_path, 'SPM.mat');
    
    % Create a contrast for each regressor in this subjects SPM.mat file
    contrasts(spm_path);
    
    %-- Concatenate the spmT images
    
    % Identify the spmT images
    concat_spmT_fn  = fullfile(data_path, 'concatenated_spmT_images.nii');
    spmT_volumes    = cellstr(spm_select('FPList', data_path, 'spmT.*\.nii'));

    if ~exist(concat_spmT_fn, 'file')
        matlabbatch      = set_3D_to_4D(volumes, 'concatenate_contrasts.nii');
        spm_jobman('run', matlabbatch)
    end

  for rr = 1:length(rois)

    % ROI label
    roi_label = rois{rr}; % name of ROI mask within which we will compare pattern similarity
    
    % full path to ROI mask
    mask_fn  = fullfile(study_path, [roi_label '.nii']);

    % load concatenated spmT images
    ds_all  = cosmo_fmri_dataset(concat_spmT_fn, ...
                                          'mask', mask_fn, ...
                                          'targets', (1:length(spmT_volumes))', ...
                                          'chunks', 1);
                                      
    % Add Labels
    SPM = [];
    load(spm_path);
    ds_all.sa.labels = {SPM.Xx.name}';
                                      
    % cosmo check to make sure data in right format
    cosmo_check_dataset(ds_all);

    % get the samples. The samples are the spmT values from all of the
    % voxels for all of the trials in the current ROI
    all_ds_samples = ds_all.samples;

    % compute correlation values between the all trials, resulting
    % in a nTrials x nTrials matrix, where cell of the matrix represents
    % the correlations between the voxels for each pair of trials.
    rho = cosmo_corr(all_ds_samples');

    % Correlations are limited between -1 and +1, thus they cannot be normally
    % distributed. To make these correlations more 'normal', apply a Fisher
    % transformation and store this in a variable 'z'
    z = atanh(rho);

    %% store and save results
    % path to save results
    output_path = fullfile(study_path, subjects{ss}, 'RSA_Individual_Results');

    % create the output path if it doesn't already exist
    if ~exist(output_path, 'dir')
        mkdir(output_path)
    end
    
    %% Write rho matrix to Excel
    filename = [subjects{ss}, '_' roi_label '_rho_matix.xlsx'];
    xlswrite(fullfile(output_path, filename), rho)

    %% Write z matrix to Excel
    filename = [subjects{ss}, '_' roi_label '_z_matrix.xlsx'];
    xlswrite(fullfile(output_path, filename), z)

    %% Store result in a cell array for later calculations
    z_all{rr}   = cat(3, z_all{rr}, z); % a nTrials x nTrials x nSubjects 3-D Matrix
    rho_all{rr} = cat(3, rho_all{rr}, rho); % a nTrials x nTrials x nSubjects 3-D Matrix
    
  end

end

%% Compute T-Statistics
% Do a t-test on the average z-value for each combination of conditions to
% see if the z-value is signifcantly greater than zero (which would
% indicate no correlation)

% initalize an empty cell array to hold the stats matrices
stats_all = cell(1,length(rois));

% output path
outpath = fullfile(study_path, 'RSA_Group_Results');

% create the output path if it doesn't already exist
if ~exist(outpath, 'dir')
    mkdir(outpath)
end

% for each roi..
for rr = 1:length(rois)
    
    % how many rows and columns are there in the z matrix?
    [rows, cols] = size(z_all{rr});
    
    % for each row; for each column...
    for crow = 1:rows
        for ccol = 1:cols
            
            % do a one-sample t-test across subjects to see if the average
            % z value is signficnatly greater than 0 for this pair of
            % trials in this ROI
            [h,p,ci,stats]            = ttest(z_all{rr}(crow,ccol,:));
            stats_all{rr}(crow, ccol) = stats.tstat;
            
        end
    end
    
    % Save the t statistics to a file
    filename = [roi_label '_onesampleT_results.xlsx'];
    xlswrite(fullfile(outpath, filename), stats_all{rr})
    
end

%% Display Averaged Rho Matrices

% plot an image of the correlation matrix averaged over
% participants (one for each roi)

% allocate space for axis handles, so that later all plots can be set to
% have the same color limits
ax_handles = zeros(length(rois),1);
col_limits = zeros(length(rois),2);

for rr = 1:length(rois)
    
    % create a new figure
    figure

    % store axis handle for current figure
    ax_handles(rr) = gca;

    % compute the average correlation matrix
    rho_mean_matrix  = mean(rho_all{rr}, 3);
    
    % visualize the rho_mean_matrix using imagesc
    imagesc(rho_mean_matrix);

    % set labels, taken from the SPM.mat file
    set(gca, 'xtick', 1:numel(ds_all.sa.labels), 'xticklabel', ds_all.sa.labels)
    set(gca, 'ytick', 1:numel(ds_all.sa.labels), 'yticklabel', ds_all.sa.labels)

    % colorbar
    colorbar;
    
    % title
    desc=sprintf(['Average correlations among trials across subjects '...
                    'in mask ''%s'''], rois{rr});
    title(desc)

    % store color limits so that all graphs are on the same scale for easy
    % visual comparison
    col_limits(rr,:) = get(gca, 'clim');
end

% give all figures the same color limits such that correlations can be
% compared visually
set(ax_handles, 'clim', [min(col_limits(:)), max(col_limits(:))])