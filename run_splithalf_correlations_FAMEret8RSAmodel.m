%% Roi-based MVPA for single subject (run_split_half_correlations_single_sub)
%
% Load beta image data from one subject, apply 'vt' mask, compute difference
% of (fisher-transformed) between on- and off diagonal split-half
% correlation values.
%
% #   For CoSMoMVPA's copyright information and license terms,   #
% #   see the COPYING file distributed with CoSMoMVPA.           #

if isunix % if we are Hammer, a unix system
    addpath(genpath('/gpfs/group/n/nad12/RSA/Scripts/CoSMoMVPA-master'))
else % if not on unix, assume we are on Anvil
    addpath(genpath('S:\nad12\CoSMoMVPA-master'))
end

%% Set analysis parameters
subjects   = {'18y404','18y566','20y297'};%'20y415','20y439','20y441','20y444','20y455','21y299','21y437','21y521','21y534','22y422','23y452','23y546','25y543'};%'67o153','67o178','69o144','70o118','71o152','71o193','72o164','73o165','76o120','76o162','78o113','79o108','80o121','80o128','81o125','81o312','83o197'}; 
%Ss to add after Model is ran (02/10/17 ICT) '79o117','70o316','75o320','69o277','79o279','67o136'
rois       = {'rLTG_left'}; %add in leftHC and right HC for starters
study_path = 'S:\nad12\FAME8\Analysis_ret\FAMEret8RSA_hrf';

% Edit the SPM.mat file to use paths here on Hammer
%spm_changepath(fullfile(study_path, subjects{1}, 'SPM.mat'), 'S:\nad12\FAME8', '/gpfs/group/n/nad12/RSA')
%spm_changepath(fullfile(study_path, subjects{1}, 'SPM.mat'), '\', '/')

% initalizing the sum of weighted zs all array
sum_weighted_zs_all = cell(1,length(rois));

for ss = 1:length(subjects)
 
    %% Computations
    data_path  = fullfile(study_path, subjects{ss});
    
     % Average Betas
   average_betas(subjects{ss}, study_path, 'HREC', '.*HREC.*');
   average_betas(subjects{ss}, study_path, 'HFAM', '.*HFAM.*');
   average_betas(subjects{ss}, study_path, 'FAREC', '.*FAREC.*');
   average_betas(subjects{ss}, study_path, 'FAFAM', '.*FAFAM.*');
    
  for rr = 1:length(rois)

    roi_label = rois{rr}; % name of ROI mask used for running correlations  

    % the averaged betas
    HREC     = fullfile(data_path, 'average_beta_HREC.nii');
    HFAM     = fullfile(data_path, 'average_beta_HFAM.nii');
    FAREC    = fullfile(data_path, 'average_beta_FAREC.nii');
    FAFAM    = fullfile(data_path, 'average_beta_FAFAM.nii');
    
    mask_fn  = fullfile(study_path, [roi_label '.nii']); %second half of mask name

    % load two halves as CoSMoMVPA dataset structs.
    % Chunks = Runs  Targets = trial type conditions
    HREC_ds  = cosmo_fmri_dataset(HREC,'mask',mask_fn,... %encoding... run 1 (chunk) Target(rem) will be identified as 1
                                         'targets',(1)',... %ex: Rem
                                         'chunks',(1)); %ex: encoding (avg. of all enc runs)

    HFAM_ds  = cosmo_fmri_dataset(HFAM,'mask',mask_fn,...
                                         'targets',(2)',... %ex. Know
                                         'chunks',(1)); %encoding (avg. of all enc runs)

    FAREC_ds = cosmo_fmri_dataset(FAREC,'mask',mask_fn,...
                                         'targets',(3)',... %New
                                         'chunks',(1)); %ret
                                     
    FAFAM_ds = cosmo_fmri_dataset(FAFAM,'mask',mask_fn,...
                                         'targets',(4)',... %New
                                         'chunks',(1)); %ret


   

    % Combine files at encoding and retrieval to create two files (i.e.,
    % stacking)
    % make sure all ds_* changed from here on
    ds_all   = cosmo_stack({HREC_ds, HFAM_ds, FAREC_ds, FAFAM_ds});

    % Data set labels
    ds_all.sa.labels = {'HREC';'HFAM';'FAREC';'FAFAM'}; 
    
    % cosmo fxn to make sure data in right format
    cosmo_check_dataset(ds_all);
    
    % Some sanity checks to ensure that the data has matching features (voxels)
    % and matching targets (conditions)
   %  assert(isequal(Target_ds.fa,Lure_ds.fa)); %,New_ds.fa
   %  assert(isequal(Target_ds.sa.targets,Lure_ds.sa.targets)); %,New_ds.sa.targets

    % change if you change ds_* above
    nClasses = numel(ds_all.sa.labels);  %%why isn't this pulling all labels?(NAD10.31.16)

    % get the sample data - samples are the correlations being ran on each
    % voxel (a bunch of numbers)
    all_ds_samples = ds_all.samples;

    % compute all correlation values between the two halves, resulting
    % in a 6x6 matrix. Store this matrix in a variable 'rho'.
    % Hint: use cosmo_corr (or builtin corr, if the matlab stats toolbox
    %       is available) after transposing the samples in the two halves.
    % >@@>
    rho = cosmo_corr(all_ds_samples');
    % <@@<

    % Correlations are limited between -1 and +1, thus they cannot be normally
    % distributed. To make these correlations more 'normal', apply a Fisher
    % transformation and store this in a variable 'z'
    % (hint: use atanh).
    % >@@>
    z = atanh(rho);
    % <@@<

    % <@@<

    % Set up a contrast matrix to test whether the element in the diagonal
    % (i.e. a within category correlation) is higher than the average of all
    % other elements in the same row (i.e. the average between-category
    % correlations). For testing the split half correlation of n classes one
    % has an n x n matrix (here, n=6).
    %
    % To compute the difference between the average of the on-diagonal and the
    % average of the off-diagonal elements, consider that there are
    % n on-diagonal elements and n*(n-1) off-diagonal elements.
    % Therefore, set
    % - the on-diagonal elements to 1/n           [positive]
    % - the off-diagonal elements to -1/(n*(n-1)) [negative]
    % This results in a contrast matrix with weights for each element in
    % the correlation matrix, with positive and equal values on the diagonal,
    % negative and equal values off the diagonal, and a mean value of zero.
    %
    % Under the null hypothesis one would expect no difference between the
    % average on the on- and off-diagonal, hence correlations weighted by the
    % contrast matrix has an expected mean of zero. A postive value for
    % the weighted correlations would indicate more similar patterns for
    % patterns in the same condition (across the two halves) than in different
    % conditions.

    % Set the contrast matrix as described above and assign it to a variable
    % named 'contrast_matrix'
    % >@@>
    contrast_matrix = (eye(nClasses)-1/nClasses)/(nClasses-1);

    % alternative solution
    contrast_matrix_alt = zeros(nClasses,nClasses);
    for k = 1:nClasses
        for j = 1:nClasses
            if k == j
                value = 1/nClasses;
            else
                value = -1/(nClasses*(nClasses-1));
            end
            contrast_matrix_alt(k,j) = value;
        end
    end

    % <@@<

    % sanity check: ensure the matrix has a sum of zero
    if abs(sum(contrast_matrix(:)))>1e-14
        error('illegal contrast matrix: it must have a sum of zero');
    end

    % Weigh the values in the matrix 'z' by those in the contrast_matrix
    % and then average them (hint: use the '.*' operator for element-wise
    % multiplication).
    % Store the result in a variable 'weighted_z'.
    % >@@>
    weighted_z = z .* contrast_matrix;
    % <@@<

    % Compute the sum of all values in 'weighted_z', and store the result in
    % 'sum_weighted_z'.
    % >@@>
    sum_weighted_z = sum(weighted_z(:)); %Expected value under H0 is 0
    % <@@<

    %% store and save results
    output_path = fullfile(study_path, subjects{ss}, 'RSA_Results');

    if ~exist(output_path, 'dir')
        mkdir(output_path)
    end
    
    %% Write r matrix to Excel
    filename = ['RSAtest_', subjects{ss}, '_' roi_label '_rho_.xlsx'];
    rho %#ok<*NOPTS>
    xlswrite(fullfile(output_path, filename), rho)

    %% Write z matrix to Excel
    filename = ['RSAtest_', subjects{ss}, '_' roi_label '_z_.xlsx'];
    z
    xlswrite(fullfile(output_path, filename), z)

    %% Write wieghted_z matrix to excel
    filename = ['RSAtest_', subjects{ss}, '_' roi_label '_wieghted_z_.xlsx'];
    weighted_z
    xlswrite(fullfile(output_path, filename), weighted_z)

    %% Write sum of wieghted_z matrix to excel
    filename = ['RSAtest_', subjects{ss}, '_' roi_label '_sum_weighted_z_.xlsx'];
    sum_weighted_z
    xlswrite(fullfile(output_path, filename), sum_weighted_z)

    %% Store result in a matrix for later stats calculations

    % store the result for this subject in sum_weighted_zs_all
    % (at the i_subj-th position), so that
    % group statistics can be computed
    % >@@>
    sum_weighted_zs_all{rr}(i_subj) = sum_weighted_z;
    % <@@<
    
    
  end

end

%% Compute T-Statistics and Print Results

for rr = 1:length(rois)

    % Using matlab's stat toolbox (if present)
    if cosmo_check_external('@stats',false)

        % >@@>
        [h,p,ci,stats] = ttest(sum_weighted_zs_all{rr});
        fprintf(['correlation difference in %s at group level: '...
            '%.3f +/- %.3f, t_%d=%.3f, p=%.5f (using matlab stats '...
            'toolbox)\n'],...
            rois{i_roi},mean(sum_weighted_zs_all{rr}),...
                        std(sum_weighted_zs_all{rr}),...
            stats.df,stats.tstat,p);
        % <@@<
    else
        fprintf('Matlab stats toolbox not available\n');
    end

end
