%% Roi-based MVPA for single subject (run_split_half_correlations_single_sub)
%
% Load t-stat data from one subject, apply 'vt' mask, compute difference
% of (fisher-transformed) between on- and off diagonal split-half
% correlation values.
%
% #   For CoSMoMVPA's copyright information and license terms,   #
% #   see the COPYING file distributed with CoSMoMVPA.           #

% Preliminary
% clc
% clear all
addpath(genpath('S:\nad12\CoSMoMVPA-master'))
%addpath(genpath('/gpfs/group/n/nad12/RSA/Scripts/CoSMoMVPA-master'))

%% Set analysis parameters
rois = {'rLTG_left', 'rHC_left', 'rSMA_bilat'}; 
subjects = {
    {'18y404','18y566','20y297','20y415','20y439','20y441','20y444','20y455','21y299','21y437','21y521','21y534','22y422','23y452','23y546','25y543'}
    %{'67o136','67o153','67o178','69o144','69o277','70o118','70o316','71o152','71o193','72o164','73o165','75o320','76o120','76o162','79o108','790117','79o279','80o121','80o128','81o125','81o312','83o197'}
}; 
study_path = 'S:\nad12\FAME8\Analysis_ret\FAMEret8RSA_hrf';

for rr = 1:length(rois)
    roi_label = rois{rr};
    filename = ['Compiled_z_', roi_label, '.xlsx'];

    for sg = 1:length(subjects) %subject group
        for ss = 1:length(subjects{sg}) %subject num in group
            %Path toexcel with correlation matrix   
            data_path = fullfile(study_path, subjects{sg}{ss}, 'RSA_Results', ['RSAtest_' subjects{sg}{ss} '_' roi_label '_z_.xlsx']);

            %Determine target row and sheet in excel file
            rownum = num2str(ss);
            sheet  = strcat('Sheet', num2str(sg));

            %Gather correlation values and write to file
            %NOTE: will overwrite existing data in excel if run previously,
            %but won't delete data if more Ss were run previously
            
            %Sheet1 will be group1 (YA), Sheet2 will be grp2 (OA)
            corr1 = xlsread(data_path, 'Sheet1', 'B1');
            xlswrite(fullfile(study_path, filename), corr1, sheet, strcat('A', rownum));

            corr2 = xlsread(data_path, 'Sheet1', 'C1');
            xlswrite(fullfile(study_path, filename), corr2, sheet, strcat('B', rownum));

            corr3 = xlsread(data_path, 'Sheet1', 'D1');
            xlswrite(fullfile(study_path, filename), corr3, sheet, strcat('C', rownum));
            
            corr4 = xlsread(data_path, 'Sheet1', 'C2');
            xlswrite(fullfile(study_path, filename), corr4, sheet, strcat('D', rownum));
            
            corr5 = xlsread(data_path, 'Sheet1', 'D2');
            xlswrite(fullfile(study_path, filename), corr5, sheet, strcat('E', rownum));
            
            corr6 = xlsread(data_path, 'Sheet1', 'D3');
            xlswrite(fullfile(study_path, filename), corr6, sheet, strcat('F', rownum));
        end
    end
end