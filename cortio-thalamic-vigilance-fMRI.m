%%%
% GOAL: Analysis of functional connectivity from Thalamus (subject specific)
%       parcellation --> to cortical regions. 
% DATA: fMRI timeseries from patients with TLE and healthy controls. 


% Select hi/low vigilance window from each subject. Compare using 4 plot of
% connectivity to all brain regions contained in each of the networks, then
% look specifically in the language regions, and those regions that are
% essential for cognitive impairments, patients with Epilespy suffer from. 

% CEB
% last updated: 11/22/2023



%% 
close all
clear all
clc
restoredefaultpath
addpath Z:\shared_toolboxes\MatlabProgressBar\
addpath Z:\shared_toolboxes\Derek_functions\
addpath Z:\Camden\CEB-helper-functions\
struct_folder = "Z:\000_Data\fMRI\spm8_new_preprocessed_data\data\parcellated_and_connectivity_fmri_data\subj_specific_dkatlas_nbm_AANv2_thalamus_2minWindow_2sStride";
addpath Z:\Camden\parcellations\dk_atlas\
load("thtlamus.mat")

[pat_fpaths, pat_names] = get_files_in_dir(struct_folder, "*pat*", 1);
[con_fpaths, con_names] = get_files_in_dir(struct_folder, "*con*", 1);


% load in CSV 
% Not entirely sure where the CSV is at the momment...

X_pat = nan(2,2,109,109,size(pat_fpaths,1)); % res1/rest2 x hi/low vig x n x n x subject#
X_con = nan(2,2,109,109,size(con_fpaths,1)); % res1/rest2 x hi/low vig x n x n x subject#
dsNotCnt_vig_pat = ones(size(pat_fpaths,1),1);
dsNotCnt_vig_con = ones(size(con_fpaths,1),1);

%%
%%% PATIENT fMRI DOWNLOAD %%%

for pat = progress(1:size(pat_fpaths, 1))

    contained_mat_dir = dir(fullfile(pat_fpaths(pat),'*.m*'));
    F = importdata([contained_mat_dir.folder + "\" + contained_mat_dir.name]);

    % CHECK FOR VIGILANCE INFO IN STRUCT
    if isfield(F,"rest1_minimum_vigilance_window")
    else
        dsNotCnt_vig_pat(pat) = 0;
        continue
    end

    % CHECK FOR ISPI/CONTRA FLIPPED AGE-CORRECTED
    if isfield(F, "age_corrected_rest1_zscore_ic")
    else
        dsNotCnt_vig_pat(pat) = 2;
        continue
    end
    X_pat(1,1,:,:,pat) = F.age_corrected_rest1_zscore_ic(:,:,F.rest1_minimum_vigilance_window);
    X_pat(1,2,:,:,pat) = F.age_corrected_rest1_zscore_ic(:,:,F.rest1_maximum_vigilance_window);
    X_pat(2,1,:,:,pat) = F.age_corrected_rest2_zscore_ic(:,:,F.rest2_minimum_vigilance_window);
    X_pat(2,2,:,:,pat) = F.age_corrected_rest2_zscore_ic(:,:,F.rest2_maximum_vigilance_window);
end

% OMIT all data that doesn't pass the "if" loop
X_pat = X_pat(:,:,:,:,dsNotCnt_vig_pat == 1);

%%

%%% CONTROL fMRI DOWNLOAD %%%
% INITIALIZE 75% FLIP VECTOR
random_flip = rand(1, size(con_fpaths, 1)); 
random_flip = (random_flip <= 0.75);

for con = progress(1:size(con_fpaths, 1))

    contained_mat_dir = dir(fullfile(con_fpaths(con),'*.m*'));
    S = importdata([contained_mat_dir.folder + "\" + contained_mat_dir.name]);

    % CHECK FOR VIGILANCE 
    if isfield(S,"rest1_minimum_vigilance_window")
    else
        dsNotCnt_vig_con(con) = 0;
        continue
    end

    % CHECK FOR AGE_CORRECTION & FLIPPING
    if isfield(S,"age_corrected_rest1_zscore_flipped")
    else
        dsNotCnt_vig_con(con) = 2;
        continue
    end
    
    % RANDOMLY FLIP THE CONTORLS SUBJECTS RIGHT/LEFT HEMISPHERES
    if     random_flip(con) == 0
        bip  = S.age_corrected_rest1_zscore;
        bip2 = S.age_corrected_rest2_zscore;
    elseif random_flip(con) == 1
        bip  = S.age_corrected_rest1_zscore_flipped;
        bip2 = S.age_corrected_rest2_zscore_flipped;
    end
   
    X_con(1,1,:,:,con) =  bip(:,:,S.rest1_minimum_vigilance_window);
    X_con(1,2,:,:,con) =  bip(:,:,S.rest1_maximum_vigilance_window);
    X_con(2,1,:,:,con) = bip2(:,:,S.rest2_minimum_vigilance_window);
    X_con(2,2,:,:,con) = bip2(:,:,S.rest2_maximum_vigilance_window);
end

% OMIT all data that doesn't pass the "if" loop
X_con = X_con(:,:,:,:,dsNotCnt_vig_con == 1);

%% 
%%% CREATE THALAMUS NODE DESIGNATIONS INDEX
%%% CREATE ALL OTHER NETWORK DESIGNATIONS

reg_names = S.region_names;

thal_idx = cell(size(ThalsubSpef));

for jj = 1:size(ThalsubSpef)
    thal_idx{jj} = find(contains(reg_names, ThalsubSpef(jj)));     
end
thal_idx = cell2mat(thal_idx); 


% CREATE A CSV OF P-VALUES THAT SHOW THE DIFFERENCES IN PATIENTS AND
% CONTROLS AT THE SAME VIGILANCE STATES FROM  THALAMIC NODES
% TO ALL OTHER NODES.

n = 109; 

cnt_corticothal_pats = squeeze(mean(squeeze(X_pat(:, :, thal_idx, :, :)),3));
cnt_corticothal_cons = squeeze(mean(squeeze(X_con(:, :, thal_idx, :, :)),3));

p_val = nan(2,2,n);

forCSVRest1 = nan(n,2); 
forCSVRest2 = nan(n,2);

for rest = 1:2
    for hilo = 1:2
        for reg = 1:n
            [~, p_val(rest,hilo,reg)] = ttest2(squeeze(cnt_corticothal_pats(rest,hilo,reg,:)), ...
            squeeze(cnt_corticothal_cons(rest,hilo,reg,:))); 
        end
        
    end
    forCSVRest1(:,rest) = squeeze(p_val(1,rest,:))';
    forCSVRest2(:,rest) = squeeze(p_val(2,rest,:))';
end

export_csv = [forCSVRest1, forCSVRest2, reg_names]; 

writematrix(export_csv, 'SubjectSpec-Pvals-Cortio-Thalamus.csv')






