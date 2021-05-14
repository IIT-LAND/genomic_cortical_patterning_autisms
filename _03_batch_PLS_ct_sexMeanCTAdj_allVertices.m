function result = pls_03_batch_PLS_all_ct_sexMeanCTAdj_allVertices

RUNONSERVER = 1;
if RUNONSERVER
    addpath /media/DATA/mlombardo/ace_smri_gex_pls/code/plscmd;
    rootpath = '/media/DATA/mlombardo/ace_smri_gex_pls';
else
    addpath /Users/mlombardo/Dropbox/matlab/Pls/plscmd/;
    rootpath = '/Users/mlombardo/Dropbox/ACE_sMRI_GEX_PLS/reproAnalysis';
end % if RUNONSERVER

n_perm = 10000;
n_boot = 10000;
n_splithalf = 0;

datapath = fullfile(rootpath,'data','tidy');
MEfiles = {fullfile(datapath,'ASDGood_all_MEdata.txt'), ...
    fullfile(datapath,'ASDPoor_all_MEdata.txt'), ...
    fullfile(datapath,'TD_all_MEdata.txt')};

lh_regions = 149953;
rh_regions = 149926;

fname2save = fullfile(rootpath,'pls','results','plsres_subgrp_all_allVertices_cortical_ct_me_sexMeanCTAdj.mat');


%% read in text files
% read in ME files
for ifile = 1:length(MEfiles)
%     fname = fullfile(datadir,MEfiles{ifile});
    fname = MEfiles{ifile};
    tab2use = readtable(fname);
    ME{ifile} = table2array(tab2use(:,2:end));
    me_subids{ifile} = tab2use(:,1);
end % for ifile

me_names = tab2use.Properties.VariableNames(2:end);

mri_data = readtable(fullfile(datapath,'mri_data_adj_sexMeanCT_all_cortical_ct_allVertices.csv'));
label_data = readtable(fullfile(datapath,'labelData_all_MRI.csv'));
grp_labels = {'Good','Poor','TD'};
for igrp = 1:length(grp_labels)
    mask = ismember(label_data.subgrp,grp_labels{igrp});
    BrainVols{igrp} = table2array(mri_data(mask,3:end));
    brainvols_subids{igrp} = mri_data(mask,2);
end % for igrp

brainreg_names = mri_data.Properties.VariableNames(3:end);

tmp_array = [];
for ime = 1:length(ME)
    tmp_array = [tmp_array;ME{ime}];
end % for ime
MEstacked = tmp_array;

%% Prepare input arguments for pls_analysis.m
datamat_lst = BrainVols; % initialize datamat_lst
for icell = 1:length(datamat_lst)
    num_subj_lst(icell) = size(datamat_lst{icell},1); % initialize num_subj_lst
end % for icell

num_cond = 1; %{ones(1,length(datamat_lst))};  % initialize num_cond

% specify option structure
% option.progress_hdl = [];%( user interface handle )
option.method = 3; %[1] | 2 | 3 | 4 | 5 | 6
option.num_perm = n_perm; %( single non-negative integer )
option.is_struct = 0;%[0] | 1
option.num_split = n_splithalf; %( single non-negative integer )
option.num_boot = n_boot; %( single non-negative integer )
option.clim = 95; %( [95] single number between 0 and 100 )
% option.bscan = ( subset of  1:num_cond )
% option.stacked_designdata = ( 2-D numerical matrix )
option.stacked_behavdata = MEstacked;
% option.meancentering_type = [0] | 1 | 2 | 3
option.cormode = 0; %[0] | 2 | 4 | 6
option.boot_type = 'strat'; %['strat'] | 'nonstrat'


%% run pls_analysis.m
result = pls_analysis(datamat_lst, num_subj_lst, num_cond, option);
result.me_names = me_names;

%% Compute percentage of cross-block covariance
result.crossblockCovPercent = result.s.^2/sum(result.s.^2);

%% fix p-values
result.perm_result.sprob = (result.perm_result.sp+1)./(result.perm_result.num_perm+1);

%% find the significant LVs
result.sigLVs = find(result.perm_result.sprob<=0.05);
result.sigLVs_pvals = result.perm_result.sprob(result.sigLVs);
disp(sprintf('Significant LVs: LV %d',result.sigLVs));


%% compute brain BSR
for i = 1:size(result.boot_result.compare_u,2)
    result.brain_bsr(:,i) = result.boot_result.compare_u(:,i)./result.boot_result.u_se(:,i);
end
result.brainreg_names = brainreg_names';

save(fname2save,'result');


%% compute bootstrap CIs
cis2use = {[2.5,97.5]};
for i = 1:length(result.sigLVs)
    LVnum = result.sigLVs(i);
    for j = 1:length(cis2use)
        ci_bounds = cis2use{j};
        clear asd_good_bootres asd_poor_bootres td_bootres asd_good_ci asd_poor_ci td_ci asd_good_corr asd_poor_corr td_corr asd_good asd_poor td all_data;
        ORDERBY = 'TD';
        mod_names_orig = me_names;
        
        mod_names = repmat(mod_names_orig',length(MEfiles),1);
        
        asd_good_idx = 1:length(mod_names_orig);
        asd_poor_idx = (length(mod_names_orig)+1):(length(mod_names_orig)*2);
        td_idx = ((length(mod_names_orig)*2)+1):(length(mod_names_orig)*3);
        
        asd_good_bootres = squeeze(result.boot_result.distrib(asd_good_idx,LVnum,:));
        asd_poor_bootres = squeeze(result.boot_result.distrib(asd_poor_idx,LVnum,:));
        td_bootres = squeeze(result.boot_result.distrib(td_idx,LVnum,:));
        
        asd_good_ci = prctile(asd_good_bootres',ci_bounds)';
        asd_poor_ci = prctile(asd_poor_bootres',ci_bounds)';
        td_ci = prctile(td_bootres',ci_bounds)';
        
        asd_good_corr = result.boot_result.orig_corr(asd_good_idx,LVnum);
        asd_poor_corr = result.boot_result.orig_corr(asd_poor_idx,LVnum);
        td_corr = result.boot_result.orig_corr(td_idx,LVnum);
        
        if strcmp(ORDERBY,'TD')
            [idx, plot_order] = sort(td_corr,'ascend');
            plot_order(plot_order) = 1:length(mod_names_orig);
        elseif strcmp(ORDERBY,'Good')
            [idx, plot_order] = sort(asd_good_corr,'ascend');
            plot_order(plot_order) = 1:length(mod_names_orig);
        elseif strcmp(ORDERBY,'Poor')
            [idx, plot_order] = sort(asd_poor_corr,'ascend');
            plot_order(plot_order) = 1:length(mod_names_orig);
        end
        
        asd_good = [asd_good_corr asd_good_ci];
        asd_poor = [asd_poor_corr asd_poor_ci];
        td = [td_corr td_ci];
        
        all_data = [asd_good; asd_poor; td];
        all_data(:,4) = (sign(all_data(:,1)) == sign(all_data(:,2))) & (sign(all_data(:,1)) == sign(all_data(:,3)));
        all_data(:,5) = repmat(plot_order,3,1);
        
        group_labels = [repmat({'Good'},length(mod_names_orig),1);repmat({'Poor'},length(mod_names_orig),1);repmat({'TD'},length(mod_names_orig),1)];
        
        tab2write = cell2table([group_labels, mod_names num2cell(all_data)], ...
            'VariableNames',{'Grp','ModName','corr','lo_lim','up_lim','nonzero','plot_order'});
        
        [f_path, f_name, f_ext] = fileparts(fname2save);
        file2save = fullfile(f_path,sprintf('plsres_all_ct_allVertices_sexMeanCTAdj_MEcorr_bootCI4plotting_LV%d_ci%d.csv',LVnum,ci_bounds(2)-ci_bounds(1)));
        writetable(tab2write,file2save,'FileType','text','delimiter',',');
    end % for j
end % for i

%% output brain BSR
for i = 1:length(result.sigLVs)
    bsr_names{i} = sprintf('BSR_LV%d',result.sigLVs(i));
end
feature_names = [repmat({'Cortical Thickness'},lh_regions,1); repmat({'Cortical Thickness'},rh_regions,1)];
hemi_names = [repmat({'LH'},lh_regions,1);repmat({'RH'},rh_regions,1)];
tab2write = cell2table([brainreg_names', feature_names, hemi_names, num2cell(result.brain_bsr(:,result.sigLVs))],'VariableNames',[{'brainreg','feature','hemisphere'},bsr_names]);
writetable(tab2write, fullfile(f_path,'plsres_all_ct_allVertices_sexMeanCTAdj_brainBSR4plotting.csv'),'FileType','text','Delimiter',',');

%% save result
save(fname2save,'result');







%% Correct file for bad vertices
RUNONSERVER = 1;
if RUNONSERVER
    rootpath = '/media/DATA/mlombardo/ace_smri_gex_pls';
else
    rootpath = '/Users/mlombardo/Dropbox/ACE_sMRI_GEX_PLS/reproAnalysis';
end
datapath = fullfile(rootpath,'data','tidy');
plspath = fullfile(rootpath,'pls','results');

%%
feature = 'ct';
n_regions = 163842;

brain_bsr_tab = readtable(fullfile(plspath,sprintf('plsres_all_%s_allVertices_sexMeanCTAdj_brainBSR4plotting.csv',feature)),'delimiter',',');
new_bsr_tab = readtable(fullfile(datapath,sprintf('AllGroups_all_cortical_%s_badVertices.csv',feature)),'delimiter',',');

new_bsr_tab.feature = repmat({'Cortical Thickness'},length(new_bsr_tab.vertex_label),1);
new_bsr_tab.hemisphere = [repmat({'LH'}, n_regions,1);repmat({'RH'}, n_regions,1)];
new_bsr_tab.BSR_LV1 = zeros(length(new_bsr_tab.vertex_label),1);

for i = 1:length(new_bsr_tab.vertex_label)
    disp(i);
    mask = ismember(brain_bsr_tab.brainreg, new_bsr_tab.vertex_label(i));
    if sum(mask)~=0
        new_bsr_tab.BSR_LV1(i) = brain_bsr_tab.BSR_LV1(mask);
    end % if
end % for

writetable(new_bsr_tab, fullfile(plspath,sprintf('plsres_all_%s_allVerticesGood+Bad_sexMeanCTAdj_brainBSR4plotting.csv',feature)),'FileType','text','Delimiter',',');



brain_bsr_tab = readtable(fullfile(plspath,sprintf('plsres_all_%s_allVertices_sexMeanCTAdj_brainBSR4plotting.csv',feature)),'delimiter',',');
new_bsr_tab = readtable(fullfile(datapath,sprintf('AllGroups_all_cortical_%s_badVertices.csv',feature)),'delimiter',',');

new_bsr_tab.feature = repmat({'Cortical Thickness'},length(new_bsr_tab.vertex_label),1);
new_bsr_tab.hemisphere = [repmat({'LH'}, n_regions,1);repmat({'RH'}, n_regions,1)];
new_bsr_tab.BSR_LV2 = zeros(length(new_bsr_tab.vertex_label),1);

for i = 1:length(new_bsr_tab.vertex_label)
    disp(i);
    mask = ismember(brain_bsr_tab.brainreg, new_bsr_tab.vertex_label(i));
    if sum(mask)~=0
        new_bsr_tab.BSR_LV2(i) = brain_bsr_tab.BSR_LV2(mask);
    end % if
end % for

writetable(new_bsr_tab, fullfile(plspath,sprintf('plsres_all_%s_allVerticesGood+Bad_sexMeanCTAdj_brainBSR4plotting_LV2.csv',feature)),'FileType','text','Delimiter',',');

end % function batch_PLS




