library(easypackages)
libraries("here","fsbrain","freesurferformats","matlabr","RColorBrewer")

codedir = here("code")
plsdir = here("pls","plots")
subjects_dir = '/Applications/freesurfer/subjects/'

#===============================================================================
RUNMATLAB = TRUE
if (RUNMATLAB){
  code2run = sprintf("cd %s; correct_AllVertices_brainBSRfile; make_vertexwise_mgh_maps4_BSR_NEW;",codedir)
  res = run_matlab_code(code2run)
}



# #===============================================================================
# # SA
# subject_id = 'fsaverage';
# colmap = colorRampPalette(rev(brewer.pal(11, name="RdBu")))
# lh_demo_cluster_file = file.path(plsdir,'lh_SA_LV1_BSR_sexAdj.mgh')
# rh_demo_cluster_file = file.path(plsdir,'rh_SA_LV1_BSR_sexAdj.mgh')
# 
# lh_SA = read.fs.morph(lh_demo_cluster_file);   # contains a single positive cluster
# rh_SA = read.fs.morph(rh_demo_cluster_file);  
# lh_SA[lh_SA==0] = NA
# rh_SA[rh_SA==0] = NA
# 
# msh = vis.data.on.fsaverage(subjects_dir = subjects_dir, 
#                             vis_subject_id = 'fsaverage', 
#                             morph_data_lh = lh_SA, 
#                             morph_data_rh = rh_SA,
#                             surface = "pial",
#                             makecmap_options = list('colFn'=colmap,'range'=c(-35, -5)))
# img = vis.export.from.coloredmeshes(msh, 
#                                     colorbar_legend = "SA", 
#                                     output_img = file.path(plsdir,'SA_LV1_BSR_allVertices_sexAdj.png'))
# #===============================================================================
# 
# 
# #===============================================================================
# # CT
# subject_id = 'fsaverage';
# colmap = colorRampPalette(rev(brewer.pal(11, name="RdBu")))
# lh_demo_cluster_file = file.path(plsdir,'lh_CT_LV1_BSR_sexAdj.mgh')
# rh_demo_cluster_file = file.path(plsdir,'rh_CT_LV1_BSR_sexAdj.mgh')
# 
# lh_CT = read.fs.morph(lh_demo_cluster_file);   # contains a single positive cluster
# rh_CT = read.fs.morph(rh_demo_cluster_file);  
# lh_CT[lh_CT==0] = NA
# rh_CT[rh_CT==0] = NA
# 
# msh = vis.data.on.fsaverage(subjects_dir = subjects_dir, 
#                             vis_subject_id = 'fsaverage', 
#                             morph_data_lh = lh_CT, 
#                             morph_data_rh = rh_CT,
#                             surface = "pial",
#                             makecmap_options = list('colFn'=colmap,'range'=c(-10, 25)))
# img = vis.export.from.coloredmeshes(msh, 
#                                     colorbar_legend = "CT", 
#                                     output_img = file.path(plsdir,'CT_LV1_BSR_allVertices_sexAdj.png'))
# #===============================================================================
# 
# 
# 
# 
# 
# #===============================================================================
# # SA
# subject_id = 'fsaverage';
# colmap = colorRampPalette(rev(brewer.pal(11, name="RdBu")))
# lh_demo_cluster_file = file.path(plsdir,'lh_SA_log10SAsexAdj_LV1_BSR.mgh')
# rh_demo_cluster_file = file.path(plsdir,'rh_SA_log10SAsexAdj_LV1_BSR.mgh')
# 
# lh_SA = read.fs.morph(lh_demo_cluster_file);   # contains a single positive cluster
# rh_SA = read.fs.morph(rh_demo_cluster_file);  
# lh_SA[lh_SA==0] = NA
# rh_SA[rh_SA==0] = NA
# 
# msh = vis.data.on.fsaverage(subjects_dir = subjects_dir, 
#                             vis_subject_id = 'fsaverage', 
#                             morph_data_lh = lh_SA, 
#                             morph_data_rh = rh_SA,
#                             surface = "pial",
#                             makecmap_options = list('colFn'=colmap,'range'=c(-35, -10)))
# img = vis.export.from.coloredmeshes(msh, 
#                                     colorbar_legend = "SA", 
#                                     output_img = file.path(plsdir,'SA_LV1_BSR_allVertices_sexlog10SAAdj.png'))
# #===============================================================================



#===============================================================================
# SA LV1
subject_id = 'fsaverage'
colmap = colorRampPalette(c("blue","white","red"))
# colmap = colorRampPalette(rev(brewer.pal(11, name="RdBu")))
lh_demo_cluster_file = file.path(plsdir,'lh_SA_LV1_BSR_sexMeanSAAdj.mgh')
rh_demo_cluster_file = file.path(plsdir,'rh_SA_LV1_BSR_sexMeanSAAdj.mgh')

lh_SA = read.fs.morph(lh_demo_cluster_file);   # contains a single positive cluster
rh_SA = read.fs.morph(rh_demo_cluster_file);  
lh_SA[lh_SA==0] = NA
rh_SA[rh_SA==0] = NA

msh = vis.data.on.fsaverage(subjects_dir = subjects_dir, 
                            vis_subject_id = 'fsaverage', 
                            morph_data_lh = lh_SA, 
                            morph_data_rh = rh_SA,
                            surface = "pial",
                            makecmap_options = list('colFn'=colmap,'range'=c(-10, 10)))
img = vis.export.from.coloredmeshes(msh, 
                                    colorbar_legend = "SA", 
                                    output_img = file.path(plsdir,'SA_LV1_BSR_allVertices_sexMeanSAAdj.png'))
#===============================================================================



#===============================================================================
# CT LV1
subject_id = 'fsaverage'
colmap = colorRampPalette(c("blue","white","red"))
# colmap = colorRampPalette(rev(brewer.pal(11, name="RdBu")))
lh_demo_cluster_file = file.path(plsdir,'lh_CT_LV1_BSR_sexMeanCTAdj.mgh')
rh_demo_cluster_file = file.path(plsdir,'rh_CT_LV1_BSR_sexMeanCTAdj.mgh')

lh_CT1 = read.fs.morph(lh_demo_cluster_file);   # contains a single positive cluster
rh_CT1 = read.fs.morph(rh_demo_cluster_file);  
lh_CT1[lh_CT1==0] = NA
rh_CT1[rh_CT1==0] = NA

msh = vis.data.on.fsaverage(subjects_dir = subjects_dir, 
                            vis_subject_id = 'fsaverage', 
                            morph_data_lh = lh_CT1, 
                            morph_data_rh = rh_CT1,
                            surface = "pial",
                            makecmap_options = list('colFn'=colmap,'range'=c(-10, 10)))
img = vis.export.from.coloredmeshes(msh, 
                                    colorbar_legend = "CT", 
                                    output_img = file.path(plsdir,'CT_LV1_BSR_allVertices_sexMeanCTAdj.png'))
#===============================================================================

#===============================================================================
# CT LV2
subject_id = 'fsaverage'
colmap = colorRampPalette(c("blue","white","red"))
# colmap = colorRampPalette(rev(brewer.pal(11, name="RdBu")))
lh_demo_cluster_file = file.path(plsdir,'lh_CT_LV2_BSR_sexMeanCTAdj.mgh')
rh_demo_cluster_file = file.path(plsdir,'rh_CT_LV2_BSR_sexMeanCTAdj.mgh')

lh_CT2 = read.fs.morph(lh_demo_cluster_file);   # contains a single positive cluster
rh_CT2 = read.fs.morph(rh_demo_cluster_file);  
lh_CT2[lh_CT2==0] = NA
rh_CT2[rh_CT2==0] = NA

msh = vis.data.on.fsaverage(subjects_dir = subjects_dir, 
                            vis_subject_id = 'fsaverage', 
                            morph_data_lh = lh_CT2, 
                            morph_data_rh = rh_CT2,
                            surface = "pial",
                            makecmap_options = list('colFn'=colmap,'range'=c(-10, 10)))
img = vis.export.from.coloredmeshes(msh, 
                                    colorbar_legend = "CT", 
                                    output_img = file.path(plsdir,'CT_LV2_BSR_allVertices_sexMeanCTAdj.png'))
#===============================================================================





#===============================================================================
# scaling factor computed on same dataset
subject_id = 'fsaverage'
colmap = colorRampPalette(rev(c("blue","white","red")))
# colmap = colorRampPalette(brewer.pal(11, name="RdBu"))
lh_demo_cluster_file = here("data","tidy",'lh_scalingFactorBeta.mgh')
rh_demo_cluster_file = here("data","tidy",'rh_scalingFactorBeta.mgh')

lh_sf = read.fs.morph(lh_demo_cluster_file);   # contains a single positive cluster
rh_sf = read.fs.morph(rh_demo_cluster_file);  
lh_sf[lh_sf==0] = NA
rh_sf[rh_sf==0] = NA
# lh_sf = lh_sf*-1
# rh_sf = rh_sf*-1

msh = vis.data.on.fsaverage(subjects_dir = subjects_dir, 
                            vis_subject_id = 'fsaverage5', 
                            morph_data_lh = lh_sf, 
                            morph_data_rh = rh_sf,
                            surface = "pial",
                            makecmap_options = list('colFn'=colmap,'range'=c(0.5, 1.3)))
img = vis.export.from.coloredmeshes(msh, 
                                    colorbar_legend = "SA", 
                                    output_img = file.path(plsdir,'scalingFactorBeta.png'))
#===============================================================================




# #===============================================================================
# # scaling factor from Reardon et al.,
# subject_id = 'fsaverage5'
# colmap = colorRampPalette(c("blue","white","red"))
# # colmap = colorRampPalette(brewer.pal(11, name="RdBu"))
# lh_demo_cluster_file = file.path(codedir,"Brain_Organization","AllometricScaling","lh.AllometricScaling_fsaverage5.func.gii")
# rh_demo_cluster_file = file.path(codedir,"Brain_Organization","AllometricScaling","rh.AllometricScaling_fsaverage5.func.gii")
# 
# lh_sf_reardon = read.fs.morph(lh_demo_cluster_file, format = "gii");   # contains a single positive cluster
# rh_sf_reardon = read.fs.morph(rh_demo_cluster_file, format = "gii");  
# lh_sf_reardon[lh_sf_reardon==0] = NA
# rh_sf_reardon[rh_sf_reardon==0] = NA
# # lh_sf_reardon = lh_sf_reardon*-1
# # rh_sf_reardon = rh_sf_reardon*-1
# 
# msh = vis.data.on.fsaverage(subjects_dir = subjects_dir, 
#                             vis_subject_id = 'fsaverage5', 
#                             morph_data_lh = lh_sf_reardon, 
#                             morph_data_rh = rh_sf_reardon,
#                             surface = "pial",
#                             makecmap_options = list('colFn'=colmap,'range'=c(0.5, 1.5)))
# img = vis.export.from.coloredmeshes(msh, 
#                                     colorbar_legend = "SA", 
#                                     output_img = file.path(plsdir,'scalingFactorBeta_Reardon.png'))
# #===============================================================================


#===============================================================================
# Principal Gradient,
subject_id = 'fsaverage5'
colmap = colorRampPalette(rev(c("blue","white","red")))
# colmap = colorRampPalette(brewer.pal(11, name="RdBu"))
lh_demo_cluster_file = file.path(codedir,"Brain_Organization","PrincipleGradient","Gradients.lh.fsaverage5.func.gii")
rh_demo_cluster_file = file.path(codedir,"Brain_Organization","PrincipleGradient","Gradients.rh.fsaverage5.func.gii")

lh_pcg = read.fs.morph(lh_demo_cluster_file, format = "gii");   # contains a single positive cluster
rh_pcg = read.fs.morph(rh_demo_cluster_file, format = "gii");  
lh_pcg[lh_pcg==0] = NA
rh_pcg[rh_pcg==0] = NA
# lh_sf_reardon = lh_sf_reardon*-1
# rh_sf_reardon = rh_sf_reardon*-1

msh = vis.data.on.fsaverage(subjects_dir = subjects_dir, 
                            vis_subject_id = 'fsaverage5', 
                            morph_data_lh = lh_pcg, 
                            morph_data_rh = rh_pcg,
                            surface = "pial",
                            makecmap_options = list('colFn'=colmap,'range'=c(-10, 10)))
img = vis.export.from.coloredmeshes(msh, 
                                    colorbar_legend = "SA", 
                                    output_img = file.path(plsdir,'principalGradient.png'))
#===============================================================================


# #===============================================================================
# # Evolutionary Expansion,
# subject_id = 'fsaverage5'
# colmap = colorRampPalette(brewer.pal(11, name="RdBu"))
# # lh_demo_cluster_file = file.path(codedir,"Brain_Organization","EvolutionaryExpansion","Gradients.lh.fsaverage5.func.gii")
# rh_demo_cluster_file = file.path(codedir,"Brain_Organization","EvolutionaryExpansion","rh.Hill2010_evo_fsaverage5.func.gii")
# 
# # lh_pcg = read.fs.morph(lh_demo_cluster_file, format = "gii");   # contains a single positive cluster
# rh_evo = read.fs.morph(rh_demo_cluster_file, format = "gii");  
# # lh_pcg[lh_pcg==0] = NA
# rh_evo[rh_evo==0] = NA
# # lh_sf_reardon = lh_sf_reardon*-1
# # rh_sf_reardon = rh_sf_reardon*-1
# 
# msh = vis.data.on.fsaverage(subjects_dir = subjects_dir, 
#                             vis_subject_id = 'fsaverage5', 
#                             morph_data_lh = NULL, 
#                             morph_data_rh = rh_evo,
#                             surface = "pial",
#                             makecmap_options = list('colFn'=colmap)) #,'range'=c(-12, 12)))
# img = vis.export.from.coloredmeshes(msh, 
#                                     colorbar_legend = "SA", 
#                                     output_img = file.path(plsdir,'evoluationaryExpansion.png'))
# #===============================================================================


#===============================================================================
# Myelin Map,
subject_id = 'fsaverage5'
colmap = colorRampPalette(c("blue","white","red"))
# colmap = colorRampPalette(rev(brewer.pal(11, name="RdBu")))
lh_demo_cluster_file = file.path(codedir,"Brain_Organization","Myelin","MyelinMap.lh.fsaverage5.func.gii")
rh_demo_cluster_file = file.path(codedir,"Brain_Organization","Myelin","MyelinMap.rh.fsaverage5.func.gii")

lh_my = read.fs.morph(lh_demo_cluster_file, format = "gii");   # contains a single positive cluster
rh_my = read.fs.morph(rh_demo_cluster_file, format = "gii");  
lh_my[lh_my==0] = NA
rh_my[rh_my==0] = NA
# lh_sf_reardon = lh_sf_reardon*-1
# rh_sf_reardon = rh_sf_reardon*-1

msh = vis.data.on.fsaverage(subjects_dir = subjects_dir, 
                            vis_subject_id = 'fsaverage5', 
                            morph_data_lh = lh_my, 
                            morph_data_rh = rh_my,
                            surface = "pial",
                            makecmap_options = list('colFn'=colmap,'range'=c(1, 1.55)))
img = vis.export.from.coloredmeshes(msh, 
                                    colorbar_legend = "SA", 
                                    output_img = file.path(plsdir,'myelinMap.png'))
#===============================================================================


# #===============================================================================
# # Mean CBF
# subject_id = 'fsaverage5'
# colmap = colorRampPalette(rev(brewer.pal(11, name="RdBu")))
# lh_demo_cluster_file = file.path(codedir,"Brain_Organization","MeanCBF","lh.MeanCBF.fsaverage5.func.gii")
# rh_demo_cluster_file = file.path(codedir,"Brain_Organization","MeanCBF","rh.MeanCBF.fsaverage5.func.gii")
# 
# lh_cbf = read.fs.morph(lh_demo_cluster_file, format = "gii");   # contains a single positive cluster
# rh_cbf = read.fs.morph(rh_demo_cluster_file, format = "gii");  
# lh_cbf[lh_cbf==0] = NA
# rh_cbf[rh_cbf==0] = NA
# # lh_sf_reardon = lh_sf_reardon*-1
# # rh_sf_reardon = rh_sf_reardon*-1
# 
# msh = vis.data.on.fsaverage(subjects_dir = subjects_dir, 
#                             vis_subject_id = 'fsaverage5', 
#                             morph_data_lh = lh_cbf, 
#                             morph_data_rh = rh_cbf,
#                             surface = "pial",
#                             makecmap_options = list('colFn'=colmap,'range'=c(40, 80)))
# img = vis.export.from.coloredmeshes(msh, 
#                                     colorbar_legend = "SA", 
#                                     output_img = file.path(plsdir,'meanCBF.png'))
# #===============================================================================
# 


# make correlation matrix
map_names = c("SA LV1", "CT LV1", "CT LV2", "Scaling", "Functional","Myelin")
r_mat = matrix(nrow = length(map_names), ncol = length(map_names))

r_mat[1,1] = NA
r_mat[1,2] = cor(c(lh_SA[1:length(lh_CT1)], rh_SA[1:length(rh_CT1)]), c(lh_CT1, rh_CT1), use = "pairwise.complete.obs")
r_mat[2,1] = r_mat[1,2]
r_mat[1,3] = cor(c(lh_SA[1:length(lh_CT2)], rh_SA[1:length(rh_CT2)]), c(lh_CT2, rh_CT2), use = "pairwise.complete.obs")
r_mat[3,1] = r_mat[1,3]
r_mat[1,4] = cor(c(lh_SA[1:length(lh_sf)], rh_SA[1:length(rh_sf)]), c(lh_sf, rh_sf), use = "pairwise.complete.obs")
r_mat[4,1] = r_mat[1,4]
r_mat[1,5] = cor(c(lh_SA[1:length(lh_pcg)], rh_SA[1:length(rh_pcg)]), c(lh_pcg, rh_pcg), use = "pairwise.complete.obs")
r_mat[5,1] = r_mat[1,5]
r_mat[1,6] = cor(c(lh_SA[1:length(lh_my)], rh_SA[1:length(rh_my)]), c(lh_my, rh_my), use = "pairwise.complete.obs")
r_mat[6,1] = r_mat[1,6]

r_mat[2,3] = cor(c(lh_CT1[1:length(lh_CT2)], rh_CT1[1:length(rh_CT2)]), c(lh_CT2, rh_CT2), use = "pairwise.complete.obs")
r_mat[3,2] = r_mat[2,3]
r_mat[2,4] = cor(c(lh_CT1[1:length(lh_sf)], rh_CT1[1:length(rh_sf)]), c(lh_sf, rh_sf), use = "pairwise.complete.obs")
r_mat[4,2] = r_mat[2,4]
r_mat[2,5] = cor(c(lh_CT1[1:length(lh_pcg)], rh_CT1[1:length(rh_pcg)]), c(lh_pcg, rh_pcg), use = "pairwise.complete.obs")
r_mat[5,2] = r_mat[2,5]
r_mat[2,6] = cor(c(lh_CT1[1:length(lh_my)], rh_CT1[1:length(rh_my)]), c(lh_my, rh_my), use = "pairwise.complete.obs")
r_mat[6,2] = r_mat[2,6]

r_mat[3,4] = cor(c(lh_CT2[1:length(lh_sf)], rh_CT2[1:length(rh_sf)]), c(lh_sf, rh_sf), use = "pairwise.complete.obs")
r_mat[4,3] = r_mat[3,4]
r_mat[3,5] = cor(c(lh_CT2[1:length(lh_pcg)], rh_CT2[1:length(rh_pcg)]), c(lh_pcg, rh_pcg), use = "pairwise.complete.obs")
r_mat[5,3] = r_mat[3,5]
r_mat[3,6] = cor(c(lh_CT2[1:length(lh_my)], rh_CT2[1:length(rh_my)]), c(lh_my, rh_my), use = "pairwise.complete.obs")
r_mat[6,3] = r_mat[3,6]

r_mat[4,5] = cor(c(lh_sf[1:length(lh_pcg)], rh_sf[1:length(rh_pcg)]), c(lh_pcg, rh_pcg), use = "pairwise.complete.obs")
r_mat[5,4] = r_mat[4,5]
r_mat[4,6] = cor(c(lh_sf[1:length(lh_my)], rh_sf[1:length(rh_my)]), c(lh_my, rh_my), use = "pairwise.complete.obs")
r_mat[6,4] = r_mat[4,6]

r_mat[5,6] = cor(c(lh_pcg[1:length(lh_my)], rh_pcg[1:length(rh_my)]), c(lh_my, rh_my), use = "pairwise.complete.obs")
r_mat[6,5] = r_mat[5,6]

r_mat
