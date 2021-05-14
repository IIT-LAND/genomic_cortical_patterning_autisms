# load libraries
library(easypackages)
libraries("here","limma","readxl")
options(stringsAsFactors = FALSE)

RUNONSERVER = TRUE
if(RUNONSERVER){
  rootpath = "/media/DATA/mlombardo/ace_smri_gex_pls"
} else {
  rootpath = here()
}

# load in MRI data
print("Loading data")
mri_data = read.csv(file.path(rootpath,"data","tidy","AllGroups_all_cortical_ct_allVertices.csv"))
rownames(mri_data) = mri_data[,"subjectId"]
features = colnames(mri_data)[4:ncol(mri_data)]

mri_data$globalMean = rowMeans(mri_data[,features], na.rm=TRUE)

# load in pheno data for MRI
pheno_data = read.csv(file.path(rootpath,"data","tidy","labelData_all_MRI.csv"))
df2use = merge(mri_data, pheno_data, by="subjectId", sort = FALSE)

mri_data_adj = data.frame(matrix(nrow = dim(mri_data)[1], ncol = length(features)+1))
rownames(mri_data_adj) = rownames(mri_data)
colnames(mri_data_adj) = c("subjectId",features)
mri_data_adj$subjectId = rownames(mri_data)

for (feature in features){

  print(feature)

  # formula
  form2use = sprintf("%s ~ subgrp + sex + globalMean",feature)

  # model
  mod2use = lm(formula = form2use, data = df2use)

  # remove variation from covariate
  covname2use = c("sexMale","globalMean")
  beta1 = mod2use$coefficients[covname2use, drop = FALSE]
  beta1[is.na(beta1)] = 0

  # construct model
  cov_columns = c("subgrp","sex","globalMean")
  full_model = model.matrix(~0+as.factor(subgrp) + as.factor(sex) + globalMean, data = df2use)
  colnames(full_model) = c("Good","Poor","TD","sex","globalMean")

  # remove sex covariate from mri_data
  covname2use = c("sex","globalMean")
  mri_data_adj[,feature] = t(mri_data[,feature] - beta1 %*% t(full_model[,covname2use]))
}

# write adjusted data out to a file
write.csv(mri_data_adj, file = file.path(rootpath,"data","tidy","mri_data_adj_sexMeanCT_all_cortical_ct_allVertices.csv"))
