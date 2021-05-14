

library(here)
library(ggplot2)
source(here("code","cohens_d.R"))

plspath = here("pls")
resultpath = here("pls","results")

fontSize = 20
dotSize = 10

#------------------------------------------------------------------------------
# SA LV1
data2plot = read.csv(file.path(resultpath,"AllVertices_sexMeanSAAdj_SA_LV1_withGCLUST4plotting.csv"))
data2plot$bsr[data2plot$bsr==0] = NA

p = ggplot(data = data2plot, aes(x = cluster, y = bsr, fill = cluster, group=cluster)) + 
  geom_violin() + geom_boxplot(fill=NA) + scale_fill_manual(values = c("blue","red")) +  
  xlab("Clusters") + ylab("Bootstrap Ratio (BSR)") + ggtitle("Anterior-Posterior Clusters SA PLS LV1") +
  guides(fill=FALSE) + 
  theme(text = element_text(size=fontSize), 
        axis.text.x = element_text(size=fontSize),
        axis.text.y = element_text(size=fontSize),
        plot.title = element_text(size=fontSize,hjust=0.5))
ggsave(filename = here("plots","AllVertices_sexMeanSAAdj_SA_LV1_BSRviolinPlot_byAPclusters.pdf"),width=6,height=6)
p

# t-stats
t_res = t.test(formula = bsr ~ cluster, data = data2plot)
t_res 
# effect size
d_res = cohens_d(data2plot$bsr[data2plot$cluster=="anterior"], data2plot$bsr[data2plot$cluster=="posterior"])
d_res

p = ggplot(data = data2plot, aes(x = grad_order, y = bsr, fill = grad_order, group=grad_order)) + 
  geom_violin() + geom_boxplot(fill=NA) + xlim(0,12) +
  scale_x_continuous(breaks=seq(1,12,1)) + 
  scale_fill_gradientn(colors = colorRampPalette(c("blue","white","red"))(100)) +
  ylab("Bootstrap Ratio (BSR)") + xlab("Brain Region") +  ggtitle("Genetic Similarity Gradient SA PLS LV1") +
  guides(fill=FALSE) + 
  theme(text = element_text(size=fontSize), 
        axis.text.x = element_text(size=fontSize),
        axis.text.y = element_text(size=fontSize),
        plot.title = element_text(size=fontSize,hjust=0.5))
ggsave(filename = here("plots","AllVertices_sexMeanSAAdj_SA_LV1_BSRviolinplot_byGeneSimGradient.pdf"),width=6,height=6)
p


lh_df4plot = data.frame(matrix(nrow = 12, ncol = 2))
colnames(lh_df4plot) = c("grad_order", "medBSR","hemi")
for (ireg in 1:12){
  lh_df4plot[ireg,"hemi"] = "left"
  lh_df4plot[ireg,"grad_order"] = ireg
  lh_df4plot[ireg,"medBSR"] = median(data2plot$bsr[data2plot$grad_order==ireg & data2plot$hemi=="left"], na.rm = TRUE)
}

rh_df4plot = data.frame(matrix(nrow = 12, ncol = 2))
colnames(rh_df4plot) = c("grad_order", "medBSR","hemi")
for (ireg in 1:12){
  rh_df4plot[ireg,"hemi"] = "right"
  rh_df4plot[ireg,"grad_order"] = ireg
  rh_df4plot[ireg,"medBSR"] = median(data2plot$bsr[data2plot$grad_order==ireg & data2plot$hemi=="right"], na.rm = TRUE)
}
df4plot = rbind(lh_df4plot, rh_df4plot)

p = ggplot(data = df4plot, aes(x = grad_order, y = medBSR, fill=grad_order))
p = p + geom_point(size=dotSize, colour="black", pch=21) + 
  geom_smooth(method=lm, colour="black") + xlim(0,12) +
  scale_x_continuous(breaks=seq(1,12,1)) + 
  scale_fill_gradientn(colors = colorRampPalette(c("blue","white","red"))(100)) +
  guides(colour=FALSE, fill=FALSE) +  
  xlab("Brain Region") + 
  ylab("Bootstrap Ratio (BSR)") + 
  ggtitle("Genetic Similarity Gradient SA PLS LV1") +
  theme(text = element_text(size=fontSize), 
        axis.text.x = element_text(size=fontSize),
        axis.text.y = element_text(size=fontSize),
        plot.title = element_text(size=fontSize,hjust=0.5))
ggsave(filename = here("plots","AllVertices_sexMeanSAAdj_SA_LV1_medianBSRscatterplot_byGeneSimGradient.pdf"),width=6,height=6)
p

cor.test(df4plot$grad_order,df4plot$medBSR)




#------------------------------------------------------------------------------
# CT LV1
data2plot = read.csv(file.path(resultpath,"AllVertices_sexMeanCTAdj_CT_LV1_withGCLUST4plotting.csv"))
data2plot$bsr[data2plot$bsr==0] = NA

p = ggplot(data = data2plot, aes(x = cluster, y = bsr, fill = cluster, group=cluster)) + 
  geom_violin() + geom_boxplot(fill=NA) + scale_fill_manual(values = c("red","blue")) +  
  xlab("Clusters") + ylab("Bootstrap Ratio (BSR)") + ggtitle("Dorsal-Ventral Clusters CT PLS LV1") +
  guides(fill=FALSE) + 
  theme(text = element_text(size=fontSize), 
        axis.text.x = element_text(size=fontSize),
        axis.text.y = element_text(size=fontSize),
        plot.title = element_text(size=fontSize,hjust=0.5))
ggsave(filename = here("plots","AllVertices_sexMeanCTAdj_CT_LV1_BSRviolinPlot_byDVclusters.pdf"),width=6,height=6)
p

# t-stats
t_res = t.test(formula = bsr ~ cluster, data = data2plot)
t_res 
# effect size
d_res = cohens_d(data2plot$bsr[data2plot$cluster=="dorsal"], data2plot$bsr[data2plot$cluster=="ventral"])
d_res


p = ggplot(data = data2plot, aes(x = grad_order, y = bsr, fill = grad_order, group=grad_order)) + 
  geom_violin() + geom_boxplot(fill=NA) + xlim(0,12) +
  scale_x_continuous(breaks=seq(1,12,1)) + 
  scale_fill_gradientn(colors = colorRampPalette(c("red","white","blue"))(100)) +
  ylab("Bootstrap Ratio (BSR)") + xlab("Brain Region") +  ggtitle("Genetic Similarity Gradient CT PLS LV1") +
  guides(fill=FALSE) + 
  theme(text = element_text(size=fontSize), 
        axis.text.x = element_text(size=fontSize),
        axis.text.y = element_text(size=fontSize),
        plot.title = element_text(size=fontSize,hjust=0.5))
ggsave(filename = here("plots","AllVertices_sexMeanCTAdj_CT_LV1_BSRviolinplot_byGeneSimGradient.pdf"),width=6,height=6)
p


lh_df4plot = data.frame(matrix(nrow = 12, ncol = 2))
colnames(lh_df4plot) = c("grad_order", "medBSR","hemi")
for (ireg in 1:12){
  lh_df4plot[ireg,"hemi"] = "left"
  lh_df4plot[ireg,"grad_order"] = ireg
  lh_df4plot[ireg,"medBSR"] = median(data2plot$bsr[data2plot$grad_order==ireg & data2plot$hemi=="left"], na.rm = TRUE)
}

rh_df4plot = data.frame(matrix(nrow = 12, ncol = 2))
colnames(rh_df4plot) = c("grad_order", "medBSR","hemi")
for (ireg in 1:12){
  rh_df4plot[ireg,"hemi"] = "right"
  rh_df4plot[ireg,"grad_order"] = ireg
  rh_df4plot[ireg,"medBSR"] = median(data2plot$bsr[data2plot$grad_order==ireg & data2plot$hemi=="right"], na.rm = TRUE)
}
df4plot = rbind(lh_df4plot, rh_df4plot)

p = ggplot(data = df4plot, aes(x = grad_order, y = medBSR, fill=grad_order))
p = p + geom_point(size=dotSize, colour="black", pch=21) + 
  geom_smooth(method=lm, colour="black") + xlim(0,12) +
  scale_x_continuous(breaks=seq(1,12,1)) + 
  scale_fill_gradientn(colors = colorRampPalette(c("red","white","blue"))(100)) +
  guides(colour=FALSE, fill=FALSE) +  
  xlab("Brain Region") + 
  ylab("Bootstrap Ratio (BSR)") + 
  ggtitle("Genetic Similarity Gradient CT PLS LV1") +
  theme(text = element_text(size=fontSize), 
        axis.text.x = element_text(size=fontSize),
        axis.text.y = element_text(size=fontSize),
        plot.title = element_text(size=fontSize,hjust=0.5))
ggsave(filename = here("plots","AllVertices_sexMeanCTAdj_CT_LV1_medianBSRscatterplot_byGeneSimGradient.pdf"),width=6,height=6)
p

cor.test(df4plot$grad_order,df4plot$medBSR)


#------------------------------------------------------------------------------
# CT LV2
data2plot = read.csv(file.path(resultpath,"AllVertices_sexMeanCTAdj_CT_LV2_withGCLUST4plotting.csv"))
data2plot$bsr[data2plot$bsr==0] = NA

p = ggplot(data = data2plot, aes(x = cluster, y = bsr, fill = cluster, group=cluster)) + 
  geom_violin() + geom_boxplot(fill=NA) + scale_fill_manual(values = c("red","blue")) +  
  xlab("Clusters") + ylab("Bootstrap Ratio (BSR)") + ggtitle("Dorsal-Ventral Clusters CT PLS LV2") +
  guides(fill=FALSE) + 
  theme(text = element_text(size=fontSize), 
        axis.text.x = element_text(size=fontSize),
        axis.text.y = element_text(size=fontSize),
        plot.title = element_text(size=fontSize,hjust=0.5))
ggsave(filename = here("plots","AllVertices_sexMeanCTAdj_CT_LV2_BSRviolinPlot_byDVclusters.pdf"),width=6,height=6)
p

# t-stats
t_res = t.test(formula = bsr ~ cluster, data = data2plot)
t_res 
# effect size
d_res = cohens_d(data2plot$bsr[data2plot$cluster=="dorsal"], data2plot$bsr[data2plot$cluster=="ventral"])
d_res


p = ggplot(data = data2plot, aes(x = grad_order, y = bsr, fill = grad_order, group=grad_order)) + 
  geom_violin() + geom_boxplot(fill=NA) + xlim(0,12) +
  scale_x_continuous(breaks=seq(1,12,1)) + 
  scale_fill_gradientn(colors = colorRampPalette(c("red","white","blue"))(100)) +
  ylab("Bootstrap Ratio (BSR)") + xlab("Brain Region") +  ggtitle("Genetic Similarity Gradient CT PLS LV2") +
  guides(fill=FALSE) + 
  theme(text = element_text(size=fontSize), 
        axis.text.x = element_text(size=fontSize),
        axis.text.y = element_text(size=fontSize),
        plot.title = element_text(size=fontSize,hjust=0.5))
ggsave(filename = here("plots","AllVertices_sexMeanCTAdj_CT_LV2_BSRviolinplot_byGeneSimGradient.pdf"),width=6,height=6)
p


lh_df4plot = data.frame(matrix(nrow = 12, ncol = 2))
colnames(lh_df4plot) = c("grad_order", "medBSR","hemi")
for (ireg in 1:12){
  lh_df4plot[ireg,"hemi"] = "left"
  lh_df4plot[ireg,"grad_order"] = ireg
  lh_df4plot[ireg,"medBSR"] = median(data2plot$bsr[data2plot$grad_order==ireg & data2plot$hemi=="left"], na.rm = TRUE)
}

rh_df4plot = data.frame(matrix(nrow = 12, ncol = 2))
colnames(rh_df4plot) = c("grad_order", "medBSR","hemi")
for (ireg in 1:12){
  rh_df4plot[ireg,"hemi"] = "right"
  rh_df4plot[ireg,"grad_order"] = ireg
  rh_df4plot[ireg,"medBSR"] = median(data2plot$bsr[data2plot$grad_order==ireg & data2plot$hemi=="right"], na.rm = TRUE)
}
df4plot = rbind(lh_df4plot, rh_df4plot)

p = ggplot(data = df4plot, aes(x = grad_order, y = medBSR, fill=grad_order))
p = p + geom_point(size=dotSize, colour="black", pch=21) + 
  geom_smooth(method=lm, colour="black") + xlim(0,12) +
  scale_x_continuous(breaks=seq(1,12,1)) + 
  scale_fill_gradientn(colors = colorRampPalette(c("red","white","blue"))(100)) +
  guides(colour=FALSE, fill=FALSE) +  
  xlab("Brain Region") + 
  ylab("Bootstrap Ratio (BSR)") + 
  ggtitle("Genetic Similarity Gradient CT PLS LV2") +
  theme(text = element_text(size=fontSize), 
        axis.text.x = element_text(size=fontSize),
        axis.text.y = element_text(size=fontSize),
        plot.title = element_text(size=fontSize,hjust=0.5))
ggsave(filename = here("plots","AllVertices_sexMeanCTAdj_CT_LV2_medianBSRscatterplot_byGeneSimGradient.pdf"),width=6,height=6)
p

cor.test(df4plot$grad_order,df4plot$medBSR)


#------------------------------------------------------------------------------
# plot crossblock covariance explained for GCLUST versus All Vertices
fontSize = 14

cov_expl = data.frame(atype = c("GCLUST","Vertices"),SA_LV1 = c(0.3647,0.1780),
           CT_LV1 = c(0.3779,0.1786), CT_LV2 = c(0.1963,0.1125))
cov_expl = melt(cov_expl)

p = ggplot(data = cov_expl, aes(x = atype, y = value, fill = atype)) + facet_grid(. ~ variable) + 
  geom_bar(stat = "identity") + xlab("Model") + ylab("Percentage Covariance Explained") + guides(fill=FALSE) + 
  ggtitle("Covariance Explained GCLUST vs All Vertices") +
  theme(text = element_text(size=fontSize), 
        axis.text.x = element_text(size=fontSize),
        axis.text.y = element_text(size=fontSize),
        plot.title = element_text(size=fontSize,hjust=0.5))
ggsave(filename = here("plots","PLS_CovarianceExplained_GCLUST_AllVertices.pdf"),width=6,height=6)
p 
