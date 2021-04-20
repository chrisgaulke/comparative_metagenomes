library(vegan)
library(ape)
library(phangorn)
library(phytools)
library(ggplot2)
library(qvalue)
library(gplots)
library(RColorBrewer)



# IMPORT PHYLO DATA: Mammalian Tree ------------------------------------------

#This tree is a large tree of many mammals so I need to select the ones I care about.
mammal.tre <- read.nexus("/Users/cgaulke/unsynced_projects/comp_metagenomes/data/mammal_tree/mammalianSuperTreeRevisedVersion2.txt")
mammal.keeps <- c("Macaca_mulatta",
                  "Macaca_fuscata",
                  "Macaca_fascicularis",
                  "Homo_sapiens",
                  "Bos_taurus",
                  "Sus_scrofa",
                  "Canis_lupus",
                  "Mus_musculus",
                  "Rattus_norvegicus"
)
mammal.drops <- mammal.tre$tip.label[-which(mammal.tre$tip.label %in% mammal.keeps)]
mammal.tre <- drop.tip(mammal.tre,mammal.drops )
plot(mammal.tre)
mammal.dist <- cophenetic(mammal.tre)



# IMPORT MICROBIOME DATA: KO -------------------------------------------------

#prep microbiome

ko_gm.df <- read.table("/Users/cgaulke/unsynced_projects/comp_metagenomes/analysis/2019_07_10_flat_files/2019_07_10_ko_gm.txt",
           sep ="\t")

#Don't want env samples or fish for this analysis
phylo_ko.gm <- ko_gm.df[c(1:8,10),]

#change to match names and order to the tree
rownames(phylo_ko.gm)[5] <- "Canis_lupus"

phylo_ko.gm <-
  phylo_ko.gm[mammal.tre$tip.label,]

phylo_ko.gm <- phylo_ko.gm[,which(colSums(phylo_ko.gm) > 0),drop=F]

#make a distnace matrix
phylo_ko.dist <- vegdist(phylo_ko.gm,
                          method = "bray")

phylo_ko_dist.mat <-
  as.matrix(phylo_ko.dist)

#now sort to match the mammal tree

phylo_ko_dist.mat <-
  phylo_ko_dist.mat[rownames(mammal.dist),
  colnames(mammal.dist),
 drop = F]


# IMPORT MICROBIOME DATA: Module ---------------------------------------------

module_gm.df <- read.table("/Users/cgaulke/unsynced_projects/comp_metagenomes/analysis/2019_07_10_flat_files/2019_07_10_module_gm.txt",
                       sep ="\t")

#Don't want env samples or fish for this analysis
phylo_module.gm <- module_gm.df[c(1:8,10),]

#change to match names and order to the tree
rownames(phylo_module.gm)[5] <- "Canis_lupus"

phylo_module.gm <-
  phylo_module.gm[mammal.tre$tip.label,]

#make a distnace matrix
phylo_module.dist <- vegdist(phylo_module.gm,
                         method = "bray")

phylo_module_dist.mat <-
  as.matrix(phylo_module.dist)

#now sort to match the mammal tree

phylo_module_dist.mat <-
  phylo_module_dist.mat[rownames(mammal.dist),
                    colnames(mammal.dist),
                    drop = F]

# FUNCTIONAL PHYLOSYMBIOSIS: KO Mantel ---------------------------------------

#Calculate the distances
dist.df <- data.frame(micro=phylo_ko_dist.mat[upper.tri(phylo_ko_dist.mat)],
                         tree=mammal.dist[upper.tri(mammal.dist)])

dist.plot <- ggplot(dist.df, aes(x = tree, y = micro)
       )+
  geom_point(size =3, alpha = .7) +
  theme(text = element_text(size=20, color = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.margin = ggplot2::margin(.5,.5,.5,.5 , "cm"),
        axis.line = element_line(colour = "black"),
        axis.text = element_text(color = "black")) +
  xlab("Cophenetic Distance")+
  ylab("Bray Distance")



#png("/Users/cgaulke/unsynced_projects/comp_metagenomes/analysis/Figures/ko_phylosym.png",
#    res = 300, width = 8, height = 7,units = "in")
dist.plot
#dev.off()
#now mantel test

#pearson
set.seed(731)
mammal_ko_pearson.mantel <-
  mantel(phylo_ko_dist.mat,
         mammal.dist,
         permutations = 10000)
mammal_ko_pearson.mantel$statistic # 0.5992446
mammal_ko_pearson.mantel$signif #0.00839916

#spearman
mammal_ko_spearman.mantel <-
  mantel(phylo_ko_dist.mat,
         mammal.dist,
         permutations = 10000,
         method = "spearman")

mammal_ko_spearman.mantel$statistic # 0.2966302
mammal_ko_spearman.mantel$signif # 0.06139386

#kendal
mammal_ko_kendall.mantel <-
  mantel(phylo_ko_dist.mat,
         mammal.dist,
         permutations = 10000,
         method = "kendall")

mammal_ko_kendall.mantel$statistic # 0.2457998
mammal_ko_kendall.mantel$signif # 0.04819518


# FUNCTIONAL PHYLOSYMBIOSIS: Pagels Module ----------------------------------

phylo_module.gm <-
  phylo_module.gm[mammal.tre$tip.label,]

pname.vec   <- NULL
plambda.vec <- NULL
ppval.vec   <- NULL
plogl.vec   <- NULL
plogl0.vec  <- NULL

for(i in 1:ncol(phylo_module.gm)){

  tree.data <- as.numeric(phylo_module.gm[,i])
  names(tree.data) <- mammal.tre$tip.label

  tree.test <- phylosig(mammal.tre,tree.data, method = "lambda", test = TRUE)

  pname.vec   <- c(pname.vec, colnames(phylo_module.gm)[i])
  plambda.vec <- c(plambda.vec, tree.test$lambda)
  ppval.vec   <- c(ppval.vec,   tree.test$P)
  plogl.vec   <- c(plogl.vec,   tree.test$logL)
  plogl0.vec  <- c(plogl0.vec,  tree.test$logL0)
}

pagels_module.df <- NULL
pagels_module.df$names  <- pname.vec
pagels_module.df$lambda <- as.numeric(plambda.vec)
pagels_module.df$pval   <- as.numeric(ppval.vec)
pagels_module.df$logl   <- as.numeric(plogl.vec)
pagels_module.df$logl0  <- as.numeric(plogl0.vec)
pagels_module.df$qval   <- as.numeric(qvalue::qvalue(as.numeric(pagels_module.df$pval))$qvalue)
pagels_module.df$fdr    <- as.numeric(p.adjust(pagels_module.df$pval,method = "fdr"))

pagels_module.df <- as.data.frame(pagels_module.df)
pagels_module.df$names <- gsub(pattern = "md\\.", replacement = "",pagels_module.df$names)
rownames(pagels_module.df) <- pagels_module.df$names

keeps_modules.pagels <- rownames(pagels_module.df[which(pagels_module.df$fdr < .05),])

#read in module names
module_names.df <- read.table("/Users/cgaulke/unsynced_projects/comp_metagenomes/data/Kegg_data/module_names_parsed.txt",
                              sep = "\t", row.names = 1, quote = "")


keeps_modules_pagel.df <- pagels_module.df[keeps_modules.pagels, ]

keeps_modules_pagel.df$long_name <- module_names.df[keeps_modules.pagels,]


#make denro for the heatmap
hc <- as.hclust.phylo(mammal.tre)
hc <- as.dendrogram(hc)


cols <- brewer.pal(5, "RdBu")
pal <- colorRampPalette(cols)
colnames(phylo_module.gm) <- gsub(pattern = "md\\.", replacement = "",colnames(phylo_module.gm) )

#Make heatmap
tpagel_module.heatmap <- heatmap.2(
  (as.matrix(phylo_module.gm[,keeps_modules.pagels])),
  trace = "none",
  col = pal(100),
  scale = "column",
  Rowv = hc,
  dendrogram = "row",
  margins = c(17,10),
  labCol = FALSE
)



