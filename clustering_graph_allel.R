library("igraph")
library("data.table")
library("stringdist")
library("Biostrings")
library("dplyr")
library("tibble")
library("ggplot2")
library(RColorBrewer)

#### choose distance - Levestein or Hamming ####
method_clustering = "hamming" # lv" #### or "hamming"
igraph_from_data_cdr3.2<-function(df,max_errs=1) {
  df = df %>% group_by(cdr3aa,v,patient) %>% as.data.table()
  seqs = df$cdr3aa
  graph<-graph.empty(n = length(seqs), directed=F)
  tmp<-stringdistmatrix(seqs,seqs,method=method_clustering)
  graph<-add.edges(graph, t(which(tmp<=max_errs,arr.ind=T)))
  graph = graph %>% set_edge_attr('weight', value = tmp[tmp<=max_errs])
  graph<-igraph::simplify(graph)
  graph<-set.vertex.attribute(graph, 'label', V(graph), seqs)
  graph<-set.vertex.attribute(graph, 'vgene', V(graph), df$v)
  graph<-set.vertex.attribute(graph, 'patient', V(graph), df$patient)
  return(graph)
}

df <- data.frame(CD=character(),
                 Allel=character(),
                 Distance=character(),
                 Method=character(),
                 Number_all_clones=integer(),
                 Number_clones_in_clusters=integer(),
                 Fract_clustered=integer(),
                 Number_of_clusters=integer(),
                 Average_size_of_cluster=character(),
                 Unique_CDR3aa=logical(),
                 Members=character(),
                 stringsAsFactors=FALSE)
cds = c('CD4', 'CD8', 'UNDEF')
uniqe_or_not = c(TRUE, FALSE)
allel_list_cd4 = c('drb1_07_01', 'drb1_15_01', 'drb1_11_01', 'drb3_02_02', 'drb4_01_03', 'drb5_01_01', 'dqa1_01_02', 'dqa1_02_01', 'dqa1_05_05', 'dqb1_03_01',
                   'dqb1_06_02', 'dqb1_02_02', 'dqb1_05_01', 'dpa1_01_03', 'dpb1_04_01', 'dpb1_04_02')
allel_list_cd8 = c('a01', 'a02','a03','a25','a24', 'b07', 'c04','c06', 'c07')
allel_list_undef = c(allel_list_cd4, allel_list_cd8)
distances = c(0,1,2)
for (cd in cds) {
  if (cd == 'CD4') {
    allel_list = allel_list_cd4
  } else if (cd == 'CD8') {
    allel_list = allel_list_cd8
  } else if (cd == 'UNDEF') {
    allel_list = allel_list_undef
  }
  for (allel in allel_list) {
    res_path = paste('~/../../Volumes/LTI/TeamFolder/Data/Petrovax/data_analysis/final_results/venn_used_6m_info/homology/allel_specific/',tolower(cd),'/', sep='')
    fordb = paste(res_path ,tolower(allel),'.txt', sep='')
    db = fread(fordb)
    print(nrow(db))
    for (unique in uniqe_or_not) {
      if (unique == TRUE) {
        ant_ori = db %>% select(cdr3aa,v,patient) %>% distinct_at(vars(cdr3aa,patient), .keep_all = TRUE)
      } else {
        ant_ori = db %>% select(cdr3aa,v,patient) #%>% unique()
      }
      print(nrow(ant_ori))
      for (max_dist in distances) {

        for_gr = igraph_from_data_cdr3.2(ant_ori,max_dist)
        for_gr_gs = delete.vertices(simplify(for_gr), degree(for_gr)==0)
        E(for_gr_gs)$weight[E(for_gr_gs)$weight==0]<-0.5
        comp_gs = components(for_gr_gs)
        df.allel <- data.frame(cd, allel, max_dist,    method_clustering, nrow(ant_ori), gorder(for_gr_gs), gorder(for_gr_gs)/nrow(ant_ori),
                               comp_gs$no, mean(comp_gs$csize),
                               unique)
        names(df.allel) <- c('CD', 'Allel', 'Distance', 'Method', 'Number_all_clones', 'Number_clones_in_clusters', 'Fract_clustered',
                             'Number_of_clusters', 'Average_size_of_cluster',
                                     'Unique_CDR3aa')
        memb.value <- as.numeric(c(comp_gs$membership))
        df.allel$Members <- list(memb.value)
        df <- rbind(df, df.allel)
        df.cluster <- as_long_data_frame(for_gr_gs)
        res_tables <- paste(res_path,'tables/',sep='')
        dir.create(file.path(res_tables), showWarnings = FALSE)
        outfile_df_cluster = paste(res_tables,allel, '_', max_dist,'_', unique, '.txt', sep = "")

        fwrite(df.cluster,file = outfile_df_cluster,sep="\t",quote = F)

      }
    }
  }
}

outfile_df = '~/../../Volumes/LTI/TeamFolder/Data/Petrovax/data_analysis/final_results/venn_used_6m_info/homology/allel_specific/all_clusters_data_with_members.txt'
fwrite(df,file = outfile_df,sep="\t",quote = F)

cd='CD4'
res_path = paste('~/../../Volumes/LTI/TeamFolder/Data/Petrovax/data_analysis/final_results/venn_used_6m_info/homology/allel_specific/',tolower(cd),'/', sep='')

allel = 'DPB1_04_02'
to_tit = paste(cd, allel, sep=' ')
max_dist = 2
type_dist= 'max Hamming 2' #'Public clones'#

#### table should be in working directory in VDJdb format, source (eg patient's id) of CDR3 should be in 'meta.subject.id' column ####
fordb = paste(res_path ,tolower(allel),'.txt', sep='')
db = fread(fordb)

#### select layout for graphs ####
layout_for_graphs = layout_nicely# layout_with_kk #layout_on_sphere #layout_nicely #lay#out_with_fr #layout_with_gem #layout_with_mds #layout_with_graphopt #layout_as_star #layout_on_grid #layout_as_tree

#### select your target antigen (it should be in 'antigen.gene' column) ####
ant_name = "FRDYVDRFYKTLRAEQASQE"

#### alpha ####

ant_ori = db %>% select(cdr3aa,v,patient) #%>% unique()
forpng = paste(res_path,allel, '_', max_dist, '_', method_clustering, '_V_gene.png', sep = "")
forpng2 = paste(res_path,allel, '_', max_dist, '_', method_clustering, '_source.png', sep = "")
#fortitle = paste(ant_name, 'beta vgene')
#fortitle1 = paste(ant_name, 'beta source')
### selects top 8 V-genes ###
top_vg = ant_ori %>% group_by(v) %>% summarise(v_n = n()) %>% as.data.table() %>% arrange(desc(v_n)) %>% head(14)
if (nrow(ant_ori %>% group_by(v) %>% summarise(v_n = n()) %>% as.data.table()) > 15) {
  ant = ant_ori %>% mutate(v = ifelse(v %in% top_vg$v,v,"others")) %>% as.data.table()
} else {
  ant = ant_ori
}

### max Levenstein/Hamming distance = second argument in igraph function ###
gs_ant = igraph_from_data_cdr3.2(ant,max_dist)
gs = delete.vertices(simplify(gs_ant), degree(gs_ant)==0)
E(gs)$weight[E(gs)$weight==0]<-0.5

for_gr = igraph_from_data_cdr3.2(ant_ori,max_dist)
for_gr_gs = delete.vertices(simplify(for_gr), degree(for_gr)==0)
E(for_gr_gs)$weight[E(for_gr_gs)$weight==0]<-0.5

### plot V-genes ###

coul2 <- rainbow(length(unique(V(gs)$vgene))-1, alpha=.5)
clr = c('#A9A9A9', coul2)
my_color=ifelse(V(gs)$vgene == 'others', "#A9A9A9", clr[as.numeric(as.factor(V(gs)$vgene))])

V(gs)$color = my_color#[as.numeric(as.factor(V(gs)$vgene))]
legend_cats = data.frame(attr = unique(vertex_attr(gs, "vgene")),
                         color = unique(V(gs)$color))
legend_cats  = legend_cats[order(legend_cats$attr), c(1, 2)]
png(forpng, 700, 700)
set.seed(42)
simplify(gs)
if (max_dist == 0) {
  e_label = NA
  e_width = 5/(E(gs)$weight+1)
} else {
  e_label = E(gs)$weight/2
  e_width = 15/(E(gs)$weight+1)
}
plot.igraph(gs, layout = layout_for_graphs, asp = 1,
            #main = fortitle,
            ## nodes =======================================
            vertex.label = NA,#V(gs)$vgene,
            vertex.size=7,
            ## edges =======================================
            edge.label = e_label,
            edge.label.cex = 1.3,
            edge.color = "lightgray",
            edge.width=e_width,
            edge.lty = "solid",           ## linetype: blank, solid, dashed, dotted,
            ## dotdash, longdash, or twodash
            edge.curved = 0.05#,            ## 0 to 1 or TRUE (0.5)
)
title(main=paste(to_tit, type_dist,sep='\n'))
legend_cats$color = as.character(legend_cats$color)
legend(x = "bottomleft",      ## position, also takes x,y coordinates
       legend = legend_cats$attr,
       pch = 19,              ## legend symbols see ?points
       col = legend_cats$color,
       bty = "n",
       title = "V genes")
dev.off()

### plot sources ###

coul2 <- rainbow(length(unique(V(gs)$patient)), alpha=.5)
clr = c(coul2)
my_color=ifelse(V(gs)$patient == 'others', "#A9A9A9", clr[as.numeric(as.factor(V(gs)$patient))])

V(gs)$color = my_color#[as.numeric(as.factor(V(gs)$vgene))]
legend_cats = data.frame(attr = unique(vertex_attr(gs, "patient")),
                         color = unique(V(gs)$color))
legend_cats  = legend_cats[order(legend_cats$attr), c(1, 2)]
set.seed(42)
png(forpng2, 700, 700)
plot.igraph(gs, layout = layout_for_graphs, asp = 1,
            #main = fortitle1,
            ## nodes =======================================
            vertex.label = NA,#V(gs)$patient,
            vertex.size=7,
            ## edges =======================================
            edge.label = e_label,
            edge.label.cex = 1.3,
            edge.color = "lightgray",
            edge.width=e_width,
            edge.lty = "solid",           ## linetype: blank, solid, dashed, dotted,
            ## dotdash, longdash, or twodash
            edge.curved = 0.05#,            ## 0 to 1 or TRUE (0.5)
            #vertex.color=coul2[as.numeric(as.factor(V(gs)$vgene))]
)
title(main=paste(to_tit, type_dist,sep='\n'))
legend_cats$color = as.character(legend_cats$color)
legend(x = "bottomleft",      ## position, also takes x,y coordinates
       legend = legend_cats$attr,
       pch = 19,              ## legend symbols see ?points
       col = legend_cats$color,
       bty = "n",
       title = "Sources")
dev.off()

df_clusters = as_long_data_frame(for_gr_gs)
outfile = paste(res_path,allel, '_', max_dist, '_', method_clustering, '.txt', sep = "")
fwrite(df_clusters,file = outfile,sep="\t",quote = F)

