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
  df = df %>% group_by(cdr3aa,v,patient,epi) %>% as.data.table()
  seqs = df$cdr3aa
  graph<-graph.empty(n = length(seqs), directed=F)
  tmp<-stringdistmatrix(seqs,seqs,method=method_clustering)
  graph<-add.edges(graph, t(which(tmp<=max_errs,arr.ind=T)))
  graph = graph %>% set_edge_attr('weight', value = tmp[tmp<=max_errs])
  graph<-igraph::simplify(graph)
  graph<-set.vertex.attribute(graph, 'label', V(graph), seqs)
  graph<-set.vertex.attribute(graph, 'vgene', V(graph), df$v)
  graph<-set.vertex.attribute(graph, 'patient', V(graph), df$patient)
  graph<-set.vertex.attribute(graph, 'epi', V(graph), df$epi)
  return(graph)
}


igraph_from_data_cdr3.2_no_epi<-function(df,max_errs=1) {
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
                 Epitope=character(),
                 Distance=character(),
                 Method=character(),
                 Number_all_clones=integer(),
                 Number_clones_in_clusters=integer(),
                 Fract_clustered=integer(),
                 Number_of_clusters=integer(),
                 Average_size_of_cluster=character(),
                 Unique_CDR3aa=logical(),
                 stringsAsFactors=FALSE)


uniqe_or_not = c(TRUE)#, FALSE)
distances = c(2)#0,1,2)
epitopes = c('FCN', 'GLR', 'IED', 'KCY', 'LPQ', 'LTD', 'QQL', 'QYI', 'RLQ', 'VRF', 'YLQ')
res_path = paste('~/../../Volumes/LTI/TeamFolder/Data/Petrovax/GE41/used 6m info/', sep='')
for (epitope in epitopes) {
fordb = paste(res_path ,'epitopes/',epitope,'.txt', sep='')
db = fread(fordb)
print(nrow(db))
if ((epitope == 'IED')|(epitope == 'QQL')) {
  cd = 'CD4'
} else {
    cd = 'CD8'
  }
if (nrow(db) > 3) {
for (unique in uniqe_or_not) {
  if (unique == TRUE) {
    ant_ori = db %>% select(cdr3aa,v,patient) %>% distinct_at(vars(cdr3aa,patient), .keep_all = TRUE)
  } else {
    ant_ori = db %>% select(cdr3aa,v,patient) #%>% unique()
  }
  print(paste('unique - ', nrow(ant_ori),sep=''))
  for (max_dist in distances) {
    for_gr = igraph_from_data_cdr3.2_no_epi(ant_ori,max_dist)
    for_gr_gs = delete.vertices(simplify(for_gr), degree(for_gr)==0)
    E(for_gr_gs)$weight[E(for_gr_gs)$weight==0]<-0.5
    comp_gs = components(for_gr_gs)
    df.allel <- data.frame(cd, epitope, max_dist,    method_clustering, nrow(ant_ori), gorder(for_gr_gs), gorder(for_gr_gs)/nrow(ant_ori),
                           comp_gs$no, mean(comp_gs$csize),
                           unique)
    names(df.allel) <- c('CD', 'Epitope', 'Distance', 'Method', 'Number_all_clones', 'Number_clones_in_clusters', 'Fract_clustered',
                         'Number_of_clusters', 'Average_size_of_cluster',
                                 'Unique_CDR3aa')
    df <- rbind(df, df.allel)


    gs_ant = igraph_from_data_cdr3.2_no_epi(ant_ori,max_dist)
    gs = delete.vertices(simplify(gs_ant), degree(gs_ant)==0)
    E(gs)$weight[E(gs)$weight==0]<-0.5

    if (gorder(gs) > 0) {
    forpng_vgen = paste(res_path,'epitopes/plots/epitope_',epitope,'_',max_dist,'_Vgen.png', sep = "")
    forpng_source = paste(res_path,'epitopes/plots/epitope_',epitope,'_',max_dist,'_Source.png', sep = "")

    top_vg = ant_ori %>% group_by(v) %>% summarise(v_n = n()) %>% as.data.table() %>% arrange(desc(v_n)) %>% head(14)
    if (nrow(ant_ori %>% group_by(v) %>% summarise(v_n = n()) %>% as.data.table()) > 15) {
      ant = ant_ori %>% mutate(v = ifelse(v %in% top_vg$v,v,"others")) %>% as.data.table()
    } else {
      ant = ant_ori
    }
    sub_tit = paste('Clustered ', gorder(gs), ' out of ', nrow(db), ' clones.')


    ### plot V-genes ###

    coul2 <- rainbow(length(unique(V(gs)$vgene))-1, alpha=.5)
    clr = c('#A9A9A9', coul2)
    my_color=ifelse(V(gs)$vgene == 'others', "#A9A9A9", clr[as.numeric(as.factor(V(gs)$vgene))])

    V(gs)$color = my_color#[as.numeric(as.factor(V(gs)$vgene))]
    legend_cats = data.frame(attr = unique(vertex_attr(gs, "vgene")),
                             color = unique(V(gs)$color))
    legend_cats  = legend_cats[order(legend_cats$attr), c(1, 2)]
    png(forpng_vgen, 700, 700)
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
    title(main=paste(epitope, type_dist, sep='\n'), sub = sub_tit)
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
    V(gs)$color = coul2[as.numeric(as.factor(V(gs)$patient))]
    legend_cats = data.frame(attr = unique(vertex_attr(gs, "patient")),
                             color = unique(V(gs)$color))
    legend_cats  = legend_cats[order(legend_cats$attr), c(1, 2)]
    set.seed(42)
    png(forpng_source, 700, 700)
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
    title(main=paste(epitope, type_dist,sep='\n'), sub = sub_tit)
    legend_cats$color = as.character(legend_cats$color)
    legend(x = "bottomleft",      ## position, also takes x,y coordinates
           legend = legend_cats$attr,
           pch = 19,              ## legend symbols see ?points
           col = legend_cats$color,
           bty = "n",
           title = "Sources")
    dev.off()
    }
  }
  }
}
}

outfile = paste(res_path,'epitopes/summary_num_clustes_avr_size.txt', sep = "")
fwrite(df,file = outfile,sep="\t",quote = F)


max_dists = c(2)
#max_dist = 0
for (max_dist in max_dists) {
if (max_dist != 0) {
  type_dist= paste('max Hamming ', max_dist, sep='')
} else {
  type_dist= 'Public clones'
}

to_tit = 'CD8'
#### table should be in working directory in VDJdb format, source (eg patient's id) of CDR3 should be in 'meta.subject.id' column ####
fordb = paste(res_path ,'all_epi_specific_with_', to_tit, '.txt', sep='')
db = fread(fordb)

#### select layout for graphs ####
layout_for_graphs = layout_nicely# layout_with_kk #layout_on_sphere #layout_nicely #lay#out_with_fr #layout_with_gem #layout_with_mds #layout_with_graphopt #layout_as_star #layout_on_grid #layout_as_tree

#### alpha ####
ant_ori = db %>% select(cdr3aa,v,patient,epi) #%>% unique()
forpng_epi = paste(res_path,'epitopes/plots/all_epitopes_',to_tit,'_',max_dist,'_Epitope.png', sep = "")
forpng_vgen = paste(res_path,'epitopes/plots/all_epitopes_',to_tit,'_',max_dist,'_Vgen.png', sep = "")
forpng_source = paste(res_path,'epitopes/plots/all_epitopes_',to_tit,'_',max_dist,'_Source.png', sep = "")

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

sub_tit = paste('Clustered ', gorder(gs), ' out of ', nrow(db), ' clones.')

### plot Epitopes-genes ###

coul2 <- rainbow(length(unique(V(gs)$epi)), alpha=.5)

V(gs)$color = ifelse(V(gs)$epi == to_tit, "#A9A9A9", clr[as.numeric(as.factor(V(gs)$epi))])#coul2[as.numeric(as.factor(V(gs)$epi))]
legend_cats = data.frame(attr = unique(vertex_attr(gs, "epi")),
                         color = unique(V(gs)$color))
legend_cats  = legend_cats[order(legend_cats$attr), c(1, 2)]
png(forpng_epi, 700, 700)
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
title(main=paste(to_tit, type_dist,sep='\n'), sub = sub_tit)
legend_cats$color = as.character(legend_cats$color)
legend(x = "bottomleft",      ## position, also takes x,y coordinates
       legend = legend_cats$attr,
       pch = 19,              ## legend symbols see ?points
       col = legend_cats$color,
       bty = "n",
       title = "Specificity")
dev.off()

### plot V-genes ###

coul2 <- rainbow(length(unique(V(gs)$vgene))-1, alpha=.5)
clr = c('#A9A9A9', coul2)
my_color=ifelse(V(gs)$vgene == 'others', "#A9A9A9", clr[as.numeric(as.factor(V(gs)$vgene))])

V(gs)$color = my_color#[as.numeric(as.factor(V(gs)$vgene))]
legend_cats = data.frame(attr = unique(vertex_attr(gs, "vgene")),
                         color = unique(V(gs)$color))
legend_cats  = legend_cats[order(legend_cats$attr), c(1, 2)]
png(forpng_vgen, 700, 700)
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
title(main=paste(to_tit, type_dist, sep='\n'), sub = sub_tit)
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
V(gs)$color = coul2[as.numeric(as.factor(V(gs)$patient))]
legend_cats = data.frame(attr = unique(vertex_attr(gs, "patient")),
                         color = unique(V(gs)$color))
legend_cats  = legend_cats[order(legend_cats$attr), c(1, 2)]
set.seed(42)
png(forpng_source, 700, 700)
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
title(main=paste(to_tit, type_dist,sep='\n'), sub = sub_tit)
legend_cats$color = as.character(legend_cats$color)
legend(x = "bottomleft",      ## position, also takes x,y coordinates
       legend = legend_cats$attr,
       pch = 19,              ## legend symbols see ?points
       col = legend_cats$color,
       bty = "n",
       title = "Sources")
dev.off()
}


df_clusters = as_long_data_frame(for_gr_gs)
outfile = paste(res_path,allel, '_', max_dist, '_', method_clustering, '.txt', sep = "")
fwrite(df_clusters,file = outfile,sep="\t",quote = F)

