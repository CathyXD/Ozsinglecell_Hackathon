cpg_processed <- NormalizeData(cpg_processed, ) %>% FindVariableFeatures() %>% ScaleData()
vargenes <- VariableFeatures(object = cpg_processed)[grepl(pattern = "^Ig[klh]", VariableFeatures(object = cpg_processed))]
cpg_processed <- RunPCA(cpg_processed, features = vargenes)
cpg_processed <- RunUMAP(cpg_processed, dims = 1:20)
cpg_processed <- FindNeighbors(cpg_processed)
cpg_processed <- FindClusters(cpg_processed)

ind_clone <- names(which(table(cpg_processed$clone_id)>=3))
ind_clone <- ind_clone[ind_clone != "No_contig"]

distance_calculation_pca <- function(srt, ind){
  df <- srt@reductions$PCA@cell.embeddings
  df <- df[colnames(srt)[srt$clone_id == ind],]
  temp <- sapply(1:19, function(dim1){
    x <- dim1
    y  <- dim1 +1 
    test <- sapply(seq_along(rownames(df)), function(r){
      x1 = df[,x][r]
      y1 = df[,y][r]
      vx = setdiff(df[,x], x1)
      vy = setdiff(df[,y],y1 ) 
      x_diff <- sapply(abs(vx - x1), function(x) x^2)
      y_diff <- sapply(abs(vy - y1), function(x) x^2)
      sqr_dist = x_diff + y_diff
      return(sum(sqrt(sqr_dist))/length(vx))
    })
    return(structure(test, names = rownames(df)))
  })
  return(temp)
}


perclone_test <- lapply(ind_clone, function(ind){
  temp <- distance_calculation_pca(cpg_processed, ind)
  return(rowMeans(temp))
})
names(perclone_test) <- ind_clone


df <- rbindlist(lapply(perclone_test, function(x) as.data.frame(x) %>% `colnames<-`("dist")%>% mutate(rep = substr(names(x), 1, 3))), idcol = "clonetype")
head(df)
## ignore replicates 
df <- df %>% group_by(clonetype) %>% mutate(median = median(dist))
df$clonetype <- factor(df$clonetype, levels  = unique(df$clonetype[order(df$median, decreasing = T)]))
df <- df[order(df$clonetype, decreasing = F), ]

## count in replicates
df <- df %>% group_by(clonetype, rep) %>%  mutate(median = median(dist))
df$clonetype <- factor(df$clonetype, levels  = unique(df$clonetype[order(df$median, decreasing = T)]))
df <- df[order(df$clonetype, decreasing = F), ]

median_df <- df[, c("clonetype", "median", "rep")] %>% distinct()
median_df$size <- sapply(perclone_test, length)[median_df$clonetype]
median_df$category <- ifelse(median_df$size<8, "3-7", ifelse(median_df$size>7 & median_df$size<13 ,"8-12", "13-21") )

ggplot(median_df,aes(x = median, colour = category)) + geom_density(linetype = rep)


DimPlot(cpg_processed, reduction = "UMAP", cells.highlight = colnames(cpg_processed)[cpg_processed$clone_id=="B_107_6_7_173_1_5"])+
  ggtitle("B_107_6_7_173_1_5")

hist(as.numeric(perclone_test[["B_39_2_1_147_1_1"]]), prob = TRUE, breaks = 50, main = paste("B_39_2_1_147_1_1", "\nmedian at", median(perclone_test[["B_39_2_1_147_1_1"]])))
lines(density(perclone_test[["B_39_2_1_147_1_1"]]), col = 4, lwd = 2)

ggplot(df %>% filter(clonetype %in% head(levels(df$clonetype), 10)), aes(x = clonetype, y = dist))+ 
  geom_boxplot()+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  ggtitle("Top 10 clones with most heterogeneity")

ggplot(df %>% filter(clonetype %in% tail(levels(df$clonetype), 10)), aes(x = clonetype, y = dist))+ 
  geom_boxplot()+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  ggtitle("Top 10 clones with most homogeneity")

ggplot(rbind(df %>% filter(clonetype %in% head(levels(df$clonetype), 10)), df %>% filter(clonetype %in% tail(levels(df$clonetype), 10))), aes(x = clonetype, y = dist))+ 
  geom_boxplot()+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  ggtitle("Top and bottom 10 clones with most heterogeneity")

new_cells <- cpg_processed@meta.data %>% filter(clone_id != "No_contig") %>% rownames_to_column("cell_id") %>% select(cell_id, clone_id)
new_cells  %>% group_by(clone_id) %>% slice_sample(n = 1)

median_df <- df[, c("clonetype", "median")] %>% distinct()
median_df$size <- sapply(perclone_test, length)[median_df$clonetype]
median_df$category <- ifelse(median_df$size<8, "3-7", ifelse(median_df$size>7 & median_df$size<13 ,"8-12", "13-21") )

f1 <- ggplot(median_df,aes(x = median, colour = category)) + geom_density()



distance_calculation_pca <- function(srt, new_cells, n, perclone_test){
  df <- srt@reductions$PCA@cell.embeddings
  control_distribution <- sapply(3:21, function(n){
    controls <- new_cells  %>% group_by(clone_id) %>% slice_sample(n = 1) %>% pull(cell_id)
    controls <- sample(controls, n)
    subdf <- df[controls,]
  temp <- sapply(1:19, function(dim1){
    x <- dim1
    y  <- dim1 +1 
    test <- sapply(seq_along(rownames(subdf)), function(r){
      x1 = subdf[,x][r]
      y1 = subdf[,y][r]
      vx = setdiff(subdf[,x], x1)
      vy = setdiff(subdf[,y],y1 ) 
      x_diff <- sapply(abs(vx - x1), function(x) x^2)
      y_diff <- sapply(abs(vy - y1), function(x) x^2)
      sqr_dist = x_diff + y_diff
      return(sum(sqrt(sqr_dist))/length(vx))
    })
    return(structure(test, names = rownames(subdf)))
  })
  return(rowMeans(temp))
  })
}

df <- srt@reductions$PCA@cell.embeddings
control_distribution <- sapply(as.numeric(sapply(perclone_test, length)), function(n){
  controls <- sample(colnames(cpg_processed)[cpg_processed$clone_id!="No_contig"], n)
  subdf <- df[controls,]
  temp <- sapply(1:19, function(dim1){
    x <- dim1
    y  <- dim1 +1 
    test <- sapply(seq_along(rownames(subdf)), function(r){
      x1 = subdf[,x][r]
      y1 = subdf[,y][r]
      vx = setdiff(subdf[,x], x1)
      vy = setdiff(subdf[,y],y1 ) 
      x_diff <- sapply(abs(vx - x1), function(x) x^2)
      y_diff <- sapply(abs(vy - y1), function(x) x^2)
      sqr_dist = x_diff + y_diff
      return(sum(sqrt(sqr_dist))/length(vx))
    })
    return(structure(test, names = rownames(subdf)))
  })
  return(rowMeans(temp))
})


names(control_distribution) <- paste("control", names(perclone_test), sep = "_")
contr_df <- rbindlist(lapply(control_distribution, function(x) as.data.frame(x) %>% `colnames<-`("dist")), idcol = "control")
contr_df <- contr_df %>% group_by(control) %>% mutate(median = median(dist))
contr_df$control <- factor(contr_df$control, levels  = unique(contr_df$control[order(contr_df$median, decreasing = T)]))
contr_df <- contr_df[order(contr_df$control, decreasing = F), ]

ggplot(head(contr_df,10), aes(x = control, y = dist))+ 
  geom_boxplot()+ 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  ggtitle("control clone distance distribution")


median_contr_df <- contr_df[, c("control", "median")] %>% distinct()
median_contr_df$size <- sapply(perclone_test, length)[median_contr_df$control]
median_contr_df$category <- ifelse(median_contr_df$size<8, "3-7", ifelse(median_contr_df$size>7 & median_contr_df$size<13 ,"8-12", "13-21") )

median_df$group <- "real"
median_contr_df$group <- "control"
colnames(median_contr_df) <- colnames(median_df)
all_median_df <- rbind(median_df, median_contr_df)


f2 <- ggplot(all_median_df,aes(x = median, colour = category, linetype = group)) + geom_density()
f2+ggtitle("real vs control intraclone heterogeneity")
