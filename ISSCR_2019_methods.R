# package loading, installation----
library(devtools)
devtools::install_github("dynverse/dyno", force=TRUE, ref = 'master')
devtools::install_github("dynverse/dynmethods", force=TRUE, ref = 'master')
devtools::install_github("dynverse/dynplot", force=TRUE, ref = 'master')
devtools::install_github("dynverse/dynutils", force=TRUE, ref = 'master')
devtools::install_github("dynverse/dynwrap", force=TRUE, ref = 'master')
devtools::install_github("dynverse/dynguidelines", force=TRUE, ref = 'master')
devtools::install_github("dynverse/dynfeature", force=TRUE, ref = 'master')


library(dyno)
library(tidyverse)
library(data.table)
options(rgl.useNULL=TRUE)
.rs.restartR()
library("plot3Drgl")


# load raw data ----
load("/Volumes/ncatssctl/NGS_related/ddSeq/Analysis_across_runs/R_scripts/Seurat/Aggr6.Seurat.Object.RData")

# set up task ----

mycolanno <- data.table::fread('/Volumes/ncatssctl/NGS_related/ddSeq/Analysis_across_runs/R_scripts/ddseq-mycolanno.csv')
mycolanno <- subset(mycolanno, mycolanno$Barcode %in% names(as.data.frame(as.matrix(aggr6@data))))
names(mycolanno) <- c('cell_ids', 'cell_info')

grouping <- as.data.frame(mycolanno)
row.names(grouping) <- mycolanno$cell_ids
grouping <- subset(grouping, select=c('cell_info'))
head(grouping)
grouping <- as.vector(grouping)
grouping <- setNames(as.character(grouping$cell_info), row.names(grouping))

subsample.genes <- aggr6@var.genes

subsample.genes <- unique(c(subsample.genes, 'ESRG', 'TPBG', 'PAX6', 'DUSP6', 'NANOG', 'POU5F1', 'PRTG'))

aggr6.subset.raw <- as.matrix(aggr6@raw.data)[subsample.genes,]
aggr6.subset.data <- as.matrix(aggr6@data)[subsample.genes,]



task <- wrap_expression(
  counts = t(aggr6.subset.raw),
  expression= t(aggr6.subset.data),
  cell_info = mycolanno
)

# run TI methods -----
model.projectedslingshot <- infer_trajectory(task, 'projected_slingshot')
model.celltree <- infer_trajectory(task, 'celltree_maptpx')
model.scorpius <- infer_trajectory(task, 'scorpius')
model.tscan <- infer_trajectory(task, 'tscan') # works but it labels day 2 as endpoint NPC
model.paga <- infer_trajectory(task, 'paga', verbose=T)

model.embeddr <- infer_trajectory(task, 'embeddr') # ok, it is linear

# plot from each method----

model <- model.embeddr
model <- model.projectedslingshot
model <- model.tscan

model <- model %>% add_root_using_expression(c("DUSP6"), task$expression)

model <- label_milestones_markers(
  model,
  markers = list(
    Neuroprogenitor = c("PAX6", "PRTG", "TPBG"),
    Pluripotent=c('DUSP6','NANOG', 'ESRG')
  ),
  task$expression
)

task$grouping <- grouping
model$grouping <- grouping
model <- model %>% add_dimred(dyndimred::dimred_mds, expression_source = task$expression)

plot_dimred(
  model, 
  expression_source = task$expression, 
  grouping = grouping,
  groups=grouping,
  color_cells="auto"
)

dynplot::plot_topology(model)

plot_dendro(
  model, 
  expression_source = task$expression,
  grouping = grouping,
  groups=grouping,
  color_cells="auto"
)

plot_dimred(
  model, 
  expression_source = task$expression,
  grouping = grouping,
  color_cells="milestone"
)

plot_dimred(
  model, 
  expression_source = task$expression,
  grouping = grouping,
  color_cells="pseudotime"
)

plot_heatmap(
  model,
  expression_source = task$expression,
  grouping = task$grouping,
  features_oi = 50
)

# writing out pseudotime values for GSEA use----
pseudotime <- as.data.frame(model$pseudotime)
pseudotime <- as.data.frame(t(pseudotime))

fwrite(pseudotime, '/Volumes/ncatssctl/NGS_related/ddSeq/Analysis_across_runs/GSEA/Aggr6.Embeddr.phenotype.pseudotime.cls', row.names = F, col.names = F, sep='\t', quote=F)
#

all(row.names(task$counts) == row.names(model$pseudotime))
all(row.names(as.matrix(aggr6@raw.data)) == row.names(model$pseudotime))
task$counts[1:5,1:5]
as.matrix(aggr6@raw.data)[1:5,1:5]


# make combined table for bootstrapping pseudotimes -----

# weight pseudotime from each method----

# input to randomforest regression -----


