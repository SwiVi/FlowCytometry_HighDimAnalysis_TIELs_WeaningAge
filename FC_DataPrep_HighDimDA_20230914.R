# Load packages ----

# Will also have to download several package dependencies, such as data.table, flowCore....

library(Spectre)
library(Rtsne)
library(cowplot)
library(ggplot2)
library(SingleCellExperiment)
library(miloR)
library(Seurat)
library(edgeR)
library(dplyr)
library(viridis)
library(data.table)

###################### 
# DATA PREPROCESSING # ----
###################### 

# Set a working directory for output plots ----

dir.create(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'))  # create a directory; Sys.Date() will tack today's date to the end of the directory name
setwd(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_')) # set working directory
dir.create('Plots') # create a subdirectory
setwd('Plots') # set working directory; we have steps in workflow where plots will be stored by default; this is where they will now go

# Import & prep data ----

# Using the simple discovery workflow from Spectre: https://immunedynamics.io/spectre/simple-discovery/#2_Import_and_prep_data

# .... Read in files ----

data.list <- read.files(file.loc = '/home/Jayne.Wiarda/FS24/FCS_CD3Epos',
                        file.type = ".fcs",
                        do.embed.file.names = TRUE)

# .... Do a simple QC check ----

check <- do.list.summary(data.list)
check$name.table # check column names & values
check$ncol.check # check number of columns per sample
check$nrow.check # check number of rows (cells) per sample
data.list[[1]] # see first and last few rows of first sample

rm(check)

# .... Merge files ----

cell.dat <- do.merge.files(dat = data.list) # merge all files
cell.dat # view some cells from first and last file

rm(data.list)

# .... Rename columns ----

colnames(cell.dat) # note that column names just give channel designation and don't tell us about our cell marker IDs
colnames(cell.dat) <- c('FSC_A', 'SSC_A', 'CD8a_APC', 'CD27_FITC', 
                        'CD4_PerCPCy5.5', 'CD16_BUV395', 'MHCII_BUV496',
                        'CD45RC_BV421', 'gdTCR_BV510', 'CD2_BV711',
                        'CD8b_PE', 'CD3e_PECy7', 'FileName', 'FileNo') # rename some of the columns to incorporate cell marker ID + fluorophore ID
colnames(cell.dat) # check new column names

# Add metadata ----

# .... Read in the meta data ----

meta.dat <- fread('/home/Jayne.Wiarda/FS24/metadata/FS24meta.csv') # read in metadata file
meta.dat # look at first and last rows of metadata
colnames(meta.dat) # show all column names
meta.dat <- subset(meta.dat, select = -c(V1,V16)) # remove any unwanted metadata columns

# .... Add meta data to cell data frame ----

cell.dat <- do.add.cols(cell.dat, "FileName", meta.dat, "FileName", rmv.ext = TRUE) # add meta data
cell.dat # check that meta data was added

rm(meta.dat)

# Transform data ----

# Using the simple discovery workflow from Spectre: https://immunedynamics.io/spectre/simple-discovery/#2_Import_and_prep_data
# Using the data transformation workflow from Spectre: https://wiki.centenary.org.au/display/SPECTRE/Data+transformation

# .... Identify columns of data to transform ----

as.matrix(names(cell.dat)) # identify column numbers we want to apply data transformations to
to.asinh <- names(cell.dat)[c(3:12)] # specify column numbers we want to apply data transformations to
to.asinh # list columns we want to apply data transformations to

# .... Perform asinh transformations with a range of cofactors ----

cell.dat <- do.asinh(cell.dat,  # perform arcsinh transformation
                     to.asinh,# specify columns to transform
                     append.cf = TRUE, # appends cofactor value to end of transformation suffix
                     cofactor = 100) # specify a cofactor; flow cytometry (NOT mass cytometry) data typically needs cofactor between 100 to 1000...test different cofactors on different markers/channels to optimize

cell.dat <- do.asinh(cell.dat,  # perform arcsinh transformation
                     to.asinh,# specify columns to transform
                     append.cf = TRUE, # appends cofactor value to end of transformation suffix
                     cofactor = 250) # specify a cofactor; flow cytometry (NOT mass cytometry) data typically needs cofactor between 100 to 1000...test different cofactors on different markers/channels to optimize

cell.dat <- do.asinh(cell.dat,  # perform arcsinh transformation
                     to.asinh,# specify columns to transform
                     append.cf = TRUE, # appends cofactor value to end of transformation suffix
                     cofactor = 500) # specify a cofactor; flow cytometry (NOT mass cytometry) data typically needs cofactor between 100 to 1000...test different cofactors on different markers/channels to optimize

cell.dat <- do.asinh(cell.dat,  # perform arcsinh transformation
                     to.asinh,# specify columns to transform
                     append.cf = TRUE, # appends cofactor value to end of transformation suffix
                     cofactor = 750) # specify a cofactor; flow cytometry (NOT mass cytometry) data typically needs cofactor between 100 to 1000...test different cofactors on different markers/channels to optimize

cell.dat <- do.asinh(cell.dat,  # perform arcsinh transformation
                     to.asinh,# specify columns to transform
                     append.cf = TRUE, # appends cofactor value to end of transformation suffix
                     cofactor = 1000) # specify a cofactor; flow cytometry (NOT mass cytometry) data typically needs cofactor between 100 to 1000...test different cofactors on different markers/channels to optimize

# Batch correction ----

# I want to emphasize that not all batch correction methods will work well on every dataset, and not every dataset is designed properly to use every batch correction method! 
# I show the batch correction method I found to work well for this dataset, but DO YOUR RESEARCH AND TRY SEVERAL METHODS! 
# Also use preplanning to minimize batch effects as much as possible when running experiments and try to always utilize proper batch controls!!!
# Here I test multiple asinh transformations and compare results both before and after batch correction. The method I use (CytofBatchAdjust) performs batch correction on each parameter individually so I can do this; this is not the case for every batch correction method, so be aware that you may need to optimize to a single transformation for each parameter before using certain batch correction methods!

# ... Identify a reference sample ----

# Extract data for only a single reference sample and make into another dataset; we will use Ctrl03.

ref.dat <- do.filter(cell.dat, 
                     use.col = 'Sample', 
                     values = 'Ctrl03') # extract Ctrl3 samples and create new dataset

# .... Run coarse alignment ----

# Use coarse alignment protocol that is part of CytofBatchAdjust with function from Spectre: https://wiki.centenary.org.au/pages/viewpage.action?pageId=157272386#id-?2.Batchalignment-Initial(pre-alignment)plots

crs.dat <- cell.dat
as.matrix(names(crs.dat))
crs.dat <- run.align(ref.dat = ref.dat, # specify batch control dataset
                     target.dat = crs.dat, # specify target dataset
                     batch.col = 'Date', # specify column with batch labels
                     align.cols = colnames(crs.dat[,c(1:2, 29:78)]), # specify columns with cell markers to align; we opt to also batch-normalize FSC/SSC values even though they weren't transformed
                     method = '90p', # specify a method to use for normalization; should try multiple methods out to optimize (e.g.: SD, quantile, 50p, 75p, 90p, 95p tried on this dataset)
                     append.name = '_aligned', # specify suffix to append to new batch-aligned column names
                     dir = paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_')) # set a directory
as.matrix(names(crs.dat)) # check for new column names
crs.dat # check data

# .... Save data ----

dir.create(paste(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'), 'Data', sep = '/')) # create a directory to store data in
fwrite(crs.dat, # write data into the new directory as .csv file
       paste(paste(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'), 'Data', sep = '/'), 'BatchCorrectedData_AllTransformations.csv', sep = '/'))

cell.dat <- crs.dat # reassign batch-aligned data to cell.dat for further work
rm(crs.dat, ref.dat)

# .... Plot data transformations with different cofactors ----

# Showing plotting example for just one parameter. Also viewed transformations for additional parameters in order to optimize each.

var <- to.asinh[1] # change the [1] to view transformation of other parameters found in to.asinh
a <- make.colour.plot(do.subsample(cell.dat, 
                                   targets = rep(100, length(unique(cell.dat$FileName))), # plot only 100 cells from each sample to minimize plotting time
                                   divide.by = 'FileName'), 
                      paste0(var, '_asinh_cf100'), # x-axis parameter; could add '_aligned' suffix to look at distribution of batch-corrected data if desired
                      'SSC_A') # y-axis parameter
b <- make.colour.plot(do.subsample(cell.dat, 
                                   targets = rep(100, length(unique(cell.dat$FileName))), 
                                   divide.by = 'FileName'), 
                      paste0(var, '_asinh_cf250'),
                      'SSC_A')
c <- make.colour.plot(do.subsample(cell.dat, 
                                   targets = rep(100, length(unique(cell.dat$FileName))), 
                                   divide.by = 'FileName'), 
                      paste0(var, '_asinh_cf500'),
                      'SSC_A')
d <- make.colour.plot(do.subsample(cell.dat, 
                                   targets = rep(100, length(unique(cell.dat$FileName))), 
                                   divide.by = 'FileName'), 
                      paste0(var, '_asinh_cf750'),
                      'SSC_A')
e <- make.colour.plot(do.subsample(cell.dat, 
                                   targets = rep(100, length(unique(cell.dat$FileName))), 
                                   divide.by = 'FileName'), 
                      paste0(var, '_asinh_cf1000'),
                      'SSC_A')
plot_grid(a, b, c, d, e)

rm(var, a, b, c, d, e)

# Plot data transformations with different cofactors before and after batch correction 

subCtrl <- cell.dat[grepl("^Ctrl", cell.dat$Sample),]
colnames(subCtrl)
sub <- do.subsample(subCtrl, # create subsample of data
                    targets = rep(5000, length(unique(cell.dat$FileName))), # plot only 5000 cells from each sample to minimize plotting time & balance samples
                    divide.by = 'FileName')

setwd(paste(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'), 'Plots', sep = '/')) # set working directory
a <- ggplot(sub, aes(x= CD8b_PE_asinh_cf100, color = Date)) + # to look at other parameters, simply sub out 'CD8b_PE' for other parameters stored in to.asinh
        geom_density() + 
        xlim(c(-2,8)) +
        theme_bw() 
b <- ggplot(sub, aes(x= CD8b_PE_asinh_cf100_aligned, color = Date)) +
        geom_density() + 
        xlim(c(-2,8)) +
        theme_bw() 
c <- ggplot(sub, aes(x= CD8b_PE_asinh_cf250, color = Date)) +
        geom_density() + 
        xlim(c(-2,8)) +
        theme_bw() 
d <- ggplot(sub, aes(x= CD8b_PE_asinh_cf250_aligned, color = Date)) +
        geom_density() + 
        xlim(c(-2,8)) +
        theme_bw() 
e <- ggplot(sub, aes(x= CD8b_PE_asinh_cf500, color = Date)) +
        geom_density() + 
        xlim(c(-2,8)) +
        theme_bw() 
f <- ggplot(sub, aes(x= CD8b_PE_asinh_cf500_aligned, color = Date)) +
        geom_density() + 
        xlim(c(-2,8)) +
        theme_bw() 
g <- ggplot(sub, aes(x= CD8b_PE_asinh_cf750, color = Date)) +
        geom_density() + 
        xlim(c(-2,8)) +
        theme_bw() 
h <- ggplot(sub, aes(x= CD8b_PE_asinh_cf750_aligned, color = Date)) +
        geom_density() + 
        xlim(c(-2,8)) +
        theme_bw() 
i <- ggplot(sub, aes(x= CD8b_PE_asinh_cf1000, color = Date)) +
        geom_density() + 
        xlim(c(-2,8)) +
        theme_bw() 
j <- ggplot(sub, aes(x= CD8b_PE_asinh_cf1000_aligned, color = Date)) +
        geom_density() + 
        xlim(c(-2,8)) +
        theme_bw() 

plot_grid(a,b, c, d, e, f, g, h, i, j)

rm(a, b, c, d, e, f, g, h, i, j, to.asinh, sub, subCtrl)

# .... Select optimized cofactors to transform each parameter ----

as.matrix(names(cell.dat)) # list all column names
keep <- colnames(cell.dat[,c(1:2, 13:28, 35, 61, 69, 70, 72:74, 76:78, 79:80, 87, 113, 121, 122, 124:126, 128:130)]) # identify columns to keep; only keep one transformed column per cell marker
keep # check that we identified the proper columns to keep
cell.dat <- cell.dat[,c(1:2, 13:28, 35, 61, 69, 70, 72:74, 76:78, 79:80, 87, 113, 121, 122, 124:126, 128:130)] # subset to our optimized transformations
colnames(cell.dat) # see what column names are now left

rm(keep)

# Visualize data pre-alignment & post-alignment using t-SNE with only batch control samples ----

# .... Subset data ----

subCtrl <- cell.dat[grepl("^Ctrl", cell.dat$Sample),]
colnames(subCtrl)
sub <- do.subsample(subCtrl, # create subsample of data
                    targets = rep(750, length(unique(subCtrl$FileName))), # plot only 750 cells from each sample to minimize plotting time
                    divide.by = 'FileName')

# .... Calculate t-SNE coordinates ----

# Pre-alignment:

colnames(sub)
set.seed(123) # set a seed immediately before running Rtsne so that results are reproducible
tsne.out <- Rtsne(sub[,c(19:27)], # specify transformed cell markers to use for contructing t-SNE dimensionality reduction; note we opted to exclude CD3e since we have already gated cells based on CD3e expression (all cells CD3e+); however, CD3e expression level may be biologically relevant, so it wouldn't necessarily be wrong to also include CD3e as a clustering marker if we think the level has biological relevance. Regardless, we chose to omit but can still see if CD3e expression level corresponds to any other patterns in our dimensionally reduced data once plotted.
                  perplexity = 15) # toggle with perplexity parameter to adjust arrangement/spread of cells in plot; larger values will reduce cell spread; perplexity should be: 3 * perplexity < (# of subsampled cells) - 1; default starting value is 30
# alternative option is run.umap() to generate UMAP coordinates instead of t-SNE
sub$tsne1_PreAlignment <- tsne.out$Y[,1] # assign tSNE coordinates to cells in data subsample
sub$tsne2_PreAlignment <- tsne.out$Y[,2] # assign tSNE coordinates to cells in data subsample

rm(tsne.out)

# Post-alignment:

tsne.out <- Rtsne(sub[,c(31:39)], # specify transformed cell markers to use for contructing t-SNE dimensionality reduction; note we opted to exclude CD3e since we have already gated cells based on CD3e expression (all cells CD3e+); however, CD3e expression level may be biologically relevant, so it wouldn't necessarily be wrong to also include CD3e as a clustering marker if we think the level has biological relevance. Regardless, we chose to omit but can still see if CD3e expression level corresponds to any other patterns in our dimensionally reduced data once plotted.
                  perplexity = 15) # toggle with perplexity parameter to adjust arrangement/spread of cells in plot; larger values will reduce cell spread; perplexity should be: 3 * perplexity < (# of subsampled cells) - 1; default starting value is 30
# alternative option is run.umap() to generate UMAP coordinates instead of t-SNE
sub$tsne1_PostAlignment <- tsne.out$Y[,1] # assign tSNE coordinates to cells in data subsample
sub$tsne2_PostAlignment <- tsne.out$Y[,2] # assign tSNE coordinates to cells in data subsample

# .... Plot parameters on t-SNE ----

# Plot cell markers pre-alignment:

setwd(paste(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'), 'Plots', sep = '/')) # set working directory
as.matrix(names(sub))
make.multi.plot(sub, 
                'tsne1_PreAlignment', 
                'tsne2_PreAlignment', 
                colnames(sub[,c(1:2, 19:28)]), # Plot cell markers. Note: some of these parameters were not used to perform t-SNE dimensionality reduction, but we can still see where they show up in the plot
                colours = 'inferno',
                dot.size = 0.25,
                figure.title = 'TransformedCellMarkers_PreAlignment_750downsample_CtrlsOnly')

# Plot cell markers post-alignment:

setwd(paste(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'), 'Plots', sep = '/')) # set working directory
as.matrix(names(sub))
make.multi.plot(sub, 
                'tsne1_PostAlignment', 
                'tsne2_PostAlignment', 
                colnames(sub[,c(29:40)]), # Plot cell markers. Note: some of these parameters were not used to perform t-SNE dimensionality reduction, but we can still see where they show up in the plot
                colours = 'inferno',
                dot.size = 0.25,
                figure.title = 'TransformedCellMarkers_PostAlignment_750downsample_CtrlsOnly')

# Plot experimental parameters pre-alignment:

setwd(paste(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'), 'Plots', sep = '/')) # set working directory
as.matrix(names(sub))
make.multi.plot(sub, 
                'tsne1_PreAlignment', 
                'tsne2_PreAlignment', 
                c('Date', 'Sample'), # Plot selected experimental variables
                col.type = 'factor', # set all as qualitative over quantitative variables
                dot.size = 0.25,
                figure.title = 'Batches_PreAlignment_750downsample_CtrlsOnly')
# Can look at this plot to get an impression of if we have batch effects

# Plot experimental parameters pre-alignment:

setwd(paste(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'), 'Plots', sep = '/')) # set working directory
as.matrix(names(sub))
make.multi.plot(sub, 
                'tsne1_PostAlignment', 
                'tsne2_PostAlignment', 
                c('Date', 'Sample'), # Plot selected experimental variables
                col.type = 'factor', # set all as qualitative over quantitative variables
                dot.size = 0.25,
                figure.title = 'Batches_PostAlignment_750downsample_CtrlsOnly')
# Can look at this plot to get an impression of if we still have batch effects

rm(subCtrl, sub, tsne.out)

# Visualize data post-alignment using t-SNE with only experimental samples ----

# .... Subset data ----

cell.dat <- cell.dat[!grepl("Ctrl", cell.dat$Sample),] # remove ctrl samples from dataset for further work
cell.dat

sub <- do.subsample(cell.dat, # create subsample of data
                    targets = rep(500, length(unique(cell.dat$FileName))), # plot only 500 cells from each sample to minimize plotting time
                    divide.by = 'FileName')

# .... Calculate t-SNE coordinates ----

colnames(sub)
set.seed(123) # set a seed immediately before running Rtsne so that results are reproducible
tsne.out <- Rtsne(sub[,c(31:39)], # specify transformed cell markers to use for contructing t-SNE dimensionality reduction; note we opted to exclude CD3e since we have already gated cells based on CD3e expression (all cells CD3e+); however, CD3e expression level may be biologically relevant, so it wouldn't necessarily be wrong to also include CD3e as a clustering marker if we think the level has biological relevance. Regardless, we chose to omit but can still see if CD3e expression level corresponds to any other patterns in our dimensionally reduced data once plotted.
                  perplexity = 15) # toggle with perplexity parameter to adjust arrangement/spread of cells in plot; larger values will reduce cell spread; perplexity should be: 3 * perplexity < (# of subsampled cells) - 1; default starting value is 30
# alternative option is run.umap() to generate UMAP coordinates instead of t-SNE
sub$tsne1_PostAlignment <- tsne.out$Y[,1] # assign tSNE coordinates to cells in data subsample
sub$tsne2_PostAlignment <- tsne.out$Y[,2] # assign tSNE coordinates to cells in data subsample

# .... Plot parameters on t-SNE ----

# Plot cell markers

setwd(paste(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'), 'Plots', sep = '/')) # set working directory
as.matrix(names(sub))
make.multi.plot(sub, 
                'tsne1_PostAlignment', 
                'tsne2_PostAlignment', 
                colnames(sub[,c(29:40)]), # Plot cell markers. Note: some of these parameters were not used to perform t-SNE dimensionality reduction, but we can still see where they show up in the plot
                colours = 'inferno',
                dot.size = 0.25,
                figure.title = 'TransformedCellMarkers_PostAlignment_500downsample_ExperimentSamplesOnly')

# Plot experimental parameters

setwd(paste(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'), 'Plots', sep = '/')) # set working directory
as.matrix(names(sub))
make.multi.plot(sub, 
                'tsne1_PostAlignment', 
                'tsne2_PostAlignment', 
                colnames(sub[,c(6, 8, 9, 11, 12, 17, 18)]), # Plot selected experimental variables
                col.type = 'factor', # set all as qualitative over quantitative variables
                dot.size = 0.1,
                figure.title = 'Batches_PostAlignment_500downsample_ExperimentSamplesOnly')

rm(sub, tsne.out)

# Cluster cells ----

# Use FlowSOM clustering for quick ID of cell subsets, as done in Spectre dicovery workflow: https://immunedynamics.io/spectre/simple-discovery/#5_Clustering_and_DR

as.matrix(names(cell.dat))
cell.dat <- run.flowsom(cell.dat, 
                        colnames(cell.dat[,c(31:39)]), # use same cell markers to cluster as we've been using for dimensionality reduction
                        xdim = 12, # set number of x dimensions used in self-organizing map. Fewer dimensions needed for less complex datasets. Less dimensions will yield fewer clusters, as determined by x*y dimensions.
                        ydim = 12) # set number of y dimensions used in self-organizing map. Fewer dimensions needed for less complex datasets. Less dimensions will yield fewer clusters, as determined by x*y dimensions.
# Note: we most likely severely overclustered the data with 12x12 dimensions, so many clusters would need to be merged back together for appropriate biological interpretations, but overclustering does allow us to more precisely pinpoint inappropriate cell events by cluster IDs
cell.dat

# Visualize clustered cells ----

# .... Subset data ----

sub <- do.subsample(cell.dat, # create subsample of data
                    targets = rep(500, length(unique(cell.dat$FileName))), # plot only 500 cells from each sample to minimize plotting time
                    divide.by = 'FileName')

# .... Calculate t-SNE coordinates ----

colnames(sub)
set.seed(123) # set a seed immediately before running Rtsne so that results are reproducible
tsne.out <- Rtsne(sub[,c(31:39)], # specify transformed cell markers to use for contructing t-SNE dimensionality reduction; note we opted to exclude CD3e since we have already gated cells based on CD3e expression (all cells CD3e+); however, CD3e expression level may be biologically relevant, so it wouldn't necessarily be wrong to also include CD3e as a clustering marker if we think the level has biological relevance. Regardless, we chose to omit but can still see if CD3e expression level corresponds to any other patterns in our dimensionally reduced data once plotted.
                  perplexity = 15) # toggle with perplexity parameter to adjust arrangement/spread of cells in plot; larger values will reduce cell spread; perplexity should be: 3 * perplexity < (# of subsampled cells) - 1; default starting value is 30
# alternative option is run.umap() to generate UMAP coordinates instead of t-SNE
sub$tsne1_PostAlignment <- tsne.out$Y[,1] # assign tSNE coordinates to cells in data subsample
sub$tsne2_PostAlignment <- tsne.out$Y[,2] # assign tSNE coordinates to cells in data subsample

# .... Plot parameters on t-SNE ----

# Plot cell markers

setwd(paste(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'), 'Plots', sep = '/')) # set working directory
as.matrix(names(sub))
make.multi.plot(sub, 
                'tsne1_PostAlignment', 
                'tsne2_PostAlignment', 
                colnames(sub[,c(29:40)]), # Plot cell markers. Note: some of these parameters were not used to perform t-SNE dimensionality reduction, but we can still see where they show up in the plot
                colours = 'inferno',
                dot.size = 0.25,
                figure.title = 'TransformedCellMarkers_PostAlignment_500downsample_ExperimentSamplesOnly_ForClustering')

# Plot clusters

setwd(paste(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'), 'Plots', sep = '/')) # set working directory
as.matrix(names(sub))
make.multi.plot(sub, 
                'tsne1_PostAlignment', 
                'tsne2_PostAlignment', 
                'FlowSOM_cluster',
                col.type = 'factor',
                dot.size = 0.25,
                figure.title = 'ClustersFlowSOM_PostAlignment_500downsample_ExperimentSamplesOnly')
make.colour.plot(sub,
                 'tsne1_PostAlignment', 
                 'tsne2_PostAlignment', 
                 'FlowSOM_cluster',
                 col.type = 'factor',
                 dot.size = 0.25,
                 save.to.disk = FALSE, # don't save file
                 add.label = TRUE) # add labels to the clusters

rm(sub, tsne.out)

# .... Plot parameters on heatmap ----

setwd(paste(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'), 'Plots', sep = '/')) # set working directory
as.matrix(names(cell.dat))
exp <- do.aggregate(dat = cell.dat, 
                    use.cols = colnames(cell.dat[,c(29:40)]),
                    by = 'FlowSOM_cluster')
make.pheatmap(dat = exp, 
              sample.col = 'FlowSOM_cluster', 
              plot.cols = colnames(cell.dat[,c(29:40)]),
              standard.colours = 'inferno')
dev.off()

rm(exp)

# .... Save data ----

fwrite(cell.dat, 
       paste(paste(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'), 'Data', sep = '/'), 'BatchCorrectedData_ExperimentSamplesOnly_CellClustersPrefiltering.csv', sep = '/'))

# Filter cells ----

# .... Filter out non-TIEL events ----

# T-IEL = intraepithelial T cell

# cluster 109 appears to be CD2- and gdTCR-, indicating it is potentially non-T cell debri. This is supported by it having exceptionally low values for all other fluorescent intensities and forward/side scatter as well.

cell.dat <- subset(cell.dat, !(FlowSOM_cluster %in% c('109'))) # specify clusters to remove
table(cell.dat$FlowSOM_cluster) # check appropriate clusters were removed

# .... Save data ----

fwrite(cell.dat, 
       paste(paste(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'), 'Data', sep = '/'), 'BatchCorrectedData_ExperimentSamplesOnly_CellsFiltered.csv', sep = '/'))

# .... Save individual FCS files ----

as.matrix(names(cell.dat))
keep <- colnames(cell.dat[, c(3, 29:40) ])
keep
cell.dat <- cell.dat[, keep, with=FALSE] # whittle down to only the corrected cell markers and remove all meta data except FileName
as.matrix(names(cell.dat))

setwd(paste(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'), 'Data', sep = '/'))
dir.create('FCS_processed')
setwd('FCS_processed')
write.files(cell.dat, # save as individual FCS files in new directory created above
            file.prefix = "Processed",
            divide.by = 'FileName',
            write.csv = FALSE,
            write.fcs = TRUE)

rm(keep, cell.dat)

# Prep data for differential abundance analysis ----

# .... Read in data ----

# If starting from saved files, bring in last saved .csv with fread() function:
cell.dat <- fread(paste(paste(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'), 'Data', sep = '/'), 'BatchCorrectedData_ExperimentSamplesOnly_CellsFiltered.csv', sep = '/'))

# .... Subset data ----

min(table(cell.dat$FileName)) # see how many cells are left in smallest sample
#max(table(cell.dat$FileName)) # see how many cells are left in largest sample
sub <- do.subsample(cell.dat, # create subsample of data that we can more easily analyze; too many cells will inhibit some downstream analyses in miloR
                    targets = rep(min(table(cell.dat$FileName)), length(unique(cell.dat$FileName))), # extract as many cells from each sample as are found in our smallest sample
                    divide.by = 'FileName')
table(sub$FileName) # see same number of cells in each sample

rm(cell.dat) 

# .... Create Seurat object ----

# convert to Seurat object to calculate dimensionality reductions and some visualizations since the Seurat package is very user friendly for tinkering with single-cell data

as.matrix(names(sub)) # identify columns we want to subset
dat <- as.data.frame(t(sub[,c(31:39)])) # take only the cell markers we will use in downstream analyses
dim(dat)
rownames(dat) <- sub("^(([^_]*_){1}[^_]*).*", "\\1", rownames(dat)) # shorten row names
colnames(dat) <- paste((paste0('cell', 1:ncol(dat))), (paste0('sample', sub$Sample)), sep = '_') # add unique cell identifiers as column names
dat[1:9,1:9]
meta <- as.data.frame(sub[,c(3:18)]) # extract meta data
rownames(meta) <- colnames(dat)
head(meta)
seu <- CreateSeuratObject(dat)
dim(seu@assays$RNA@counts) # look at number of cell markers x number of total cells
seu <- AddMetaData(seu, meta) # add meta data

rm(sub, dat, meta)

# .... Calculate dimensionality reductions ----

seu <- ScaleData(seu) # scale the data so we can run PCA
seu <- RunPCA(seu, # run PCA
              features = rownames(seu)) # specify to use all markers for PCA (specify subset of rownames if that's instead what's desired)
ElbowPlot(seu) # make elbow plot of PCs
seu <- RunUMAP(seu, # create UMAP coordinates
               dims = 1:(ncol(seu[['pca']])), # use all PCs since we have low dimensionality
               reduction = 'pca',
               seed.use = 123, # set a seed for reproducibility
               n.neighbors = 15, # set number of neighbors between 5 and 50 to emphasize global (larger values) versus local (smaller values) structure
               spread = 1, # scale of embedded points that determines clumpiness of graph in conjunction with min.dist
               min.dist = 0.5) # minimum distance between points (larger values = greater spread of points)
# Check out https://jlmelville.github.io/uwot/abparams.html for an idea of how to toggle with spread & min.dist parameters (there's pictures!)

# Plot data ----

seu$treat <- paste(seu$WeaningGroup, seu$DaysPostWeaning_Timepoint, sep = '_') # make new metadata variable to plot
unique(seu$treat) # see what levels of variable are
seu$treat <- factor(seu$treat, # rearrange order of the new factor levels
                    levels = c('Standard_0', 'Standard_3', 'Standard_7', 'Standard_21',
                               'Late_0', 'Late_3', 'Late_7', 'Late_21'))
DimPlot(seu, # show specific treatment groups on UMAP
        group.by = 'treat', # alter this slot to group by a specific meta data variables
        shuffle = TRUE,
        reduction = 'umap')
DimPlot(seu, # plot each treatment group individually on UMAP
        split.by = 'treat', # alter this slot to split by a specific meta data variables
        shuffle = TRUE,
        reduction = 'umap',
        cols = rep('black', length(unique(seu$treat))),
        ncol = 4) & # set the number of columns to plot
        NoAxes() & NoLegend()
DimPlot(seu, # plot each sample individually on UMAP
        split.by = 'Sample', # alter this slot to split by a specific meta data variables
        shuffle = TRUE,
        reduction = 'umap',
        cols = rep('black', length(unique(seu$Sample))),
        pt.size = 0.01, # reduce point size
        ncol = 8) & NoAxes() & NoLegend()
Idents(seu) <- seu$treat
DimPlot(seu, # plot samples of only a single treatment group on UMAP
        group.by = 'Sample',
        cells = WhichCells(seu, idents = "Late_0"), # to view only cells from this treatment group
        shuffle = TRUE,
        reduction = 'umap')
FeaturePlot(seu, # plot cell markers on UMAP; this quartile scale should work for most markers
            features = c(rownames(seu)), # plot all cell markers we stored
            min.cutoff = 'q5', # set minimum gradient value to 5th quantile
            max.cutoff = 'q95') & # set maximum gradient value to 95th quantile
        scale_color_viridis(option = 'inferno') # use viridis inferno gradient scale
FeaturePlot(seu, # try this scale for markers that have very few true negative cells (like CD2)
            features = c(rownames(seu)), # plot all cell markers we stored
            min.cutoff = 'q0', # set minimum gradient value to minimum value
            max.cutoff = 'q70') & # set maximum gradient value to 95th quantile
        scale_color_viridis(option = 'inferno') # use viridis inferno gradient scale
FeaturePlot(seu, # try this scale for markers that have very few true positive cells (like CD4)
            features = c(rownames(seu)), # plot all cell markers we stored
            min.cutoff = 'q90', # set minimum gradient value to minimum value
            max.cutoff = 'q99') & # set maximum gradient value to 95th quantile
        scale_color_viridis(option = 'inferno') # use viridis inferno gradient scale

# Make another meta data slot for deviation from standard timepoint ----

seu$DOAvariation <- seu$DaysOfAge - seu$DaysOfAge_Timepoint

# Save Seurat object ----

saveRDS(seu, 
       paste(paste(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'), 'Data', sep = '/'), 'SeuratObject.rds', sep = '/'))

# Seurat plots for pub ----
FeaturePlot(seu, 
            features = 'gdTCR-BV510', 
            min.cutoff = 'q10', 
            max.cutoff = 'q90') & 
  scale_color_viridis(option = 'inferno') & NoLegend()
FeaturePlot(seu, 
            features = 'CD4-PerCPCy5.5', 
            min.cutoff = 'q90', 
            max.cutoff = 'q99') & 
  scale_color_viridis(option = 'inferno') & NoLegend()
FeaturePlot(seu, 
            features = 'CD8a-APC', 
            min.cutoff = 'q00', 
            max.cutoff = 'q85') & 
  scale_color_viridis(option = 'inferno') & NoLegend()
FeaturePlot(seu, 
            features = 'CD27-FITC', 
            min.cutoff = 'q10', 
            max.cutoff = 'q90') & 
  scale_color_viridis(option = 'inferno') & NoLegend()
FeaturePlot(seu, 
            features = 'CD16-BUV395', 
            min.cutoff = 'q05', 
            max.cutoff = 'q95') & 
  scale_color_viridis(option = 'inferno') & NoLegend()
FeaturePlot(seu, 
            features = 'MHCII-BUV496', 
            min.cutoff = 'q05', 
            max.cutoff = 'q85') & 
  scale_color_viridis(option = 'inferno') & NoLegend()
FeaturePlot(seu, 
            features = 'CD45RC-BV421', 
            min.cutoff = 'q10', 
            max.cutoff = 'q85') & 
  scale_color_viridis(option = 'inferno') & NoLegend()
FeaturePlot(seu, 
            features = 'CD2-BV711', 
            min.cutoff = 'q00', 
            max.cutoff = 'q98') & 
  scale_color_viridis(option = 'inferno') & NoLegend()
FeaturePlot(seu, 
            features = 'CD8b-PE', 
            min.cutoff = 'q10', 
            max.cutoff = 'q90') & 
  scale_color_viridis(option = 'inferno') & NoLegend()

DimPlot(seu, # plot samples of only a single treatment group on UMAP
        group.by = 'Sample',
        cells = WhichCells(seu, idents = "Late_0"), # to view only cells from this treatment group
        shuffle = TRUE,
        reduction = 'umap') & NoLegend() & NoAxes()

DimPlot(seu, # plot samples of only a single treatment group on UMAP
        group.by = 'Sample',
        cells = WhichCells(seu, idents = "Late_3"), # to view only cells from this treatment group
        shuffle = TRUE,
        reduction = 'umap') & NoLegend() & NoAxes()

DimPlot(seu, # plot samples of only a single treatment group on UMAP
        group.by = 'Sample',
        cells = WhichCells(seu, idents = "Late_7"), # to view only cells from this treatment group
        shuffle = TRUE,
        reduction = 'umap') & NoLegend() & NoAxes()

DimPlot(seu, # plot samples of only a single treatment group on UMAP
        group.by = 'Sample',
        cells = WhichCells(seu, idents = "Late_21"), # to view only cells from this treatment group
        shuffle = TRUE,
        reduction = 'umap') & NoLegend() & NoAxes()

DimPlot(seu, # plot samples of only a single treatment group on UMAP
        group.by = 'Sample',
        cells = WhichCells(seu, idents = "Standard_0"), # to view only cells from this treatment group
        shuffle = TRUE,
        reduction = 'umap') & NoLegend() & NoAxes()

DimPlot(seu, # plot samples of only a single treatment group on UMAP
        group.by = 'Sample',
        cells = WhichCells(seu, idents = "Standard_3"), # to view only cells from this treatment group
        shuffle = TRUE,
        reduction = 'umap') & NoLegend() & NoAxes()

DimPlot(seu, # plot samples of only a single treatment group on UMAP
        group.by = 'Sample',
        cells = WhichCells(seu, idents = "Standard_7"), # to view only cells from this treatment group
        shuffle = TRUE,
        reduction = 'umap') & NoLegend() & NoAxes()

DimPlot(seu, # plot samples of only a single treatment group on UMAP
        group.by = 'Sample',
        cells = WhichCells(seu, idents = "Standard_21"), # to view only cells from this treatment group
        shuffle = TRUE,
        reduction = 'umap') & NoLegend() & NoAxes()

# Create Milo object ----

set.seed(123)
milo <- as.SingleCellExperiment(seu, # convert to SingleCellExperiment object
                                assay = 'RNA') # use default RNA assay, which is where we stored our cell marker expression data
milo # see what's in our SingleCellExperiment object
milo <- Milo(milo) # convert to Milo object
milo # see what's in our Milo object; note graph names(0), meaning no nearest neighbor graph is present

rm(seu)

# Create cell neighborhoods ----

milo <- buildGraph(milo, # create nearest neighbor graph used to contruct cell neighborhoods
                   k = 500, # specify how many nearest neighbors to use. Increasing k will ensure more cells from each sample are found in a neighborhood, but increasing k too much will limit neighborhood resolution and potentially end up over-smoothing the data so significant differences aren't detected. Rule of thumb is to have between k > 3*(# of samples) to k > 5*(# of samples). Can expand further if wanting to be cautious in regards to lingering batch effects, but again could cause oversmoothing. I went ultra-conservative with high k value since only considering 1/4 of cells in any pairwise treatment comparison.
                   d = 8, # the number of dimensions to use from data; can determine by looking at scree plot if unsure... here we opt to use all dims we have (n = 8)
                   transposed = TRUE,
                   reduced.dim = "PCA") # dimensionality reduction to use
milo <- makeNhoods(milo, # create cell neighborhoods based on the kNN graph contructed created and stored in the Milo object in the preceding step
                   prop = 0.05, # proportion of cells to sample to start analysis. For scRNA-seq, recommend using 0.1 for datasets <30k. For larger datasets >100k, can decrease down to 0.05 to speed up computational time. Increasing prop will increase required computational time. May need p > 0.1 if data has rare, disconnected populations.
                   k = 500, # set k to same as for buildGraph()
                   d = 8, # set d to same as for buildGraph()
                   reduced_dims = 'PCA', # use same reduced_dims as for buildGraph()
                   refined = TRUE, # always use refined unless you use graph-based data batch correction, then consider either-or
                   refinement_scheme='graph') # use graph-based approach so bottlenecking step calcNhoodDistance() isn't required
# see Jan 27, 2022 comment on thread: https://github.com/MarioniLab/miloR/issues/108

plotNhoodSizeHist(milo) # plot neighborhood sizes

# Count cells in each neighborhood ----

milo <- countCells(milo, # count cells coming from each sample within each neighborhood; these are counts we will perform differential abundance analysis on
                   meta.data = as.data.frame(colData(milo)), 
                   sample="Sample")
head(nhoodCounts(milo)) # see counts in neighborhoods; check that we don't have lots of data sparsity here

# Visualize cell neighborhoods ----

milo <- buildNhoodGraph(milo)
plotNhoodGraph(milo, 
               layout = 'UMAP')

# Save milo object ----

dir.create(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'))
dir.create(paste(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'), 'Data', sep = '/'))
saveRDS(milo, 
        paste(paste(paste('/home/Jayne.Wiarda/FS24/HighDimAnalysis', Sys.Date(), sep = '_'), 'Data', sep = '/'), 'miloObjectK500p5d8.rds', sep = '/'))

# Create experimental design & identify contrasts ----

## Standard vs Late weaned at DPW-matched timepoints ----

### 0dpw ----
names(colData(milo)) # check out available meta data variables to use for design matrix
design <- data.frame(colData(milo))[,c("Sample", "WeaningGroup", "DaysPostWeaning_Timepoint", "Gender", "FarrowingPen", "Room", "DOAvariation", "PreweaningWeight")]
design$treat <- paste(design$WeaningGroup, design$DaysPostWeaning_Timepoint, sep = '_')
design <- design[,c('Sample', 'treat', 'Gender', 'FarrowingPen', 'Room', "DOAvariation", "PreweaningWeight")]
head(design)
design <- subset(design, treat == 'Late_0' | treat == 'Standard_0') # since we also want to check potential confounding factors (e.g. gender, litter, room) that are also nested within some weaning groups, we take the subset approach to DA testing shown here: https://github.com/MarioniLab/miloR/pull/109
design$treat <- as.factor(design$treat) 
design$Gender <- as.factor(design$Gender) 
design$FarrowingPen <- as.factor(design$FarrowingPen)
design$Room <- as.factor(design$Room)
design$PreweaningWeight <- as.numeric(design$PreweaningWeight)
design$DOAvariation <- as.numeric(design$DOAvariation)
design <- distinct(design)
rownames(design) <- design$Sample
head(design)
model <- model.matrix(~0+treat, 
                      data=design)
#head(model) # view this to see what variables we can use to set contrasts for DA testing (based on column names)
da_results <- testNhoods(milo, 
                         design = ~0+treat, 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_0 - treatLate_0')

#head(da_results)
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Late'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Standard'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is gender a confounding variable? --> not testing since no significant DA

### 3dpw ----
names(colData(milo)) # check out available meta data variables to use for design matrix
design <- data.frame(colData(milo))[,c("Sample", "WeaningGroup", "DaysPostWeaning_Timepoint", "Gender", "FarrowingPen", "Room", "DOAvariation", "PreweaningWeight")]
design$treat <- paste(design$WeaningGroup, design$DaysPostWeaning_Timepoint, sep = '_')
design <- design[,c('Sample', 'treat', 'Gender', 'FarrowingPen', 'Room', "DOAvariation", "PreweaningWeight")]
head(design)
design <- subset(design, treat == 'Late_3' | treat == 'Standard_3') # since we also want to check potential confounding factors (e.g. gender, litter, room) that are also nested within some weaning groups, we take the subset approach to DA testing shown here: https://github.com/MarioniLab/miloR/pull/109
design$treat <- as.factor(design$treat) 
design$Gender <- as.factor(design$Gender) 
design$FarrowingPen <- as.factor(design$FarrowingPen)
design$Room <- as.factor(design$Room)
design$PreweaningWeight <- as.numeric(design$PreweaningWeight)
design$DOAvariation <- as.numeric(design$DOAvariation)
design <- distinct(design)
rownames(design) <- design$Sample
head(design)
model <- model.matrix(~0+treat, 
                      data=design)
#head(model) # view this to see what variables we can use to set contrasts for DA testing (based on column names)
da_results <- testNhoods(milo, 
                         design = ~0+treat, 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_3 - treatLate_3')

#head(da_results)
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Late'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Standard'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is gender a confounding variable?
da_results <- testNhoods(milo, 
                         design = ~0+treat+Gender, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_3 - treatLate_3')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Late'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Standard'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is variation from targeted age a confounding variable?
da_results <- testNhoods(milo, 
                         design = ~0+treat+DOAvariation, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_3 - treatLate_3')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Late'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Standard'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# combined confounding variables:
da_results <- testNhoods(milo, 
                         design = ~0+treat+Gender+DOAvariation, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_3 - treatLate_3')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Late'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Standard'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

### 7dpw ----
names(colData(milo)) # check out available meta data variables to use for design matrix
design <- data.frame(colData(milo))[,c("Sample", "WeaningGroup", "DaysPostWeaning_Timepoint", "Gender", "FarrowingPen", "Room", "DOAvariation", "PreweaningWeight")]
design$treat <- paste(design$WeaningGroup, design$DaysPostWeaning_Timepoint, sep = '_')
design <- design[,c('Sample', 'treat', 'Gender', 'FarrowingPen', 'Room', "DOAvariation", "PreweaningWeight")]
head(design)
design <- subset(design, treat == 'Late_7' | treat == 'Standard_7') # since we also want to check potential confounding factors (e.g. gender, litter, room) that are also nested within some weaning groups, we take the subset approach to DA testing shown here: https://github.com/MarioniLab/miloR/pull/109
design$treat <- as.factor(design$treat) 
design$Gender <- as.factor(design$Gender) 
design$FarrowingPen <- as.factor(design$FarrowingPen)
design$Room <- as.factor(design$Room)
design$PreweaningWeight <- as.numeric(design$PreweaningWeight)
design$DOAvariation <- as.numeric(design$DOAvariation)
design <- distinct(design)
rownames(design) <- design$Sample
head(design)
model <- model.matrix(~0+treat, 
                      data=design)
#head(model) # view this to see what variables we can use to set contrasts for DA testing (based on column names)
da_results <- testNhoods(milo, 
                         design = ~0+treat, 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_7 - treatLate_7')

#head(da_results)
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Late'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Standard'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is gender a confounding variable? --> not testing since no significant DA

### 21dpw ----
names(colData(milo)) # check out available meta data variables to use for design matrix
design <- data.frame(colData(milo))[,c("Sample", "WeaningGroup", "DaysPostWeaning_Timepoint", "Gender", "FarrowingPen", "Room", "DOAvariation", "PreweaningWeight")]
design$treat <- paste(design$WeaningGroup, design$DaysPostWeaning_Timepoint, sep = '_')
design <- design[,c('Sample', 'treat', 'Gender', 'FarrowingPen', 'Room', "DOAvariation", "PreweaningWeight")]
head(design)
design <- subset(design, treat == 'Late_21' | treat == 'Standard_21') # since we also want to check potential confounding factors (e.g. gender, litter, room) that are also nested within some weaning groups, we take the subset approach to DA testing shown here: https://github.com/MarioniLab/miloR/pull/109
design$treat <- as.factor(design$treat) 
design$Gender <- as.factor(design$Gender) 
design$FarrowingPen <- as.factor(design$FarrowingPen)
design$Room <- as.factor(design$Room)
design$PreweaningWeight <- as.numeric(design$PreweaningWeight)
design$DOAvariation <- as.numeric(design$DOAvariation)
design <- distinct(design)
rownames(design) <- design$Sample
head(design)
model <- model.matrix(~0+treat, 
                      data=design)
#head(model) # view this to see what variables we can use to set contrasts for DA testing (based on column names)
da_results <- testNhoods(milo, 
                         design = ~0+treat, 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_21 - treatLate_21')

#head(da_results)
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Late'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Standard'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is gender a confounding variable? --> not testing since no significant DA

## Standard vs Late weaned at age-matched timepoints ----

### 28doa (Standard_7 v Late_0) ----
names(colData(milo)) # check out available meta data variables to use for design matrix
design <- data.frame(colData(milo))[,c("Sample", "WeaningGroup", "DaysPostWeaning_Timepoint", "Gender", "FarrowingPen", "Room", "DOAvariation", "PreweaningWeight")]
design$treat <- paste(design$WeaningGroup, design$DaysPostWeaning_Timepoint, sep = '_')
design <- design[,c('Sample', 'treat', 'Gender', 'FarrowingPen', 'Room', "DOAvariation", "PreweaningWeight")]
head(design)
design <- subset(design, treat == 'Late_0' | treat == 'Standard_7') # since we also want to check potential confounding factors (e.g. gender, litter, room) that are also nested within some weaning groups, we take the subset approach to DA testing shown here: https://github.com/MarioniLab/miloR/pull/109
design$treat <- as.factor(design$treat) 
design$Gender <- as.factor(design$Gender) 
design$FarrowingPen <- as.factor(design$FarrowingPen)
design$Room <- as.factor(design$Room)
design$PreweaningWeight <- as.numeric(design$PreweaningWeight)
design$DOAvariation <- as.numeric(design$DOAvariation)
design <- distinct(design)
rownames(design) <- design$Sample
head(design)
model <- model.matrix(~0+treat, 
                      data=design)
#head(model) # view this to see what variables we can use to set contrasts for DA testing (based on column names)
da_results <- testNhoods(milo, 
                         design = ~0+treat, 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_7 - treatLate_0')

#head(da_results)
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Late'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Standard'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is gender a confounding variable? 
da_results <- testNhoods(milo, 
                         design = ~0+treat+Gender, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_7 - treatLate_0')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Late'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Standard'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is variation from targeted age a confounding variable?
da_results <- testNhoods(milo, 
                         design = ~0+treat+DOAvariation, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_7 - treatLate_0')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Late'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Standard'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# combined confounding variables:
da_results <- testNhoods(milo, 
                         design = ~0+treat+Gender+DOAvariation, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_7 - treatLate_0')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Late'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Standard'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

## Comparing sequential timepoints within weaning group ----

### Late 0v3dpw ----
names(colData(milo)) # check out available meta data variables to use for design matrix
design <- data.frame(colData(milo))[,c("Sample", "WeaningGroup", "DaysPostWeaning_Timepoint", "Gender", "FarrowingPen", "Room", "DOAvariation", "PreweaningWeight")]
design$treat <- paste(design$WeaningGroup, design$DaysPostWeaning_Timepoint, sep = '_')
design <- design[,c('Sample', 'treat', 'Gender', 'FarrowingPen', 'Room', "DOAvariation", "PreweaningWeight")]
head(design)
design <- subset(design, treat == 'Late_0' | treat == 'Late_3') # since we also want to check potential confounding factors (e.g. gender, litter, room) that are also nested within some weaning groups, we take the subset approach to DA testing shown here: https://github.com/MarioniLab/miloR/pull/109
design$treat <- as.factor(design$treat) 
design$Gender <- as.factor(design$Gender) 
design$FarrowingPen <- as.factor(design$FarrowingPen)
design$Room <- as.factor(design$Room)
design$PreweaningWeight <- as.numeric(design$PreweaningWeight)
design$DOAvariation <- as.numeric(design$DOAvariation)
design <- distinct(design)
rownames(design) <- design$Sample
head(design)
model <- model.matrix(~0+treat, 
                      data=design)
#head(model) # view this to see what variables we can use to set contrasts for DA testing (based on column names)
da_results <- testNhoods(milo, 
                         design = ~0+treat, 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatLate_0 - treatLate_3')

#head(da_results)
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is gender a confounding variable? 
da_results <- testNhoods(milo, 
                         design = ~0+treat+Gender, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatLate_0 - treatLate_3')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is farrowing pen a confounding variable? 
da_results <- testNhoods(milo, 
                         design = ~0+treat+FarrowingPen, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatLate_0 - treatLate_3')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is preweaning weight a confounding variable?
da_results <- testNhoods(milo, 
                         design = ~0+treat+PreweaningWeight, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatLate_0 - treatLate_3')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Late'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Standard'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is room a confounding variable? --> can't test since no room assignments at 0dpw

# combined confounding variables:
da_results <- testNhoods(milo, 
                         design = ~0+treat+Gender+FarrowingPen+PreweaningWeight, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatLate_0 - treatLate_3')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

### Late 3v7dpw ----
names(colData(milo)) # check out available meta data variables to use for design matrix
design <- data.frame(colData(milo))[,c("Sample", "WeaningGroup", "DaysPostWeaning_Timepoint", "Gender", "FarrowingPen", "Room", "DOAvariation", "PreweaningWeight")]
design$treat <- paste(design$WeaningGroup, design$DaysPostWeaning_Timepoint, sep = '_')
design <- design[,c('Sample', 'treat', 'Gender', 'FarrowingPen', 'Room', "DOAvariation", "PreweaningWeight")]
head(design)
design <- subset(design, treat == 'Late_3' | treat == 'Late_7') # since we also want to check potential confounding factors (e.g. gender, litter, room) that are also nested within some weaning groups, we take the subset approach to DA testing shown here: https://github.com/MarioniLab/miloR/pull/109
design$treat <- as.factor(design$treat) 
design$Gender <- as.factor(design$Gender) 
design$FarrowingPen <- as.factor(design$FarrowingPen)
design$Room <- as.factor(design$Room)
design$PreweaningWeight <- as.numeric(design$PreweaningWeight)
design$DOAvariation <- as.numeric(design$DOAvariation)
design <- distinct(design)
rownames(design) <- design$Sample
head(design)
model <- model.matrix(~0+treat, 
                      data=design)
#head(model) # view this to see what variables we can use to set contrasts for DA testing (based on column names)
da_results <- testNhoods(milo, 
                         design = ~0+treat, 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatLate_3 - treatLate_7')

#head(da_results)
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# not testing for confounding variables since no significant DA

### Late 7v21dpw ----
names(colData(milo)) # check out available meta data variables to use for design matrix
design <- data.frame(colData(milo))[,c("Sample", "WeaningGroup", "DaysPostWeaning_Timepoint", "Gender", "FarrowingPen", "Room", "DOAvariation", "PreweaningWeight")]
design$treat <- paste(design$WeaningGroup, design$DaysPostWeaning_Timepoint, sep = '_')
design <- design[,c('Sample', 'treat', 'Gender', 'FarrowingPen', 'Room', "DOAvariation", "PreweaningWeight")]
head(design)
design <- subset(design, treat == 'Late_7' | treat == 'Late_21') # since we also want to check potential confounding factors (e.g. gender, litter, room) that are also nested within some weaning groups, we take the subset approach to DA testing shown here: https://github.com/MarioniLab/miloR/pull/109
design$treat <- as.factor(design$treat) 
design$Gender <- as.factor(design$Gender) 
design$FarrowingPen <- as.factor(design$FarrowingPen)
design$Room <- as.factor(design$Room)
design$PreweaningWeight <- as.numeric(design$PreweaningWeight)
design$DOAvariation <- as.numeric(design$DOAvariation)
design <- distinct(design)
rownames(design) <- design$Sample
head(design)
model <- model.matrix(~0+treat, 
                      data=design)
#head(model) # view this to see what variables we can use to set contrasts for DA testing (based on column names)
da_results <- testNhoods(milo, 
                         design = ~0+treat, 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatLate_7 - treatLate_21')

#head(da_results)
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# not testing for confounding variables since no significant DA

### Standard 0v3dpw ----
names(colData(milo)) # check out available meta data variables to use for design matrix
design <- data.frame(colData(milo))[,c("Sample", "WeaningGroup", "DaysPostWeaning_Timepoint", "Gender", "FarrowingPen", "Room", "DOAvariation", "PreweaningWeight")]
design$treat <- paste(design$WeaningGroup, design$DaysPostWeaning_Timepoint, sep = '_')
design <- design[,c('Sample', 'treat', 'Gender', 'FarrowingPen', 'Room', "DOAvariation", "PreweaningWeight")]
head(design)
design <- subset(design, treat == 'Standard_0' | treat == 'Standard_3') # since we also want to check potential confounding factors (e.g. gender, litter, room) that are also nested within some weaning groups, we take the subset approach to DA testing shown here: https://github.com/MarioniLab/miloR/pull/109
design$treat <- as.factor(design$treat) 
design$Gender <- as.factor(design$Gender) 
design$FarrowingPen <- as.factor(design$FarrowingPen)
design$Room <- as.factor(design$Room)
design$PreweaningWeight <- as.numeric(design$PreweaningWeight)
design$DOAvariation <- as.numeric(design$DOAvariation)
design <- distinct(design)
rownames(design) <- design$Sample
head(design)
model <- model.matrix(~0+treat, 
                      data=design)
#head(model) # view this to see what variables we can use to set contrasts for DA testing (based on column names)
da_results <- testNhoods(milo, 
                         design = ~0+treat, 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_0 - treatStandard_3')

#head(da_results)
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is gender a confounding variable? 
da_results <- testNhoods(milo, 
                         design = ~0+treat+Gender, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_0 - treatStandard_3')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is farrowing pen a confounding variable? 
da_results <- testNhoods(milo, 
                         design = ~0+treat+FarrowingPen, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_0 - treatStandard_3')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is preweaning weight a confounding variable?
da_results <- testNhoods(milo, 
                         design = ~0+treat+PreweaningWeight, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_0 - treatStandard_3')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Late'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Standard'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is room a confounding variable? --> can't test since no room assignments at 0dpw

# combined confounding variables:
da_results <- testNhoods(milo, 
                         design = ~0+treat+Gender+FarrowingPen+PreweaningWeight, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_0 - treatStandard_3')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

### Standard 3v7dpw ----
names(colData(milo)) # check out available meta data variables to use for design matrix
design <- data.frame(colData(milo))[,c("Sample", "WeaningGroup", "DaysPostWeaning_Timepoint", "Gender", "FarrowingPen", "Room", "DOAvariation", "PreweaningWeight")]
design$treat <- paste(design$WeaningGroup, design$DaysPostWeaning_Timepoint, sep = '_')
design <- design[,c('Sample', 'treat', 'Gender', 'FarrowingPen', 'Room', "DOAvariation", "PreweaningWeight")]
head(design)
design <- subset(design, treat == 'Standard_3' | treat == 'Standard_7') # since we also want to check potential confounding factors (e.g. gender, litter, room) that are also nested within some weaning groups, we take the subset approach to DA testing shown here: https://github.com/MarioniLab/miloR/pull/109
design$treat <- as.factor(design$treat) 
design$Gender <- as.factor(design$Gender) 
design$FarrowingPen <- as.factor(design$FarrowingPen)
design$Room <- as.factor(design$Room)
design$PreweaningWeight <- as.numeric(design$PreweaningWeight)
design$DOAvariation <- as.numeric(design$DOAvariation)
design <- distinct(design)
rownames(design) <- design$Sample
head(design)
model <- model.matrix(~0+treat, 
                      data=design)
#head(model) # view this to see what variables we can use to set contrasts for DA testing (based on column names)
da_results <- testNhoods(milo, 
                         design = ~0+treat, 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_3 - treatStandard_7')

#head(da_results)
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is gender a confounding variable? 
da_results <- testNhoods(milo, 
                         design = ~0+treat+Gender, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_3 - treatStandard_7')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is farrowing pen a confounding variable? 
da_results <- testNhoods(milo, 
                         design = ~0+treat+FarrowingPen, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_3 - treatStandard_7')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is preweaning weight a confounding variable?
da_results <- testNhoods(milo, 
                         design = ~0+treat+PreweaningWeight, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_3 - treatStandard_7')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Late'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Standard'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is room a confounding variable? 
da_results <- testNhoods(milo, 
                         design = ~0+treat+Room, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_3 - treatStandard_7')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Late'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Standard'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# combined confounding variables:
da_results <- testNhoods(milo, 
                         design = ~0+treat+Gender+FarrowingPen+PreweaningWeight+Room, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_3 - treatStandard_7')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

### Standard 7v21dpw ----
names(colData(milo)) # check out available meta data variables to use for design matrix
design <- data.frame(colData(milo))[,c("Sample", "WeaningGroup", "DaysPostWeaning_Timepoint", "Gender", "FarrowingPen", "Room", "DOAvariation", "PreweaningWeight")]
design$treat <- paste(design$WeaningGroup, design$DaysPostWeaning_Timepoint, sep = '_')
design <- design[,c('Sample', 'treat', 'Gender', 'FarrowingPen', 'Room', "DOAvariation", "PreweaningWeight")]
head(design)
design <- subset(design, treat == 'Standard_7' | treat == 'Standard_21') # since we also want to check potential confounding factors (e.g. gender, litter, room) that are also nested within some weaning groups, we take the subset approach to DA testing shown here: https://github.com/MarioniLab/miloR/pull/109
design$treat <- as.factor(design$treat) 
design$Gender <- as.factor(design$Gender) 
design$FarrowingPen <- as.factor(design$FarrowingPen)
design$Room <- as.factor(design$Room)
design$PreweaningWeight <- as.numeric(design$PreweaningWeight)
design$DOAvariation <- as.numeric(design$DOAvariation)
design <- distinct(design)
rownames(design) <- design$Sample
head(design)
model <- model.matrix(~0+treat, 
                      data=design)
#head(model) # view this to see what variables we can use to set contrasts for DA testing (based on column names)
da_results <- testNhoods(milo, 
                         design = ~0+treat, 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_7 - treatStandard_21')

#head(da_results)
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is gender a confounding variable? 
da_results <- testNhoods(milo, 
                         design = ~0+treat+Gender, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_7 - treatStandard_21')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is farrowing pen a confounding variable? 
da_results <- testNhoods(milo, 
                         design = ~0+treat+FarrowingPen, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_7 - treatStandard_21')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is preweaning weight a confounding variable?
da_results <- testNhoods(milo, 
                         design = ~0+treat+PreweaningWeight, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_7 - treatStandard_21')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Late'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Standard'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# is room a confounding variable? 
da_results <- testNhoods(milo, 
                         design = ~0+treat+Room, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_7 - treatStandard_21')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Late'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Standard'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)

# combined confounding variables:
da_results <- testNhoods(milo, 
                         design = ~0+treat+Gender+FarrowingPen+PreweaningWeight+Room, # don't include farrowing pen or gender since these variables are nested within standard vs late weaning groups 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_7 - treatStandard_21')
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)
ggplot(da_results, aes(SpatialFDR)) + geom_histogram(bins=50)
plotNhoodGraphDA(milo, 
                 da_results, 
                 layout="UMAP",
                 alpha=0.05) # specify alpha level to use (same as used for calling da_results above)








## Testing ----
names(colData(milo)) # check out available meta data variables to use for design matrix
design <- data.frame(colData(milo))[,c("Sample", "WeaningGroup", "DaysPostWeaning_Timepoint", "Gender", "FarrowingPen", "Room", "DOAvariation", "PreweaningWeight")]
design$treat <- paste(design$WeaningGroup, design$DaysPostWeaning_Timepoint, sep = '_')
design <- design[,c('Sample', 'treat', 'Gender', 'FarrowingPen', 'Room', "DOAvariation", "PreweaningWeight")]
head(design)
design <- subset(design, treat == 'Late_3' | treat == 'Standard_3') # since we also want to check potential confounding factors (e.g. gender, litter, room) that are also nested within some weaning groups, we take the subset approach to DA testing shown here: https://github.com/MarioniLab/miloR/pull/109
design$treat <- as.factor(design$treat) 
design$Gender <- as.factor(design$Gender) 
design$FarrowingPen <- as.factor(design$FarrowingPen)
design$Room <- as.factor(design$Room)
design$PreweaningWeight <- as.numeric(design$PreweaningWeight)
design$DOAvariation <- as.numeric(design$DOAvariation)
design <- distinct(design)
rownames(design) <- design$Sample
head(design)
model <- model.matrix(~0+treat, 
                      data=design)
#head(model) # view this to see what variables we can use to set contrasts for DA testing (based on column names)
da_results <- testNhoods(milo, 
                         design = ~0+treat, 
                         design.df = design,
                         fdr.weighting = 'graph-overlap', # make sure to use graph-overlap since we used graph-based approach for makeNhoods()
                         model.contrasts = 'treatStandard_3 - treatLate_3')

#head(da_results)
da_results$Sig[da_results$SpatialFDR >= 0.05] <- 'NS'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Late'
da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Standard'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC > 0] <- 'Earlier'
#da_results$Sig[da_results$SpatialFDR < 0.05 & da_results$logFC < 0] <- 'Later'
table(da_results$Sig)

da_results <- groupNhoods(milo, da_results, da.fdr = 0.05, max.lfc.delta = 3, overlap = 100)
plotNhoodGroups(milo, da_results, layout="UMAP") 

sessionInfo()
R version 4.2.2 Patched (2022-11-10 r83330)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 22.04.3 LTS

Matrix products: default
BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3
LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.20.so

locale:
  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
[5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8    LC_PAPER=en_US.UTF-8       LC_NAME=C                 
[9] LC_ADDRESS=C               LC_TELEPHONE=C             LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
  [1] stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
  [1] data.table_1.14.8           viridis_0.6.4               viridisLite_0.4.2          
[4] dplyr_1.1.3                 miloR_1.4.0                 edgeR_3.38.4               
[7] limma_3.52.4                SingleCellExperiment_1.18.1 SummarizedExperiment_1.26.1
[10] Biobase_2.56.0              GenomicRanges_1.48.0        GenomeInfoDb_1.32.4        
[13] IRanges_2.30.1              S4Vectors_0.34.0            BiocGenerics_0.42.0        
[16] MatrixGenerics_1.8.1        matrixStats_1.0.0           cowplot_1.1.1              
[19] Rtsne_0.16                  Spectre_1.0.0               SeuratDisk_0.0.0.9020      
[22] ggplot2_3.4.3               SeuratObject_4.1.3          Seurat_4.3.0.1             

loaded via a namespace (and not attached):
  [1] circlize_0.4.15        plyr_1.8.8             igraph_1.5.1           lazyeval_0.2.2        
[5] sp_2.0-0               splines_4.2.2          BiocParallel_1.30.4    listenv_0.9.0         
[9] scattermore_1.2        digest_0.6.33          foreach_1.5.2          htmltools_0.5.6       
[13] fansi_1.0.4            magrittr_2.0.3         ScaledMatrix_1.4.1     tensor_1.5            
[17] cluster_2.1.4          doParallel_1.0.17      ROCR_1.0-11            graphlayouts_1.0.0    
[21] ComplexHeatmap_2.12.1  globals_0.16.2         spatstat.sparse_3.0-2  colorspace_2.1-0      
[25] ggrepel_0.9.3          xfun_0.40              crayon_1.5.2           RCurl_1.98-1.12       
[29] jsonlite_1.8.7         progressr_0.14.0       spatstat.data_3.0-1    survival_3.5-7        
[33] zoo_1.8-12             iterators_1.0.14       glue_1.6.2             polyclip_1.10-4       
[37] gtable_0.3.4           zlibbioc_1.42.0        XVector_0.36.0         leiden_0.4.3          
[41] DelayedArray_0.22.0    GetoptLong_1.0.5       BiocSingular_1.12.0    future.apply_1.11.0   
[45] shape_1.4.6            abind_1.4-5            scales_1.2.1           DBI_1.1.3             
[49] spatstat.random_3.1-6  miniUI_0.1.1.1         Rcpp_1.0.11            xtable_1.8-4          
[53] clue_0.3-64            reticulate_1.32.0      rsvd_1.0.5             bit_4.0.5             
[57] htmlwidgets_1.6.2      httr_1.4.7             RColorBrewer_1.1-3     ellipsis_0.3.2        
[61] ica_1.0-3              farver_2.1.1           pkgconfig_2.0.3        uwot_0.1.16           
[65] deldir_1.0-9           locfit_1.5-9.8         utf8_1.2.3             tidyselect_1.2.0      
[69] rlang_1.1.1            reshape2_1.4.4         later_1.3.1            munsell_0.5.0         
[73] tools_4.2.2            cli_3.6.1              generics_0.1.3         ggridges_0.5.4        
[77] evaluate_0.21          stringr_1.5.0          fastmap_1.1.1          yaml_2.3.7            
[81] goftest_1.2-3          knitr_1.44             bit64_4.0.5            tidygraph_1.2.3       
[85] fitdistrplus_1.1-11    purrr_1.0.2            RANN_2.6.1             ggraph_2.1.0          
[89] pbapply_1.7-2          future_1.33.0          nlme_3.1-163           mime_0.12             
[93] hdf5r_1.3.8            compiler_4.2.2         rstudioapi_0.15.0      beeswarm_0.4.0        
[97] plotly_4.10.2          png_0.1-8              spatstat.utils_3.0-3   tweenr_2.0.2          
[101] tibble_3.2.1           stringi_1.7.12         lattice_0.21-8         Matrix_1.6-1          
[105] vctrs_0.6.3            pillar_1.9.0           lifecycle_1.0.3        spatstat.geom_3.2-5   
[109] lmtest_0.9-40          GlobalOptions_0.1.2    BiocNeighbors_1.14.0   RcppAnnoy_0.0.21      
[113] bitops_1.0-7           irlba_2.3.5.1          httpuv_1.6.11          patchwork_1.1.3       
[117] R6_2.5.1               promises_1.2.1         KernSmooth_2.23-22     gridExtra_2.3         
[121] vipor_0.4.5            parallelly_1.36.0      codetools_0.2-19       gtools_3.9.4          
[125] MASS_7.3-60            rjson_0.2.21           withr_2.5.0            sctransform_0.3.5     
[129] GenomeInfoDbData_1.2.8 parallel_4.2.2         beachmat_2.12.0        grid_4.2.2            
[133] tidyr_1.3.0            rmarkdown_2.24         ggforce_0.4.1          spatstat.explore_3.2-3
[137] shiny_1.7.5            ggbeeswarm_0.7.2    
