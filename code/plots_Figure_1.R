# Make composite plots of pi and Tajima's D obtained with popoolation v.1
# Code used to produce Figure 1 and Supplementary Figure X in Chen et al. 2021
# By: Angela Fuentes-Pardo, email: apfuentesp@gmail.com
# Uppsala University
# Date: 2021-04-25

# Clean environment space
rm(list=ls())

# Set environment variables
working_dir <- '~/path/pi_TajD_FST'
dAF_file <- '~/path/poolAF.txt'
#plot_type <- 'all'
plot_type <- 'zoom'
region <- 'chr15:8-10Mbp'
#region <- 'chr15:6-12Mbp'

# Set working dir
setwd(working_dir)

# Load dAF data
library(data.table)
dAF_df <- fread(dAF_file, data.table = FALSE, header = TRUE, sep = '\t', stringsAsFactors = FALSE)
head(dAF_df)


# --------------------------------------------------------------------------------------
# Make boxplots of pi and Tajima's D within and outside TSHR
# --------------------------------------------------------------------------------------

# TSHR region coordinates in the paper: 8.85 - 8.95 Mb

# Set min coverage threshold
cov <- 0.5

## Load pi data -----------
files1 <- list.files("./results/subsampling/", pattern = "\\.pi.txt$")
files1

# Initialize a list to store the plots of the current file
df1_list <- list()

for(x in files1) {  # loop over files
  #x="DalGeB_Atlantic_Autumn.chr15.wd.10000.st.2000.idl.scov.pi.txt"
  
  # Save sample name
  sampleName <- gsub('.chr.+$','',x)
  
  # Load pi data
  pi_df <- fread(paste0('./results/subsampling/',x), header = F, data.table = F)
  colnames(pi_df) <- c('contig','bp','SNPcount','coverage','pi')
  head(pi_df)
  #hist(pi_df$coverage)
  pi_df$pi <- gsub("na",NA,pi_df$pi)
  pi_df <- pi_df[!is.na(pi_df$pi), ]
  pi_df$pi <- as.numeric(pi_df$pi)
  #hist(pi_df$pi)
  colnames(pi_df) <- gsub('pi',sampleName,colnames(pi_df))  # assign the sample name to the pi column
  
  # Filter windows by minimum  depth of coverage %, and subset to the windows in the TSHR region
  df1_list[[sampleName]] <- pi_df[pi_df$coverage >= cov, c('contig','bp',sampleName)]
}

## Load Tajima's D data -----------
files2 <- list.files("./results/subsampling/", pattern = "\\.TajD.txt$")
files2

# Initialize a list to store the plots of the current file
df2_list <- list()

for(x in files2) {  # loop over files
  #x="DalGeB_Atlantic_Autumn.chr15.wd.10000.st.2000.idl.scov.TajD.txt"
  
  # Save sample name
  sampleName <- gsub('.chr.+$','',x)
  
  # Load TajD data
  TajD_df <- fread(paste0('./results/subsampling/',x), header = F, data.table = F)
  colnames(TajD_df) <- c('contig','bp','SNPcount','coverage','TajD')
  head(TajD_df)
  #hist(TajD_df$coverage)
  TajD_df$TajD <- gsub("na",NA,TajD_df$TajD)
  TajD_df <- TajD_df[!is.na(TajD_df$TajD), ]
  TajD_df$TajD <- as.numeric(TajD_df$TajD)
  #hist(TajD_df$TajD)
  colnames(TajD_df) <- gsub('TajD',sampleName,colnames(TajD_df))  # assign the sample name to the TajD column
  
  # Filter windows by minimum  depth of coverage %, and subset to the windows in the TSHR region
  df2_list[[sampleName]] <- TajD_df[TajD_df$coverage >= cov, c('contig','bp',sampleName)]
}

# Merge the dataframes of each of the stats
df1 <- df1_list %>% purrr::reduce(left_join, by = c('contig','bp'))
df2 <- df2_list %>% purrr::reduce(left_join, by = c('contig','bp'))
head(df1)
head(df2)
str(df2)

# Make two dfs, one with the values within TSHR and other with the values outside
df1_tshr <- df1[df1$bp >= 8850000 & df1$bp <= 8950000, ]
df1_out <- df1[df1$bp < 8850000 | df1$bp > 8950000, ]
df2_tshr <- df2[df2$bp >= 8850000 & df2$bp <= 8950000, ]
df2_out <- df2[df2$bp < 8850000 | df2$bp > 8950000, ]

nrow(df1_tshr)
#[1] 41
nrow(df1_out)
#[1] 12497
nrow(df2_tshr)
#[1] 41
nrow(df2_out)
#[1] 12497

# Store the dfs in a list for post-processing
df_list <- list(df1_tshr, df1_out, df2_tshr, df2_out)
names(df_list) <- c('df1_tshr', 'df1_out', 'df2_tshr', 'df2_out')
length(df_list)

# Create empty list to store post-processed dfs
dff_list <- list()

for(dfName in names(df_list)){
  #dfName="df1_tshr"
  dff <- df_list[[dfName]]
  head(dff)
  
  # Merge dfs by contig and bp in common
  dff$wd <- paste0(dff$contig,'-', dff$bp)
  dff$contig <- NULL
  dff$bp <- NULL
  head(dff)
  #tmp <- as.data.frame(t(dff),stringsAsFactors = F)
  #head(tmp)
  
  # Convert the long to the wide format
  library(reshape)
  dff_melt <- reshape2::melt(dff)
  dff_melt <- tidyr::separate(data = dff_melt, col = wd, into = c("contig", "bp"), sep = "-")  # split wd into contig and bp
  head(dff_melt)
  
  # Assign the color you want for each line
  dff_melt$season <- NA
  dff_melt[grepl('Spring',dff_melt$variable), 'season'] <- 'spring'
  dff_melt[grepl('Autumn',dff_melt$variable), 'season'] <- 'autumn'
  dff_melt$bp = as.integer(dff_melt$bp)
  dff_melt <- dff_melt[order(dff_melt$bp, decreasing = F), ]
  if(grepl('df1',dfName)) { dff_melt$stat_type <- "pi" }
  if(grepl('df2',dfName)) { dff_melt$stat_type <- "TajD" }
  if(grepl('tshr',dfName)) { dff_melt$category <- "TSHR" }
  if(grepl('out',dfName)) { dff_melt$category <- "background" }
  head(dff_melt)
  
  dff_list[[dfName]] <- dff_melt
}

# Make the boxplot - with significance values (Wilcoxon test, non-parametric)  -------------

library(ggpubr)

## Make the plot for pi --------
d1_in <- dff_list[["df1_tshr"]]
d1_out <- dff_list[["df1_out"]]
head(d1_in)
head(d1_out)
dd1 <- rbind(d1_in, d1_out, stringsAsFactors = F)
head(dd1)

# Specify the comparisons you want
my_comparisons <- list(c("TSHR", "background"))

pp1 <- ggboxplot(dd1, x = "category", y = "value",
                 facet.by = "season", ylim = c(0, 0.04), 
                 color = "category", palette = c("#a6611a","#018571")) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", group.by = "season"#,
                     #label = "p.signif"
  ) + # add pairwise comparisons p-value
  #stat_compare_means(label.y = 0.6) + # add global p-value, label.y denotes the value in the y-axis where to print
  theme_bw() + ylab('pi') +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    text = element_text(size=18),
    #axis.title.y=element_blank()
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
  )
pp1

## Make the plot for Tajima's D --------
d2_in <- dff_list[["df2_tshr"]]
d2_out <- dff_list[["df2_out"]]
head(d2_in)
head(d2_out)
dd2 <- rbind(d2_in, d2_out, stringsAsFactors = F)
head(dd2)

# Specify the comparisons you want
my_comparisons <- list(c("TSHR", "background"))

pp2 <- ggboxplot(dd2, x = "category", y = "value",
                 facet.by = "season", ylim = c(-2, 1.7), 
                 color = "category", palette = c("#a6611a","#018571")) + 
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", group.by = "season"#,
                     #label = "p.signif"
  ) + # add pairwise comparisons p-value
  #stat_compare_means(label.y = 0.6) + # add global p-value, label.y denotes the value in the y-axis where to print
  theme_bw() + ylab('Tajima\'s D') +
  theme(#axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
    text = element_text(size=18),
    #axis.title.y=element_blank()
    axis.title.x=element_blank(),
    axis.text.x=element_blank(),
    axis.ticks.x=element_blank()
  )
pp2


## Make a combined plot ------------
library(patchwork)
#joint_plot <- pp1 / pp2 + plot_layout(guides = "collect")
joint_plot <- pp1 / pp2 / pp3 + plot_layout(guides = "collect")
joint_plot

pdf('boxplots_pi_TajD_background_vs_TSHR_Herring_pools_pi_TajD.pdf')
#pdf('boxplots_pi_TajD_Hp.background_vs_TSHR_Herring_pools.pdf')
#pdf('boxplots_pi_TajD_Hp.background_vs_TSHR_Herring_pools.pdf',5,7)
print(joint_plot)
dev.off()

# Save plots
saveRDS(joint_plot, file='boxplots_pi_TajD_Hp.background_vs_TSHR_Herring_pools.RDS')

# --------------------------------------------------------------------------------------
# Make a combined close-up plot for the TSHR region
# --------------------------------------------------------------------------------------

# Set min coverage threshold
cov <- 0.5

if(plot_type == 'zoom' & !is.null(region)){
  
  # --------------------------------------------------------------------------------------
  # Process pi data
  # --------------------------------------------------------------------------------------
  
  # Load data
  files1 <- list.files("./results/subsampling/", pattern = "\\.pi.txt$")
  #order_files <- c(8,1,10,3,5,7,2,9,4,6)
  #files <- files[order_files]  # sort
  #files_spring <- files[grep('Spring',files)]
  #files_autumn <- files[grep('Autumn',files)]
  files1
  #files_spring
  #files_autumn
  
  # Initialize a list to store the plots of the current file
  df1_list <- list()
  
  for(x in files1) {  # loop over files
    #x="DalGeB_Atlantic_Autumn.chr15.wd.10000.st.2000.idl.scov.pi.txt"
    #x="DalInB_Atlantic_Spring.chr15.wd.10000.st.2000.idl.scov.pi.txt"
    
    # Save sample name
    sampleName <- gsub('.chr.+$','',x)
    
    # Load pi data
    pi_df <- fread(paste0('./results/subsampling/',x), header = F, data.table = F)
    colnames(pi_df) <- c('contig','bp','SNPcount','coverage','pi')
    head(pi_df)
    #hist(pi_df$coverage)
    pi_df$pi <- gsub("na",NA,pi_df$pi)
    pi_df <- pi_df[!is.na(pi_df$pi), ]
    pi_df$pi <- as.numeric(pi_df$pi)
    #hist(pi_df$pi)
    colnames(pi_df) <- gsub('pi',sampleName,colnames(pi_df))  # assign the sample name to the pi column
    
    # Create an output file name
    outputname <- x
    #outputname <- gsub('HGS36_MilfordHavenApr_Atlantic_Spring.chr15.','',x)
    outputname <- gsub('.idl.scov.pi.txt','',outputname)
    if(region == 'chr15:8-10Mbp'){ outputname <- paste0(outputname,'.chr15:8-10Mbp') }
    if(region == 'chr15:6-12Mbp'){ outputname <- paste0(outputname,'chr15:6-12Mbp') }
    outputname
    
    # Filter windows by minimum  depth of coverage %, and subset to the windows in the TSHR region
    if(region == 'chr15:8-10Mbp'){ df1_list[[sampleName]] <- pi_df[pi_df$coverage >= cov & pi_df$bp >= 8000000 & pi_df$bp <= 10000000, c('contig','bp',sampleName)] }
    if(region == 'chr15:6-12Mbp'){ df1_list[[sampleName]] <- pi_df[pi_df$coverage >= cov & pi_df$bp >= 6000000 & pi_df$bp <= 12000000, c('contig','bp',sampleName)] }
    
  }
  #head(df1_list)
  
  # Merge dfs by contig and bp in common
  library(tidyverse)
  df1 <- df1_list %>% purrr::reduce(left_join, by = c('contig','bp'))
  df1$wd <- paste0(df1$contig,'-', df1$bp)
  df1$contig <- NULL
  df1$bp <- NULL
  head(df1)
  #tmp <- as.data.frame(t(df1),stringsAsFactors = F)
  #head(tmp)
  
  # Convert the long to the wide format
  library(reshape)
  df1_melt <- reshape2::melt(df1)
  df1_melt <- tidyr::separate(data = df1_melt, col = wd, into = c("contig", "bp"), sep = "-")  # split wd into contig and bp
  head(df1_melt)
  
  # Assign the color you want for each line
  df1_melt$season <- NA
  df1_melt[grepl('Spring',df1_melt$variable), 'season'] <- 'spring'
  df1_melt[grepl('Autumn',df1_melt$variable), 'season'] <- 'autumn'
  df1_melt$bp = as.integer(df1_melt$bp)
  df1_melt <- df1_melt[order(df1_melt$bp, decreasing = F), ]
  head(df1_melt)
  
  # Make the plot
  p1 <- ggplot(data=df1_melt, aes(x=bp, y=value, group=variable, colour=season)) + 
    geom_line(size=0.3) +
    #geom_line(size=0.3) +  # zoom
    theme_classic() + 
    #scale_linetype_manual(values=c(0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7)) +
    scale_color_manual(values=c('#377eb8','#e41a1c')) +
    #ggtitle(paste0('Comparison pi among spring and autumn pools', ', filter %cov >= ',cov)) +
    #scale_x_continuous(limits = c(8700000,9010000)) +  # zoom
    #scale_y_continuous(limits = c(0.000,0.005)) +  # zoom
    #scale_x_continuous(limits = c(8780000,9050000)) +  # zoom
    #scale_x_continuous(breaks=as.character(seq(min(df_melt$bp),max(df_melt$bp), 100000)))
    #xlab('position (bp)') + 
    ylab('pi') +
    theme(text = element_text(size=18),
          #axis.title.y=element_blank()
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()
    )
  #p1
  
  #pdf('zoom_pi_drop.pdf')  # zoom
  #print(p1)
  #dev.off()
  
  # --------------------------------------------------------------------------------------
  # Process Tajima's D data
  # --------------------------------------------------------------------------------------
  
  # Load data
  files2 <- list.files("./results/subsampling/", pattern = "\\.TajD.txt$")
  #order_files <- c(8,1,10,3,5,7,2,9,4,6)
  #files <- files[order_files]  # sort
  #files_spring <- files[grep('Spring',files)]
  #files_autumn <- files[grep('Autumn',files)]
  files2
  #files_spring
  #files_autumn
  
  # Initialize a list to store the plots of the current file
  df2_list <- list()
  
  for(x in files2) {  # loop over files
    #x="DalGeB_Atlantic_Autumn.chr15.wd.10000.st.2000.idl.scov.TajD.txt"
    #x="DalInB_Atlantic_Spring.chr15.wd.10000.st.2000.idl.scov.pi.txt"
    
    # Save sample name
    sampleName <- gsub('.chr.+$','',x)
    
    # Load pi data
    TajD_df <- fread(paste0('./results/subsampling/',x), header = F, data.table = F)
    colnames(TajD_df) <- c('contig','bp','SNPcount','coverage','TajD')
    head(TajD_df)
    #hist(TajD_df$coverage)
    TajD_df$TajD <- gsub("na",NA,TajD_df$TajD)
    TajD_df <- TajD_df[!is.na(TajD_df$TajD), ]
    TajD_df$TajD <- as.numeric(TajD_df$TajD)
    #hist(TajD_df$TajD)
    colnames(TajD_df) <- gsub('TajD',sampleName,colnames(TajD_df))  # assign the sample name to the pi column
    
    # Create an output file name
    outputname <- x
    #outputname <- gsub('HGS36_MilfordHavenApr_Atlantic_Spring.chr15.','',x)
    outputname <- gsub('.idl.scov.TajD.txt','',outputname)
    if(region == 'chr15:8-10Mbp'){ outputname <- paste0(outputname,'.chr15:8-10Mbp') }
    if(region == 'chr15:6-12Mbp'){ outputname <- paste0(outputname,'chr15:6-12Mbp') }
    outputname
    
    # Filter windows by minimum  depth of coverage %, and subset to the windows in the TSHR region
    if(region == 'chr15:8-10Mbp'){ df2_list[[sampleName]] <- TajD_df[TajD_df$coverage >= cov & TajD_df$bp >= 8000000 & TajD_df$bp <= 10000000, c('contig','bp',sampleName)] }
    if(region == 'chr15:6-12Mbp'){ df2_list[[sampleName]] <- TajD_df[TajD_df$coverage >= cov & TajD_df$bp >= 6000000 & TajD_df$bp <= 12000000, c('contig','bp',sampleName)] }
    
  }
  #head(df_list)
  
  # Merge dfs by contig and bp in common
  library(tidyverse)
  df2 <- df2_list %>% purrr::reduce(left_join, by = c('contig','bp'))
  df2$wd <- paste0(df2$contig,'-', df2$bp)
  df2$contig <- NULL
  df2$bp <- NULL
  head(df2)
  #tmp <- as.data.frame(t(df2),stringsAsFactors = F)
  #head(tmp)
  
  # Convert the long to the wide format
  library(reshape)
  df2_melt <- reshape2::melt(df2)
  df2_melt <- tidyr::separate(data = df2_melt, col = wd, into = c("contig", "bp"), sep = "-")  # split wd into contig and bp
  head(df2_melt)
  
  # Assign the color you want for each line
  df2_melt$season <- NA
  df2_melt[grepl('Spring',df2_melt$variable), 'season'] <- 'spring'
  df2_melt[grepl('Autumn',df2_melt$variable), 'season'] <- 'autumn'
  df2_melt$bp = as.integer(df2_melt$bp)
  df2_melt <- df2_melt[order(df2_melt$bp, decreasing = F), ]
  head(df2_melt)
  
  # Make the plot
  p2 <- ggplot(data=df2_melt, aes(x=bp, y=value, group=variable, colour=season)) + 
    geom_line(size=0.3) +
    #geom_line(size=0.3) +  # zoom
    theme_classic() + 
    #scale_linetype_manual(values=c(0,1,2,3,4,5,6,7,0,1,2,3,4,5,6,7)) +
    scale_color_manual(values=c('#984ea3','#ff7f00')) +
    #ggtitle(paste0('Comparison pi among spring and autumn pools', ', filter %cov >= ',cov)) +
    #scale_x_continuous(limits = c(8700000,9010000)) +  # zoom
    #scale_y_continuous(limits = c(0.000,0.005)) +  # zoom
    #scale_x_continuous(limits = c(8780000,9050000)) +  # zoom
    #scale_x_continuous(breaks=as.character(seq(min(df_melt$bp),max(df_melt$bp), 100000)))
    #xlab('position (bp)') + 
    ylab('Tajima\'s D') +
    theme(text = element_text(size=18),
          #axis.title.y=element_blank()
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank()
    )
  #p2
  
  #pdf('zoom_pi_drop.pdf')  # zoom
  #print(p1)
  #dev.off()
  
  # --------------------------------------------------------------------------------------
  # Process dAF data
  # --------------------------------------------------------------------------------------
  
  # Subset dAF values for the TSHR region
  if(region == 'chr15:8-10Mbp'){
    dAF_tshr <- dAF_df[dAF_df$CHROM == "chr15" & dAF_df$POS >= 8000000 & dAF_df$POS <= 10000000, c("POS","dAF_SpringAtlantic_vs_AutumnAtlantic")]
  }
  if(region == 'chr15:6-12Mbp'){
    dAF_tshr <- dAF_df[dAF_df$CHROM == "chr15" & dAF_df$POS >= 6000000 & dAF_df$POS <= 12000000, c("POS","dAF_SpringAtlantic_vs_AutumnAtlantic")]
  }
  colnames(dAF_tshr) <- c('POS','dAF')
  dim(dAF_tshr)
  #[1] 8123    2
  #[1] 27003     2
  head(dAF_tshr)
  class(dAF_tshr)
  
  p3 <- ggplot(data=dAF_tshr, aes(x=POS, y=dAF)) + 
    geom_point(size = 0.3, color="black") + 
    theme_classic() + 
    xlab('position (bp)') +
    theme(text = element_text(size=18))
  # Add title
  #if(region == 'chr15:8-10Mbp'){ p2 <- p2 + labs(title="dAF, Spring_vs_Fall, chr15:8-10Mbp", x = "bp", y = "dAF") }
  #if(region == 'chr15:6-12Mbp'){ p2 <- p2 + labs(title="dAF, Spring_vs_Fall, chr15:6-12Mbp", x = "bp", y = "dAF") }
  
  # Make a joint plot
  joint_plot <- p1 / p2 / p3
  #joint_plot
  
  pdf('plot_pi_TajD_dAF_TSHR_Herring_pools_spring_autumn_chr15.wd.10000.st.2000.chr15:8-10Mbp.pdf')
  #pdf('plot_pi_dAF_TSHR_Herring_pools_spring_autumn_chr15.wd.10000.st.2000.chr15:8-10Mbp_Not-subsampling.pdf')
  #pdf(paste0('Plot_pi_dAF_TSHR_',outputname,'.pdf'), 7,10)
  print(joint_plot)
  dev.off()
  
} else {
  print("plot_type or region missing...")
}

# --------------------------------------------------------------------------------------
# Comparison of pi values obtained when applying and not coverage subsampling
# --------------------------------------------------------------------------------------

# (Supplementary Figure)

require(tidyverse)
require(data.table)

df_list <- list()

# Load data
sub <- fread("./results/subsampling/DalGeB_Atlantic_Autumn.chr15.wd.10000.st.2000.idl.scov.pi.txt", header = F, data.table = F)
noSub <- fread("./results/not-subsampling/DalGeB_Atlantic_Autumn.chr15.wd.10000.st.2000.idl.pi.txt", header = F, data.table = F)
colnames(sub) <- c('contig','bp','Nsnps','cov','pi.sub')
colnames(noSub) <- c('contig','bp','Nsnps','cov','pi.noSub')

df_list[[1]] <- sub
df_list[[2]] <- noSub

df <- df_list %>% purrr::reduce(left_join, by = c('contig','bp'))
df <- df[df$cov.x >= 0.5, ]
head(df)
str(df)

pdf('Comparison_pi_when_subsampling_or_not_cov_popoolation1.pdf')
ggplot(df, aes(x=as.numeric(pi.sub), y=as.numeric(pi.noSub))) + geom_point(size=0.3) + geom_smooth(method = "lm") + 
  theme_classic() + xlab('pi (subsampled dataset)') + ylab('pi (not subsampled dataset)') +
  theme(text = element_text(size=18), 
        #axis.title.y=element_blank(),
        #axis.title.x=element_blank()#,
        axis.text.x=element_blank(), axis.text.y=element_blank(), axis.ticks.x=element_blank(), axis.ticks.y=element_blank()
  )
dev.off()

