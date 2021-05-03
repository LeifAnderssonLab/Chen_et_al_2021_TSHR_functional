# Using the illumina genome version to compare the different TSHR alleles
# Code developed by Mats Pettersson, e-mail: mats.pettersson@imbim.uu.se
# Uppsala University

#TSHR_gene_model <- read.table("~/Projects/Herring/data/TSHR/TSHR_Ch_v2.blastout", stringsAsFactors=F)
#TSHR_gene_model <- TSHR_gene_model[TSHR_gene_model[,3] >= 95, ]
#TSHR_SNPs <- read.table("~/Projects/Herring/data/TSHR/TSHR_SNPs_Ch_v2.blastout", stringsAsFactors=F)
#TSHR_SNPs <- TSHR_SNPs[TSHR_SNPs[,3] > 98 & TSHR_SNPs[,2] == "Chr17", ]

# Load data
TSHR_BAC_raw <- read.table("~/Projects/Herring/data/TSHR/BAC_P19L22_PB_and_Ilu_blastout.txt", stringsAsFactors=F)
TSHR_BAC <- TSHR_BAC_raw[TSHR_BAC_raw[,2] == "scaffold1420" & TSHR_BAC_raw[,3] > 97,]

TSHR_gene_model_ilu <- read.table("~/Projects/Herring/data/TSHR/TSHR_cDNA_Ilu_PB_purged.blastout", stringsAsFactors=F)
TSHR_gene_model_ilu <- TSHR_gene_model_ilu[TSHR_gene_model_ilu[,3] >= 95 & TSHR_gene_model_ilu[,2] == "scaffold1420", ]

#TSHR_Ch_v2_GR <- GRanges(seqnames = TSHR_gene_model[,2], ranges = IRanges(start = TSHR_gene_model[,9], end = TSHR_gene_model[,10]))
#TSHR_Ch_v2_GR_red <- reduce(TSHR_Ch_v2_GR, with.revmap =T)
#revmap <- mcols(TSHR_Ch_v2_GR_red)$revmap
#TSHR_gene_model[revmap[[11]],]
#red_exon_list <- c("5'UTR, Exon 1 copy", "Exon 1", "Exon 2", "Exon 3", "Exon 4", "Exon 5", "Exon 6", "Exon 7", "Exon 8", "Exon 9", "Exon 10")
#red_exon_list <- c("Exon 1", "Exon 2", "Exon 3", "Exon 4", "Exon 5", "Exon 6", "Exon 7", "Exon 8", "Exon 9", "Exon 10")

#TSHR_Ch_v2_GR_red$exon <- red_exon_list
#TSHR_Ch_v2_df <- as.data.frame(TSHR_Ch_v2_GR_red)
#TSHR_Ch_v2_df[,"col"] <- "blue"

#pdf("~/Projects/Herring/doc/Ch_genome_v2/TSHR_BAC.pdf", width = 12, height = 7)
plot(y = 412, x = 9341554, xlim = c(0, 2e5), ylim = range(TSHR_BAC[,7:8]), type = "n", xlab = "", ylab = "BAC_P19L22")
segments(x0 = TSHR_BAC[,9], x1 = TSHR_BAC[,10], y0 = TSHR_BAC[,7], y1 = TSHR_BAC[,8])

lines( x = c(TSHR_Ch_v2_df$start[3], TSHR_Ch_v2_df$end[11]), y = c(9e4, 9e4), col = "black", lwd = 3)
rect(xleft = TSHR_Ch_v2_df$start, xright = TSHR_Ch_v2_df$end, ybottom = 8.8e4, ytop = 9.2e4, col = TSHR_Ch_v2_df$col, border = NA)
text_x_vec <- rowMeans(TSHR_Ch_v2_df[,c("start", "end")])
points(x = rowMeans(TSHR_SNPs[,9:10]), y = rep(8.6e4, dim(TSHR_SNPs)[1]), pch = 16)
text(x = rowMeans(TSHR_SNPs[,9:10]), y =c(3, 2, 1, 3,2,3,2) * 4.5e3 + 6.8e4, labels = sub("_scaff.+", "", TSHR_SNPs[,1]), cex = 0.9 )
#dev.off()


# Parallel plot of the BAC (spring haplotype) & Ilumina v1.2 (autumn haplotype)
plot(x = 0, y = 0, xlim = c(0.5e5, 2e5), ylim = c(0, 1.3e5), xlab = "Position on scaffold1420", ylab = "Position on BAC_P19L22", main = "", type  = "n")
segments(y0 = TSHR_BAC[,7], y1 = TSHR_BAC[,8], x0 = TSHR_BAC[,9], x1 = TSHR_BAC[,10], col = c("red", "blue")[as.numeric(TSHR_BAC[,10] > TSHR_BAC[,9]) + 1])

# Picking the target alignment blocks 
BAC_target_blocks <- TSHR_BAC[TSHR_BAC[,7] < 7.e4 & TSHR_BAC[,8] > 2.9e4,]
segments(y0 = BAC_target_blocks[,7], y1 = BAC_target_blocks[,8], x0 = BAC_target_blocks[,9], x1 = BAC_target_blocks[,10], col = c("darkorchid"), lwd  =3)

# Shifting the position vectors, and inverting the bac ones
chr_start_vec <- BAC_target_blocks[,9]
chr_end_vec <- BAC_target_blocks[,10]

BAC_start_vec <- -(BAC_target_blocks [,7] - min(BAC_target_blocks [,7:8])) + max(BAC_target_blocks[,9:10])
BAC_end_vec <- -(BAC_target_blocks [,8] - min(BAC_target_blocks [,7:8])) + max(BAC_target_blocks[,9:10])

# The same for the connections
connection_df <- data.frame(bac = c(39673,40901, 18241, 42679,30747,36395,46960, 88740), ilu = c(131349 + (40901 - 39673), 131349 - (40901 - 39673),148451, 128334, 135875, 135880, 123127,81670), 
                            type = rep(c("ilu_dup", "aln_block_1", "ilu_del", "aln_block_2"),each = 2 ), 
                            col = rep(c("grey50","black", "grey50", "black"), each = 2), stringsAsFactors = F)
connection_df[,"plot_bac"] <- -(connection_df[,"bac"]- min(BAC_target_blocks [,7:8])) + max(BAC_target_blocks[,9:10])

# Other features
BEL_pos  <- c(30747, 36395) 
BEL_pos  <- -(BEL_pos- min(BAC_target_blocks [,7:8])) + max(BAC_target_blocks[,9:10])
BAC_range <- c(18241,88740 )
BAC_range  <- -(BAC_range- min(BAC_target_blocks [,7:8])) + max(BAC_target_blocks[,9:10])
BAC_dup_pos  <- c(39673, 40901) 
BAC_dup_pos  <- -(BAC_dup_pos- min(BAC_target_blocks [,7:8])) + max(BAC_target_blocks[,9:10])
ILU_dup_pos <- c(131349 - (40901 - 39673), 131349, 131349 + (40901 - 39673))
ILU_range <- c(81670, 148451)

TSHR_SNPs[,"ilu_pos"] <- as.numeric(sub(".+d1420_([0-9]+)", "\\1", TSHR_SNPs[,1]))

pdf("~/Projects/Herring/doc/TSHR/BAC_P19L22_vs_ilu_v2.pdf", width = 12, height = 7)
BAC_level <- 1
box_offset <- 0.1
ILU_level <- 0.6
plot(x = 0, y = 0, xlim = c(min(BAC_target_blocks[,9:10])-1e4, max(BAC_target_blocks[,9:10]) + 0.5e4), ylim = c(0, 3), ylab = "", xlab = "Position on scaffold 1420", main = "", type  = "n", axes = F)
segments(y0 = ILU_level + box_offset, y1 = BAC_level - box_offset, x0 = connection_df[,"ilu"], x1 = connection_df[,"plot_bac"], col =  connection_df[,"col"], lwd  = 2)
rect(ybottom = BAC_level - 0.5*box_offset, ytop = BAC_level + 0.5*box_offset, xleft = BAC_range[1], xright = BAC_range[2], col = c("steelblue3"), border = "steelblue3")
rect(ybottom = ILU_level - 0.5*box_offset, ytop = ILU_level + 0.5*box_offset, xleft = ILU_range[1], xright = ILU_range[2], col = c("darkorchid3"), border = "darkorchid3")

rect(ybottom = ILU_level - box_offset, ytop = ILU_level + box_offset,  xleft = chr_start_vec, xright = chr_end_vec, col = "darkorchid", border = "darkorchid")
rect(ybottom = BAC_level - box_offset, ytop = BAC_level + box_offset, xleft = BAC_start_vec, xright = BAC_end_vec, col = c("steelblue"), border = "steelblue")

rect(ybottom = BAC_level - 1.5*box_offset, ytop = BAC_level + 1.5*box_offset, xleft = BEL_pos[1], xright = BEL_pos[2], col = c("steelblue1"))
rect(ybottom = BAC_level - 1.5*box_offset, ytop = BAC_level + 1.5*box_offset, xleft = BAC_dup_pos[1], xright = BAC_dup_pos[2], col = c("steelblue2"))
rect(ybottom = ILU_level - 1.5*box_offset, ytop = ILU_level + 1.5*box_offset, xleft = ILU_dup_pos[1:2], xright = ILU_dup_pos[2:3], col = c("darkorchid1", "darkorchid2"))
rect(xleft = TSHR_gene_model_ilu[,9], xright = TSHR_gene_model_ilu[,10], ybottom = ILU_level- 3*box_offset, ytop = ILU_level - 2*box_offset, col = "black")
lines(x = range(TSHR_gene_model_ilu[9:10]), y = rep(ILU_level - 2.5* box_offset, 2))
segments(x0 = TSHR_SNPs[,"ilu_pos"], y0 = ILU_level - 4* box_offset, y1 = ILU_level-3.5*box_offset, lwd = 4, col = "firebrick")
axis(1, pos = 0, cex.axis = 2, lwd = 3)
dev.off()

# Excising mysterious "spacer" from the two versions
Ilu_fasta <- readDNAStringSet("~/Projects/Herring/data/genomic_data/assembly_files/references/ClupeaIlu_v1.2.fasta.gz") #Genome release v1.2; Martinez-Barrio et al 2016
BAC_fasta <- readDNAStringSet("~/Projects/Herring/data/TSHR/BAC_P19L22_sequence.fa")
spacer_blocks <- DNAStringSet()
tmp_seq <- subseq(BAC_fasta[1], start = 42679, end = 46960)
names(tmp_seq) <- "BAC"
spacer_blocks <- append(spacer_blocks, tmp_seq)

tmp_seq <-  subseq(Ilu_fasta["scaffold1420"], start = 123127, end = 128334)
names(tmp_seq) <- "ILU"
spacer_blocks <- append(spacer_blocks, tmp_seq)

tmp_seq <- subseq(BAC_fasta[1], start = 18241, end = 42679)
names(tmp_seq) <- "BAC_validation"
spacer_blocks <- append(spacer_blocks, tmp_seq)

tmp_seq <- subseq(Ilu_fasta["scaffold1420"], start = 128334, end = 148451)
names(tmp_seq) <- "ILU_validation"
spacer_blocks <- append(spacer_blocks, tmp_seq)

writeXStringSet(spacer_blocks, file = "~/Projects/Herring/data/TSHR/spacers.fasta")

#segments(y0 = 5, x0 = BAC_start_vec, x1 = BAC_end_vec, col = c("red", "blue")[as.numeric(TSHR_BAC[,10] > TSHR_BAC[,9]) + 1])
#segments(y0 = 1, y1 = 5, x0 = chr_start_vec, x1 = BAC_start_vec, col = "darkorchid2")
#segments(y0 = 1, y1 = 5, x0 = chr_end_vec, x1 = BAC_end_vec, col = "firebrick")

# Forskarhjälpen
require(GenABEL)
load("~/Projects/Herring/data/SNP_chip/Forskarhjalpen/fh_LD_all_original_SNPs_data.Rdata")
load("~/Projects/Herring/data/SNP_chip/Forskarhjalpen/SNP_map.Rdata")
TSHR_SNP_map_raw<- SNP_map[SNP_map[,"scaffold"] == 1420,]
TSHR_SNP_map <- TSHR_SNP_map_raw[order(TSHR_SNP_map_raw[,"pos"]),]

TSHR_SNP_map[TSHR_SNP_map$pos %in% c(133030, 108994, 109295,100524,137485),]
TSHR_1_gt <- as.numeric(fh_LD_all_original_SNPs_data[,"AX-99143205"]@gtdata)
TSHR_1_gt[is.na(TSHR_1_gt)] <- 99

# Connected scaffolds

table(SNP_map_lo[SNP_map_lo[,"seqnames"] == "chr15" & SNP_map_lo$SNP_HiC_pos > 5.0e6 &  SNP_map_lo$SNP_HiC_pos < 15e6,]$chr)
chr15_5_15MB_scaffs <- names(table(SNP_map_lo[SNP_map_lo[,"seqnames"] == "chr15" & SNP_map_lo$SNP_HiC_pos > 5.0e6 &  SNP_map_lo$SNP_HiC_pos < 15e6,]$chr))
chr15_5_15MB_scaffs <- as.numeric(sub("scaffold", "", chr15_5_15MB_scaffs))

TSHR_SNP1 <- LD_plot_by_allele("AX-99143205", chr15_5_15MB_scaffs, plot_prefix = "~/Projects/Herring/doc/TSHR/200827", xlim = c(8.7e6, 9.1e6), ylim = c(0,0.6))
#TSHR_SNP1 <- LD_plot_by_allele("AX-99143205", 1420, plot_prefix = "~/Projects/Herring/doc/TSHR/")
TSHR_SNP2 <- LD_plot_by_allele("AX-99034443", 1420, plot_prefix = "~/Projects/Herring/doc/TSHR/")
TSHR_SNP3 <- LD_plot_by_allele("AX-99136256", 1420, plot_prefix = "~/Projects/Herring/doc/TSHR/")
TSHR_SNP4 <- LD_plot_by_allele("AX-99110137", 1420, plot_prefix = "~/Projects/Herring/doc/TSHR/")
TSHR_SNP5 <- LD_plot_by_allele("AX-99030516", 1420, plot_prefix = "~/Projects/Herring/doc/TSHR/", xlim = c(8.7e6, 9.1e6))
cor(cbind(TSHR_SNP1$gt, TSHR_SNP2$gt, TSHR_SNP3$gt, TSHR_SNP4$gt, TSHR_SNP5$gt), use = "complete")


target_scaff_vec <- c(1420) # Contains the TSHR locus
TSHR_subset_filter <- fh_LD_all_original_SNPs_data@gtdata@chromosome %in% target_scaff_vec
LD_marker_map <- as.data.frame(cbind(as.character(fh_LD_all_original_SNPs_data[,TSHR_subset_filter]@gtdata@chromosome), fh_LD_all_original_SNPs_data[,TSHR_subset_filter]@gtdata@map), stringsAsFactors = F)
names(LD_marker_map) <- c("chr", "pos")
LD_marker_map[, "chr"] <- paste("scaffold", LD_marker_map[, "chr"], sep = "")
LD_marker_map[, "pos"] <- as.numeric(LD_marker_map[, "pos"])
LD_marker_map[, "LD_idx"] <- 1:dim(LD_marker_map)[1]
LD_lo_df <- HiC_liftover(scaffold_data=LD_marker_map, liftover_df=Ch_v2.0.2_v_Ilu_satsuma, chr_size_df=Ch_v2.0.2_sizes, lo_cols="LD_idx")
LD_lo_df <- LD_lo_df[LD_lo_df$SNP_HiC_pos > 8e6,]

#TSHR_subset_filter <- (fh_LD_all_original_SNPs_data@gtdata@chromosome %in% target_scaff_vec)
GenABEL_TSHR_subset <- fh_LD_all_original_SNPs_data[,TSHR_subset_filter]
GenABEL_TSHR_ordered <- GenABEL_TSHR_subset[,LD_lo_df$LD_idx]
GenABEL_TSHR_ordered@gtdata@map <- LD_lo_df$SNP_HiC_pos
#GenABEL_TSHR_ordered@gtdata@chromosome <- factor(12)
GenABEL_TSHR_ordered@gtdata@chromosome <- factor(15)
#for comparison, from Rhodopsin_FH.R
TSHR_by_school_reorder

THSR_2_inds <- TSHR_1_gt == 2
THSR_0_inds <- TSHR_1_gt == 0
# Heterozygosity
TSHR1_type2_GTs <- as.numeric(GenABEL_TSHR_ordered[THSR_2_inds,]@gtdata)  # Autumn version
pdf("~/Projects/Herring/doc/TSHR/TSHR1_2_inds_HZ.pdf")
plot(y = colSums(TSHR1_type2_GTs == 1, na.rm = T)/dim(TSHR1_type2_GTs)[1], x = GenABEL_TSHR_ordered@gtdata@map, col = "red", pch = 16, xlab = "Position", ylab = "Fraction of heterozygous individuals", main = "TSHR_1 GT == 2", xlim = c(8.5e6, 9.2e6))
dev.off()

TSHR1_type1_GTs <- as.numeric(GenABEL_TSHR_ordered[THSR_0_inds,]@gtdata) # Spring version
pdf("~/Projects/Herring/doc/TSHR/TSHR1_0_inds_HZ.pdf")
plot(y = colSums(TSHR1_type1_GTs == 1, na.rm = T)/dim(TSHR1_type1_GTs)[1], x = GenABEL_TSHR_ordered@gtdata@map, col = "red", pch = 16, xlab = "Position", ylab = "Fraction of heterozygous individuals", main = "TSHR_1 GT == 0", xlim = c(8.5e6, 9.2e6))
dev.off()

TSHR1_all_GTs <- as.numeric(GenABEL_TSHR_ordered@gtdata)
pdf("~/Projects/Herring/doc/TSHR/TSHR1_all_inds_HZ.pdf")
plot(y = colSums(TSHR1_all_GTs == 1, na.rm = T)/dim(TSHR1_all_GTs)[1], x = GenABEL_TSHR_ordered@gtdata@map, col = "red", pch = 16, xlab = "Position", ylab = "Fraction of heterozygous individuals", main = "All individulas", xlim = c(8.5e6, 9.2e6))
dev.off()

# LD
all_TSHR_LD <- plot_subset_LD_v2(marker_lo_map = LD_lo_df, t_marker_vec = TSHR_subset_filter, png_file = "~/Projects/Herring/doc/TSHR/All_inds_LD.png", main_title = "TSHR LD; All individuals", focal_pos = LD_lo_df["AX-99143205",]$SNP_HiC_pos)
TSHR_ld_decay_all <- plot.LD.decay(data = GenABEL_TSHR_ordered, dmin = 0, dmax = 2e6, N = 100, pdf_file = "~/Projects/Herring/doc/TSHR/All_inds_LD_decay.pdf")

TSHR_0_LD <- plot_subset_LD_v2(marker_lo_map = LD_lo_df, ind_list = THSR_0_inds, t_marker_vec = TSHR_subset_filter, png_file = "~/Projects/Herring/doc/TSHR/TSHR_0_LD.png", main_title = "TSHR LD; TSHR_1 == 0", focal_pos = LD_lo_df["AX-99143205",]$SNP_HiC_pos)
TSHR_0_ld_decay <- plot.LD.decay(data = GenABEL_TSHR_ordered[THSR_0_inds,], dmin = 0, dmax = 2e6, N = 100, pdf_file = "~/Projects/Herring/doc/TSHR/TSHR_0_LD_decay.pdf")

TSHR_2_LD <- plot_subset_LD_v2(marker_lo_map = LD_lo_df, ind_list = THSR_2_inds, t_marker_vec = TSHR_subset_filter, png_file = "~/Projects/Herring/doc/TSHR/TSHR_2_LD.png", main_title = "TSHR LD; TSHR_1 == 2", focal_pos = LD_lo_df["AX-99143205",]$SNP_HiC_pos)
TSHR_2_ld_decay <- plot.LD.decay(data = GenABEL_TSHR_ordered[THSR_2_inds,], dmin = 0, dmax = 2e6, N = 100, pdf_file = "~/Projects/Herring/doc/TSHR/TSHR_2_LD_decay.pdf")


# Extracting from new annotation
TSHR_reg_gtf <- cluhar_v2.0.2_gtf[cluhar_v2.0.2_gtf@seqnames == "chr15" & cluhar_v2.0.2_gtf@ranges@start > 8.80e6 & cluhar_v2.0.2_gtf@ranges@start < 9.0e6 & cluhar_v2.0.2_gtf@elementMetadata$type == "gene" ]
TSHR_exon_gtf <- cluhar_v2.0.2_gtf[cluhar_v2.0.2_gtf@elementMetadata$gene_id == "ENSCHAG00000002528" & cluhar_v2.0.2_gtf@elementMetadata$type == "exon"]

# Extracting seqence block surrounding TSHR for conservation estimate
TSHR_block_GR <- GRanges(seqnames = "chr15", ranges = IRanges(start = 8.5e6, end = 9.5e6))
TSHR_block_chr15_8.5_to_9.5_Mb <- Ch_v2.0.2[TSHR_block_GR]
writeXStringSet(TSHR_block_chr15_8.5_to_9.5_Mb, filepath = "~/Projects/Herring/data/TSHR/TSHR_block_chr15_8.5_to_9.5_Mb.fasta")

# 200k block for Online BLAST
TSHR_block_GR <- GRanges(seqnames = "chr15", ranges = IRanges(start = 8.8e6, end = 8.9999e6))
TSHR_block_chr15_8.8_to_9_Mb <- Ch_v2.0.2[TSHR_block_GR]
writeXStringSet(TSHR_block_chr15_8.8_to_9_Mb, filepath = "~/Projects/Herring/data/TSHR/TSHR_block_chr15_8.8_to_9_Mb.fasta")

# Collecting some TSHR regions from teleosts
TSHR_regions <- readDNAStringSet("~/Projects/Herring/data/TSHR/Other_fish_TSHR_regions/TSHR_regions.fa")
TSHR_regions_DNAbin <- as.DNAbin(TSHR_regions)
require(fishtree)
TSHR_fish <- c("Oryzias latipes", "Tetraodon nigroviridis", "Gasterosteus aculeatus", "Danio rerio", "Lepisosteus oculatus", "Esox lucius", "Amphiprion percula", "Neolamprologus brichardi", "Poecilia reticulata", "Carassius auratus", "Denticeps clupeoides", "Clupea harengus")
THSR_guide_tree <- fishtree_phylogeny(species = TSHR_fish)

#phyloFit --tree "((((((((guppy,a_molly),medaka),lyretail_chiclid),o_clownfish),(stickleback,tetraodon)),northern_pike),(goldfish,(atlantic_herring,denticle_herring))),spotted_gar)Anc0;" TSHR_regions.maf > TSHR_regions.mod
#phastCons TSHR_regions.maf phyloFit.mod > TSHR_region_socres.wig &

require(rtracklayer)
TSHR_scores <- import.wig(con = "~/Projects/Herring/data/TSHR/TSHR_region_socres.wig")
TSHR_SNPs <- read.table("~/Projects/Herring/data/HiC_assemblies/FALCON_phase/candidate/blast_eval/Ch_v2_candidate_190130_TSHR_SNPs.blastout", stringsAsFactors=F)
TSHR_SNPs <- TSHR_SNPs[TSHR_SNPs[,3] > 98 & TSHR_SNPs[,2] == "Chr_15", ]

TSHR_reg_start <- 8786015
TSHR_scores_adj <- shift(TSHR_scores, (TSHR_reg_start - 1)) 
pdf("~/Projects/Herring/doc/TSHR/TSHR_PhastCons_v2.pdf", width = 10)
plot(x=TSHR_scores_adj@ranges@start, y = TSHR_scores_adj@elementMetadata$score, ylim = c(0,2), xlim = c(8.86e6, 8.915e6), xlab = "Position on Chr 15", ylab = "PhastCons score", cex = 0.6, pch = 16, col = "grey30")
segments(x0 = TSHR_reg_gtf@ranges@start, x1 = TSHR_reg_gtf@ranges@start + TSHR_reg_gtf@ranges@width, y0 = 1.5, col = "black")
segments(x0 = TSHR_exon_gtf@ranges@start, x1 = TSHR_exon_gtf@ranges@start + TSHR_exon_gtf@ranges@width, y0 = 1.5, col = "red", lwd  = 6)
segments(y1 = seq(from = 1.40, to = 1.15, length.out = dim(TSHR_SNPs)[1]), x0 = SNP_score@ranges@start, y0 = SNP_score@elementMetadata$score, col = "grey70")
points(x = rowMeans(TSHR_SNPs[,9:10]), y = seq(from = 1.15, to = 1.40, length.out = dim(TSHR_SNPs)[1]), pch = 16, col = "darkorchid")
#text(rowMeans(TSHR_SNPs[,9:10]), y = seq(from = 1.15, to = 1.40, length.out = dim(TSHR_SNPs)[1]), pos = c(3,1), cex = 0.7, labels = sub("TSHR_[0-9]_", "", TSHR_SNPs$V1))
text(rowMeans(TSHR_SNPs[2:5,9:10]), y = seq(from = 1.15, to = 1.40, length.out = dim(TSHR_SNPs)[1])[2:5], pos = c(1,3), cex = 0.7, labels = sub("TSHR_[0-9]_", "", TSHR_SNPs$V1[2:5]))
SNP_score_idx <- which((TSHR_scores_adj@ranges@start) %in% rowMeans(TSHR_SNPs[,9:10]))
SNP_score <- TSHR_scores_adj[SNP_score_idx]
points(x = SNP_score@ranges@start, y = SNP_score@elementMetadata$score, col = "darkorchid", pch = 16, cex = 1.5)
segments(x0 = TSHR_score_peaks@ranges@start, x1 = TSHR_score_peaks@ranges@start + TSHR_score_peaks@ranges@width, y0 = 1.7, col = "darkorange", lwd = 2)
dev.off()
TSHR_SNPs$PhastCons <- SNP_score@elementMetadata$score[match(rowMeans(TSHR_SNPs[,9:10]), SNP_score@ranges@start)]
names(TSHR_SNPs)[1] <- "SNP_ID"
TSHR_SNPs$v2.0.2_pos <- rowMeans(TSHR_SNPs[,9:10])
write.table(x = TSHR_SNPs[,c("SNP_ID", "PhastCons")], file = "~/Projects/Herring/doc/TSHR/TSHR_PhastCons_scores.txt", row.names = F, quote = F, sep = "\t")
save(TSHR_scores_adj, TSHR_SNPs, file = "~/Projects/Herring/doc/TSHR/TSHR_PhastCons.RData")

# Conservation peaks
TSHR_score_peaks <- reduce(TSHR_scores_adj[TSHR_scores_adj@elementMetadata$score > 0.2], min.gapwidth = 500)
seqinfo(TSHR_score_peaks,1:1) <- Seqinfo("chr15")
write.table(x = as.data.frame(TSHR_score_peaks), file = "~/Projects/Herring/doc/TSHR/TSHR_PhastCons_peaks.txt", row.names = F, quote = F, sep = "\t")


# Support functions
plot_subset_LD_v2 <- function(marker_lo_map, t_marker_vec,  ind_list = NULL, snp_data = fh_LD_all_original_SNPs_data, png_file, main_title = "", focal_pos = NULL){
  if(!is.null(ind_list)){
    tmp_LD_data <- snp_data[ind_list,t_marker_vec]
  } else {
    tmp_LD_data <- snp_data[,t_marker_vec]
  }
  sub_LD <- r2fast(tmp_LD_data)
  sub_LD[lower.tri(sub_LD)] <- t(sub_LD)[lower.tri(sub_LD)]
  png(png_file, width = 1000, height = 1000)
  image(x = marker_lo_map[,"SNP_HiC_pos"], y = marker_lo_map[,"SNP_HiC_pos"], z = 1-sub_LD[marker_lo_map[,"LD_idx"], marker_lo_map[,"LD_idx"]], main = main_title, xlab = "Position", ylab = "Position")
  if(!is.null(focal_pos)){
    abline(v = focal_pos)
    abline(h = focal_pos)
  }
  dev.off()
  return(sub_LD)
}

LD_plot_by_allele <- function(focal_marker, focal_scaffs, plot_prefix = "./", lo_satsuma = Ch_v2.0.2_v_Ilu_satsuma, lo_sizes = Ch_v2.0.2_sizes,  GenABEL_data = fh_LD_all_original_SNPs_data, ...){
  
  focal_subset_filter <- GenABEL_data@gtdata@chromosome %in% focal_scaffs
  focal_LD_marker_map <- as.data.frame(cbind(as.character(GenABEL_data[,focal_subset_filter]@gtdata@chromosome), GenABEL_data[,focal_subset_filter]@gtdata@map), stringsAsFactors = F)
  names(focal_LD_marker_map) <- c("chr", "pos")
  focal_LD_marker_map[, "chr"] <- paste("scaffold", focal_LD_marker_map[, "chr"], sep = "")
  focal_LD_marker_map[, "pos"] <- as.numeric(focal_LD_marker_map[, "pos"])
  focal_LD_marker_map[, "LD_idx"] <- 1:dim(focal_LD_marker_map)[1]
  focal_LD_lo_df <- HiC_liftover(scaffold_data=focal_LD_marker_map, liftover_df=lo_satsuma, chr_size_df=lo_sizes, lo_cols="LD_idx")
  focal_chr <- names(table(focal_LD_lo_df$seqnames)[order(table(focal_LD_lo_df$seqnames), decreasing = T)][1])
  focal_LD_lo_df <- focal_LD_lo_df[focal_LD_lo_df$seqnames == focal_chr, ]
  
  GenABEL_focal_subset <- GenABEL_data[,focal_subset_filter]
  GenABEL_focal_ordered <- GenABEL_focal_subset[,focal_LD_lo_df$LD_idx]
  GenABEL_focal_ordered@gtdata@map <- focal_LD_lo_df$SNP_HiC_pos
  #focal_chr <- sub("chr", "", table(focal_LD_lo_df$seqnames)[order(table(focal_LD_lo_df$seqnames))][1])
  GenABEL_focal_ordered@gtdata@chromosome <- factor(focal_chr)
  
  focal_gt <- as.numeric(GenABEL_focal_ordered[,focal_marker]@gtdata)
  type_2_inds <- focal_gt == 2
  type_0_inds <- focal_gt == 0
  
  type2_GTs <- as.numeric(GenABEL_focal_ordered[type_2_inds,]@gtdata)
  pdf(paste(plot_prefix, focal_marker,"_2_inds_HZ.pdf", sep = ""))
  plot(y = colSums(type2_GTs == 1, na.rm = T)/dim(type2_GTs)[1], x = GenABEL_focal_ordered@gtdata@map, col = "red", pch = 16, xlab = "Position", ylab = "Fraction of heterozygous individuals", main = paste(focal_marker, "GT == 2\n n =", sum(type_2_inds, na.rm = T)), ...)
  abline(v = focal_LD_lo_df[focal_marker,]$SNP_HiC_pos, col = "grey50")
  dev.off()
  
  type0_GTs <- as.numeric(GenABEL_focal_ordered[type_0_inds,]@gtdata)
  pdf(paste(plot_prefix, focal_marker,"_0_inds_HZ.pdf", sep = ""))
  plot(y = colSums(type0_GTs == 1, na.rm = T)/dim(type0_GTs)[1], x = GenABEL_focal_ordered@gtdata@map, col = "red", pch = 16, xlab = "Position", ylab = "Fraction of heterozygous individuals", main = paste(focal_marker, "GT == 0\n n =", sum(type_0_inds, na.rm = T)), ...)
  abline(v = focal_LD_lo_df[focal_marker,]$SNP_HiC_pos, col = "grey50")
  dev.off()
  
  all_LD <- plot_subset_LD_v2(marker_lo_map = focal_LD_lo_df, t_marker_vec = focal_subset_filter, png_file = paste(plot_prefix, focal_marker,"_all_LD.png", sep = ""), main_title = "LD; All individuals", focal_pos = focal_LD_lo_df[focal_marker,]$SNP_HiC_pos)
  #all_ld_decay <- plot.LD.decay(data = GenABEL_focal_ordered, dmin = 0, dmax = 2e6, N = 100, pdf_file = paste(plot_prefix, focal_marker,"_all_LD_decay.pdf", sep = ""))
  
  type_0_LD <- plot_subset_LD_v2(marker_lo_map = focal_LD_lo_df, ind_list = type_0_inds, t_marker_vec = focal_subset_filter, png_file = paste(plot_prefix, focal_marker,"_type_0_LD.png", sep = ""), main_title = paste("LD; ", focal_marker, "== 0"), focal_pos = focal_LD_lo_df[focal_marker,]$SNP_HiC_pos)
  #type_0_ld_decay <- plot.LD.decay(data = GenABEL_focal_ordered[type_0_inds,], dmin = 0, dmax = 2e6, N = 100, pdf_file =  paste(plot_prefix, focal_marker,"_type_0_LD_decay.pdf", sep = ""))
  
  type_2_LD <- plot_subset_LD_v2(marker_lo_map = focal_LD_lo_df, ind_list = type_2_inds, t_marker_vec = focal_subset_filter, png_file =  paste(plot_prefix, focal_marker,"_type_2_LD.png", sep = ""), main_title = paste("LD; ", focal_marker, "== 2"), focal_pos = focal_LD_lo_df[focal_marker,]$SNP_HiC_pos)
  #type_2_ld_decay <- plot.LD.decay(data = GenABEL_focal_ordered[type_2_inds,], dmin = 0, dmax = 2e6, N = 100, pdf_file = paste(plot_prefix, focal_marker,"_type_2_LD_decay.pdf", sep = ""))
  
  return(invisible(list(gt = focal_gt, all = all_LD, type_0 = type_0_LD, type_2 = type_2_LD, snp = focal_marker)))
}

plot_subset_LD_v2 <- function(marker_lo_map, t_marker_vec,  ind_list = NULL, snp_data = fh_LD_all_original_SNPs_data, png_file, main_title = "", focal_pos = NULL){
  if(!is.null(ind_list)){
    tmp_LD_data <- snp_data[ind_list,t_marker_vec]
  } else {
    tmp_LD_data <- snp_data[,t_marker_vec]
  }
  sub_LD <- r2fast(tmp_LD_data)
  sub_LD[lower.tri(sub_LD)] <- t(sub_LD)[lower.tri(sub_LD)]
  png(png_file, width = 1000, height = 1000)
  image(x = marker_lo_map[,"SNP_HiC_pos"], y = marker_lo_map[,"SNP_HiC_pos"], z = 1-sub_LD[marker_lo_map[,"LD_idx"], marker_lo_map[,"LD_idx"]], main = main_title, xlab = "Position", ylab = "Position")
  if(!is.null(focal_pos)){
    abline(v = focal_pos)
    abline(h = focal_pos)
  }
  dev.off()
  return(sub_LD)
}
