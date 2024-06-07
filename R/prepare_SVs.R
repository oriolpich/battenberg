#' function that reads in SVs from VCFs for use in battenberg
#' (KT) NOTE: I would not necessarily trust var_count and coverage because it's not clear if this is the correct collumns
#' @param GRIPPS_SV_path vector including the paths to region specific GRIPPS output vcf
#' @param col_name wheter the midpoint of a beakpoint is used when the breakpoint (not SV) is not just one bp, only takes "mid" as input right now
#' @author KT
#' @export read.filter.gripps
read.filter.gripps <- function(GRIPPS_SV_path, col_name = "mid"){
  
  suppressPackageStartupMessages(require(VariantAnnotation))
  suppressPackageStartupMessages(require(StructuralVariantAnnotation))
  suppressPackageStartupMessages(require(dplyr))

  all_df <- list()
  for(i in 1:length(GRIPPS_SV_path)){

    sr_SV <- readVcf(GRIPPS_SV_path[[i]])
    gr <- breakpointRanges(sr_SV)
    gr <- gr[gr$FILTER == "PASS", ]
    gr_df <- as.data.frame(gr)
    bedpe <- breakpointgr2bedpe(gr)
    bedpe$strands <- paste0(bedpe$strand1, bedpe$strand2)
    
    SV_sr_pass <- sr_SV[rowRanges(sr_SV)$FILTER == "PASS"]
    SV_sr_pass_fix <- rowRanges(SV_sr_pass)
    SV_sr_pass_fix <- data.frame(SV_sr_pass_fix)
    SV_sr_pass_info <- info(SV_sr_pass)
    SV_sr_pass_info <- data.frame(SV_sr_pass_info)
    SV_sr_pass_geno <- geno(SV_sr_pass)$GT
    SV_sr_pass_geno <- data.frame(SV_sr_pass_geno)
    SV_sr_all <- cbind(SV_sr_pass_fix, SV_sr_pass_info, SV_sr_pass_geno)
    SV_sr_all$name <- rownames(SV_sr_all)
    SV_sr_all$REF <- NULL
    svtypes <- unique(SV_sr_all[, c("name", "EVENTTYPE", "TAF", "VF", "REF")])
    
    bedpe <- left_join(bedpe, svtypes[, c("name", "EVENTTYPE", "TAF", "VF", "REF")])
    
    # get the midpoint of each breakpoint
    bedpe$mid1 <- bedpe$end1 - ((bedpe$end1 - bedpe$start1)/2)
    bedpe$mid2 <- bedpe$end2 - ((bedpe$end2 - bedpe$start2)/2)
    
    bedpe <- bedpe[, c("chrom1", "start1", "end1", "mid1", "chrom2", "start2", "end2", "mid2", "strands", "EVENTTYPE", "TAF", "VF", "REF")]
    colnames(bedpe) <- c("chr1", "start1", "end1", "mid1", "chr2", "start2", "end2", "mid2", "strands", "SVTYPE", "VAF", "varcount", "coverage")

    bedpe[bedpe$strands == "+-", "SVTYPE"] <- "DEL"
    bedpe[bedpe$strands == "-+", "SVTYPE"] <- "DUP"
    bedpe[bedpe$strands == "++", "SVTYPE"] <- "h2hINV"
    bedpe[bedpe$strands == "--", "SVTYPE"] <- "t2tINV"
    bedpe[bedpe$strands == "INS", "SVTYPE"] <- "INS"
    bedpe[which(bedpe$chr1 != bedpe$chr2), "SVTYPE"] <- "BND"
    
    bedpe <- bedpe[bedpe$SVTYPE %in% c("h2hINV", "t2tINV", "INS", "BND"), ]
    
    # make into long format
    bedpe1 <- bedpe[, c("chr1", "start1", "end1", "mid1", "strands", "SVTYPE", "VAF", "varcount", "coverage")]
    bedpe2 <- bedpe[, c("chr2", "start2", "end2", "mid2", "strands", "SVTYPE", "VAF", "varcount", "coverage")]
    
    colnames(bedpe1) <- c("chr", "start", "end", "mid", "strands", "SVTYPE", "VAF", "varcount", "coverage")
    colnames(bedpe2) <- c("chr", "start", "end", "mid", "strands", "SVTYPE", "VAF", "varcount", "coverage")
    
    bedpe <- rbind(bedpe1, bedpe2)
    bedpe <- bedpe[bedpe$chr %in% paste0("chr", c(1:22, "X", "Y")), ]
    bedpe <- bedpe[, c("chr", col_name)]
    
    colnames(bedpe) <- c("chromosome", "position")
    bedpe$chromosome <- sub("chr", "", bedpe$chromosome)

    all_df[[i]] <- bedpe
    
  }

  all_df <- do.call(rbind, all_df)
  all_df <- unique(all_df)

  tumour <- strsplit(GRIPPS_SV_path, "/")[[1]][length(strsplit(GRIPPS_SV_path, "/")[[1]])]
  tumour <- sub("_.*", "", tumour)

  prior_breakpoints_file <- strsplit(GRIPPS_SV_path[[1]], "/")[[1]][1:(length(strsplit(GRIPPS_SV_path[[1]], "/")[[1]])-1)]
  prior_breakpoints_file <- paste(c(prior_breakpoints_file), collapse = "/")
  prior_breakpoints_file <- paste0(prior_breakpoints_file, "/", tumour, ".gripss.filtered.somatic.txt")

  write.table(all_df, prior_breakpoints_file, col.names = T, row.names = F, quote = F, sep = "\t")
  return(prior_breakpoints_file)
}


