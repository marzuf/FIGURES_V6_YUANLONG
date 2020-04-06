# Rscript cmp_auc_rescaled.R

outFolder <- "CMP_AUC_RESCALED"
dir.create(outFolder)


plotType <- "svg"
myHeightGG <- 7
myWidthGG <- 9
plotCex <- 1.4

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../FIGURES_V2_YUANLONG/settings.R")


require(foreach)
require(doMC)
registerDoMC(40)
require(ggplot2)
require(ggpubr)

buildTable <- TRUE

fcc_auc_dt <- get(load(file.path("../FIGURES_V2_YUANLONG/BARPLOT_FCC_AUC_RATIO/all_dt.Rdata")))


auc_v0 <- fcc_auc_dt$fcc_auc[fcc_auc_dt$hicds == "ENCSR079VIJ_G401_40kb" & fcc_auc_dt$exprds == "TCGAkich_norm_kich"]
auc_obs <- get(load(file.path("../FIGURES_V2_YUANLONG/FCC_WAVE_PLOT/ENCSR079VIJ_G401_40kb_TCGAkich_norm_kich_auc_obs.Rdata")))
auc_perm <- get(load(file.path("../FIGURES_V2_YUANLONG/FCC_WAVE_PLOT/ENCSR079VIJ_G401_40kb_TCGAkich_norm_kich_auc_permutQt.Rdata")))
stopifnot(round(auc_v0, 4) == round(auc_obs/auc_perm, 4))

auc_v0 <- fcc_auc_dt$fcc_auc[fcc_auc_dt$hicds == "LG1_40kb" & fcc_auc_dt$exprds == "TCGAluad_mutKRAS_mutEGFR"]
auc_obs <- get(load(file.path("../FIGURES_V2_YUANLONG/FCC_WAVE_PLOT/LG1_40kb_TCGAluad_mutKRAS_mutEGFR_auc_obs.Rdata")))
auc_perm <- get(load(file.path("../FIGURES_V2_YUANLONG/FCC_WAVE_PLOT/LG1_40kb_TCGAluad_mutKRAS_mutEGFR_auc_permutQt.Rdata")))
stopifnot(round(auc_v0, 4) == round(auc_obs/auc_perm, 4))

inFolder <- "FCC_WAVE_PLOT_RESCALED"

pipFolder <- file.path( "../v2_Yuanlong_Cancer_HiC_data_TAD_DA/" , "PIPELINE", "OUTPUT_FOLDER")
all_hicds <- list.files(file.path(pipFolder))
all_hicds <- all_hicds[!(grepl("RANDOM", all_hicds) | grepl("PERMUT", all_hicds))]
all_exprds <- sapply(all_hicds, function(x) list.files(file.path(pipFolder, x)))

hicds="LG1_40kb"
exprds="TCGAluad_norm_luad"


if(buildTable){
  
  all_auc_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar%{
    exprds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      obsFile <- file.path(inFolder, paste0(hicds, "_", exprds, "_auc_obs.Rdata"))
      
      if(!file.exists(obsFile)) cat(obsFile, "\n")
      
      
      auc_obs <- get(load(obsFile))
      auc_permut <- get(load(file.path(inFolder, paste0(hicds, "_", exprds, "_auc_permutQt.Rdata"))))
      auc_ratio_rescaled <- auc_obs/auc_permut
      auc_ratio_v0 <- fcc_auc_dt$fcc_auc[fcc_auc_dt$hicds == hicds & fcc_auc_dt$exprds == exprds]
      stopifnot(length(auc_ratio_v0) == 1)
      
      data.frame(
        hicds = hicds,
        exprds = exprds,
        auc_ratio_rescaled = auc_ratio_rescaled,
        auc_ratio_v0 = auc_ratio_v0,
        stringsAsFactors = FALSE
      )
      
    }
    exprds_dt
  }
  
  outFile <- file.path(outFolder, "all_auc_dt.Rdata")
  save(all_auc_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "all_auc_dt.Rdata")
  all_auc_dt <- get(load(outFile))
}


all_auc_dt$cmp <- all_cols[all_cmps[paste0(all_auc_dt$exprds)]]

outFile <- file.path(outFolder, paste0("auc_ratio_rescaled_vs_V0.", plotType))
do.call(plotType, list(outFile, height=myHeightGG, width=myHeightGG))
plot(
 x =  all_auc_dt$auc_ratio_v0,
 y =  all_auc_dt$auc_ratio_rescaled,
 xlab  = "v0",
 ylab="rescaled",
 col = all_auc_dt$cmp,
 pch=16,
 main = paste0("Comp. FCC AUC ratio (PERMG2T)"),
 cex.axis=plotCex,
 cex.lab =plotCex,
 cex.main=plotCex
)

addCorr( x =  all_auc_dt$auc_ratio_v0,
         y =  all_auc_dt$auc_ratio_rescaled, bty="n", legPos = "bottomright")
# 
legend("topleft", col=all_cols, legend=names(all_cols), pch=16, bty="n")
















