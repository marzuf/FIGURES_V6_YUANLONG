# Rscript range_fcc_randommidposstrict_iterate.R

outFolder <- "RANGE_FCC_RANDOMMIDPOSSTRICT_ITERATE"
dir.create(outFolder)


plotType <- "svg"
myHeight <- 7
myWidth <- 10
plotCex <- 1.4

source("../Cancer_HiC_data_TAD_DA/utils_fct.R")
source("../Yuanlong_Cancer_HiC_data_TAD_DA/subtype_cols.R")
source("../FIGURES_V2_YUANLONG/settings.R")

require(foreach)
require(doMC)
require(ggplot2)
registerDoMC(40)
require(ggpubr)


buildData <- TRUE

all_fccThresh <- seq(from = -1, to = 1, by = 0.05)

fcc_fract <- seq(from=-1, to=1, by=0.25)
# fcc_fract_names <- paste0("FCC > ", fcc_fract[1:(length(fcc_fract)-1)], " and FCC <= ",fcc_fract[2:length(fcc_fract)])
fcc_fract_names <- paste0("FCC \u2208 ]", fcc_fract[1:(length(fcc_fract)-1)], ", ",fcc_fract[2:length(fcc_fract)], "]")
fcc_fract_names <- paste0("]", fcc_fract[1:(length(fcc_fract)-1)], ", ",fcc_fract[2:length(fcc_fract)], "]")
fcc_fract_names[fcc_fract_names == "]-1, -0.75]"] <- "[-1, -0.75]"
# fract_sort <- "FCC > 0.75 and FCC <= 1"
fract_sort <- fcc_fract_names[length(fcc_fract_names)]


rd_type <- "RANDOMMIDPOSSTRICT"

if(buildData){
  hicds = all_hicds[1]
  hicds = all_hicds[2]
  all_fccRange_dt <- foreach(hicds = all_hicds, .combine='rbind') %dopar% {
    exprds = all_exprds[[paste0(hicds)]][1]
    hicds_dt <- foreach(exprds = all_exprds[[paste0(hicds)]], .combine='rbind') %do% {
      
      cat("... start ", hicds, " - ", exprds, " \n")
      
      rd_file <- file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER", gsub("_40kb", "_RANDOMMIDPOSSTRICT_40kb", hicds), exprds, "8cOnlyFCC_runAllDown/all_obs_prodSignedRatio.Rdata")
      rd_fcc <- get(load(rd_file))
      
      rd_fcc_hist <- hist(rd_fcc, breaks=fcc_fract)$counts
      names(rd_fcc_hist) <- fcc_fract_names
      
      rd_fcc_hist <- rd_fcc_hist/length(rd_fcc)
      stopifnot(sum(rd_fcc_hist) == 1)
      
      fccFile <- file.path("../v2_Yuanlong_Cancer_HiC_data_TAD_DA/PIPELINE/OUTPUT_FOLDER", hicds, exprds, "8cOnlyFCC_runAllDown/all_obs_prodSignedRatio.Rdata")
      obs_fcc <- get(load(fccFile))
      
      obs_fcc_hist <- hist(obs_fcc, breaks=fcc_fract)$counts
      names(obs_fcc_hist) <- fcc_fract_names
      
      obs_fcc_hist <- obs_fcc_hist/length(obs_fcc)
      stopifnot(sum(obs_fcc_hist) == 1)
      
      
     dt1 <- data.frame(
        hicds=hicds,
        exprds=exprds,
        # fcc_type="random",
        intervalFCC = names(rd_fcc_hist),
        random_ratioFCC = as.numeric(rd_fcc_hist),
        stringsAsFactors = FALSE
      )
     dt2 <- data.frame(
       hicds=hicds,
       exprds=exprds,
       # fcc_type="observed",
       intervalFCC = names(obs_fcc_hist),
       obs_ratioFCC = as.numeric(obs_fcc_hist),
       stringsAsFactors = FALSE
     )
     merge(dt1, dt2, by=c("hicds", "exprds", "intervalFCC"), all=TRUE)
    }
    hicds_dt
  }
  outFile <- file.path(outFolder, "all_fccRange_dt.Rdata")
  save(all_fccRange_dt, file=outFile, version=2)
  cat(paste0("... written: ", outFile, "\n"))
} else {
  outFile <- file.path(outFolder, "all_fccRange_dt.Rdata")
  all_fccRange_dt <- get(load(outFile))
}

all_fccRange_dt$intervalFCC <- factor(all_fccRange_dt$intervalFCC, levels = fcc_fract_names)
stopifnot(!is.na(all_fccRange_dt$intervalFCC))


all_fccRange_dt$dataset <- file.path(all_fccRange_dt$hicds, all_fccRange_dt$exprds)

all_fccRange_dt$obs_over_rd <- all_fccRange_dt$obs_ratioFCC/all_fccRange_dt$random_ratioFCC



# load("N_FCC_OVER_QTILE_V2_EITHER_ITERATE/all_fccRange_dt.Rdata")

save_dt <- all_fccRange_dt



box_p <- ggboxplot(data = all_fccRange_dt,
                   x = "intervalFCC", y = "obs_over_rd",
                   xlab="FCC threshold", ylab="Ratio TADs FCC in range",
                   title = paste0(rd_type))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  geom_hline(yintercept=1, linetype=2, color="red")+
  theme(panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
        panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
        axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
        plot.title = element_text(hjust=0.5, size = 16, face="bold")
  )

outFile <- file.path(outFolder, paste0("ratioObsPermut_overThresh_boxplot.", plotType))
ggsave(box_p, file=outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


all_fccRange_dt$cmp <- all_cmps[paste0(all_fccRange_dt$exprds)]


box_p_cond <- ggboxplot(data = all_fccRange_dt,
                        x = "intervalFCC", y = "obs_over_rd", 
                        color = "cmp",fill = "cmp",
                        xlab="FCC threshold", ylab="Ratio TADs FCC in range",
                        title = paste0(rd_type))+
  scale_y_continuous(breaks = scales::pretty_breaks(n = 10))+
  geom_hline(yintercept=1, linetype=2, color="red")+
  labs(color="", fill="")+
  theme(panel.grid.major.y =  element_line(colour = "grey", size = 0.5, linetype=1),
        panel.grid.minor.y =  element_line(colour = "grey", size = 0.5, linetype=1),
        axis.text.x = element_text(size=10, angle=90, hjust=1, vjust=0.5),
        axis.title.x = element_text(size=14, hjust=0.5, vjust=0.5),
        plot.title = element_text(hjust=0.5, size = 16, face="bold")
  ) +
  scale_color_manual(values = all_cols) + 
  scale_fill_manual(values = all_cols) 

outFile <- file.path(outFolder, paste0("ratioObsPermut_overThresh_boxplot_cmpType.", plotType))
ggsave(box_p_cond, file=outFile, height=myHeightGG, width=myWidthGG)
cat(paste0("... written: ", outFile, "\n"))


























