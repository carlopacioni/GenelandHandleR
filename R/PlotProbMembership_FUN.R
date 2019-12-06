#' Generate a bar plot with probability of membership for each sample
#' 
#' @param path.mcmc Path where the outputs are stored
#' @param pal Palette to use to colour code membership
#' @param matadata Optional: data.frame with nindiv lines with metadata associated with samples 
#'       as columns. One of these can be used to sort the samples
#' @param orderBy Optional: name of the column in metadata to order the samples by   
#' @import data.table
#' @import ggplot2
#' @export  
PlotProbMembership <- function(path.mcmc, pal="Set1", metadata, orderBy) {
  probInd <- fread(file.path(path.mcmc, "proba.pop.membership.indiv.txt"))
  keep <- probInd[, lapply(.SD, sum)]
  probInd <- probInd[, as.logical(keep[1,]), with=FALSE]
  nK <- ncol(probInd) - 2
  Knames <- paste0("C", seq_len(nK))
  setnames(probInd, c("X", "Y", Knames))
  
  probInd[, SampleOrd:=seq_len(nrow(probInd))]
  if(hasArg(orderBy)) {
    if(hasArg(metadata)) {
      probInd[, NewOrder := metadata[, orderBy]]
      probInd <- probInd[order(probInd[, NewOrder]),]
    } else {
      if(orderBy %in% names(probInd)) {
        probInd <- probInd[order(probInd[, get(orderBy)]),]
      } else {
        stop("Argument provided for 'orderBy', but no values passed to 'metadata")
      }
    }
  }
  probInd[, SampleOrd:=factor(SampleOrd, levels = unique(SampleOrd))]
  suppressWarnings(
    probInd_long <- melt.data.table(probInd, id.vars="SampleOrd", measure.vars=patterns("C[0-9]+"), 
                                    variable.name="Cluster", value.name="ProbMembership")
  )
  p <- ggplot(probInd_long, aes(SampleOrd, ProbMembership, fill=Cluster)) + 
    #geom_col() +
    geom_bar(stat="identity", width=1) +
    scale_fill_brewer(palette=pal) +
    theme_classic() +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          #axis.text.x=element_text(colour="black", angle=90), 
          legend.position="none") +
    scale_y_continuous(expand=c(0,0)) +
    ggtitle(basename(path.mcmc)) #+
  #scale_x_discrete(breaks=qind$Sample[nperPop$Bks], labels=nperPop$Label) +
  #geom_vline(xintercept=nperPop[, PopLim], size=0.1)
  
  write.csv(probInd, file.path(path.mcmc, 
                               paste0("ProbMemb_K", nK, basename(path.mcmc), ".csv")),row.names = FALSE)
  write.csv(probInd_long, file.path(path.mcmc, 
                                    paste0("ProbMemb_K", nK, basename(path.mcmc), "_long", ".csv")),row.names = FALSE)
  ggsave(file.path(path.mcmc, paste0("BarPlot_ProbMemb_K", nK, basename(path.mcmc), ".pdf")), 
         plot=p, width=24, height=8, units="cm", dpi=250)
  return(p)
}