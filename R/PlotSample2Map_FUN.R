#' Plot samples to map colour-coded as for cluster asignment
#' 
#' dirIn the directory path where the mcmc results are saved
#' map a map of class sf
#' pal palette for sample colours
#' UTM Whether the coordinats provided in the analysis to Geneland are in UTM (TRUE) 
#'     or GDS84 / MGA 55 (FALSE)
#' w, h the extent of the jitter for sample locations to be plotted. NULL (default) if none
#' a set alpha for transparency
PlotSample2Map <- function(dirIn, map=NULL, txt,
                            pal="Set1", UTM=TRUE, w=NULL, h=NULL, a=0.2) {
  library(ggplot2)
  library(maps)
  library(rgdal)
  library(data.table)
  library(sf)
  library(ggrepel)
  library(ggspatial)
  dat <- fread(file.path(dirIn, "modal.pop.indiv.txt"))
  setnames(dat, c("Long", "Lat", "Cluster"))
  
  if(UTM) {
    coord <- dat[, .(Long, Lat)]
    coordinates(coord) <- c("Long", "Lat" )
    proj4string(coord) <- CRS("+proj=utm +zone=55 ellps=GDA94")
    
    coord_dec <- spTransform(coord, CRS("+proj=utm +zone=55 +south +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"))
    coord_dec <-as.data.frame(coord_dec)
    dat[, Long:=coord_dec$Long]
    dat[, Lat:=coord_dec$Lat]
  }
  dat[, Cluster:=as.factor(Cluster)]
  nK <- length(dat[, levels(Cluster)])
  
  if(is.null(map)) map <- map_data("world", region=reg)
  
  mapplot1 <- ggplot(map) + 
    geom_sf(fill=NA) +
    annotation_scale(location = "br", width_hint = 0.25) +
    annotation_north_arrow(location = "br", which_north = "true", 
                           pad_x = unit(0.75, "in"), pad_y = unit(0.5, "in"),
                           style = north_arrow_fancy_orienteering) +
    # geom_point(data=dat, aes(Long, Lat), col="black", size=2) +
    geom_jitter(data=dat, aes(Long, Lat, col=Cluster), 
                size=1, width=w, height=h, alpha=a) +
    scale_color_brewer(palette=pal) +
    #  ggtitle(paste0(gsub(pattern="\\.|/", "",  dirIn), "_K", nK)) + 
    theme_bw() +
    labs(x="Longitude", y="Latitude") +
    geom_text_repel(aes(x=Long, y=Lat, group=Location, label=Location), #col="purple",
                    data=txt, stat="identity",
                    force=28,
    )#+
  #theme(legend.position="top") +
  #guides(fill=guide_legend(title=legend_title))
  
  ggsave(file.path(dirIn, paste0("MapSample_", basename(dirIn), "_K", nK, ".pdf")), 
         plot=mapplot1 + ggtitle(paste0(basename(dirIn), "_K", nK)), 
         width=24, height=19, units="cm", dpi=250)
  
  ggsave(file.path(dirIn, paste0("MapSample_", basename(dirIn), "_K", nK, ".png")), 
         plot=mapplot1, 
         width=25, height=12, units="cm", dpi=250)
  
  return(mapplot1)
}
