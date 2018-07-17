# This program is AIM/Enduse and CGE Japan coupling analysis

library(gdxrrw)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(maps)
library(grid)
library(RColorBrewer)
OrRdPal <- brewer.pal(9, "OrRd")
YlGnBupal <- brewer.pal(9, "YlGnBu")
Redspal <- brewer.pal(9, "Reds")
pastelpal <- brewer.pal(9, "Pastel1")
pastelpal <- brewer.pal(8, "Set1")

MyThemeLine <- theme_bw() +
  theme(
    panel.border=element_rect(fill=NA),
    panel.grid.minor = element_line(color = NA), 
    #    axis.title=element_text(size=5),
    #    axis.text.x = element_text(hjust=1,size = 10, angle = 0),
    axis.line=element_line(colour="black"),
    panel.background=element_rect(fill = "white"),
    #    panel.grid.major=element_line(linetype="dashed",colour="grey",size=0.5),
    panel.grid.major=element_blank(),
    strip.background=element_rect(fill="white", colour="white"),
    strip.text.x = element_text(size=10, colour = "black", angle = 0,face="bold"),
    axis.text.x=element_text(size = 10,angle=45, vjust=0.9, hjust=1, margin = unit(c(t = 0.3, r = 0, b = 0, l = 0), "cm")),
    axis.text.y=element_text(size = 10,margin = unit(c(t = 0, r = 0.3, b = 0, l = 0), "cm")),
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10),
    axis.ticks.length=unit(-0.15,"cm")
  )


#-- data load
dir.create("../output/")
outputdir <- c("../output/")
filename <- c("global_17_emf.gdx")
linepalette <- c("#4DAF4A","#FF7F00","#377EB8","#E41A1C","#984EA3","#FFFF33","#A65628","#F781BF","#8DD3C7","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#7f878f")
landusepalette <- c("#8DD3C7","#FF7F00","#377EB8","#4DAF4A","#A65628")
scenariomap <- read.table("../data/scenariomap.map", sep="\t",header=T, stringsAsFactors=F)
region <- read.table("../data/region.txt", sep="\t",header=F, stringsAsFactors=F)
varlist_load <- read.table("../data/varlist.txt", sep="\t",header=F, stringsAsFactors=F)
varalllist <- read.table("../data/varalllist.txt", sep="\t",header=F, stringsAsFactors=F)

varlist <- left_join(varlist_load,varalllist,by="V1")

CGEload0 <- rgdx.param(paste0('../data/',filename),'EMFtemp1') %>% rename("Value"=EMFtemp1,"Variable"=VEMF) 
CGEload1 <- CGEload0 %>% left_join(scenariomap,by="SCENARIO") %>% filter(SCENARIO %in% as.vector(scenariomap[,1]) & REMF %in% region) %>% 
   select(-SCENARIO) %>% rename(Region="REMF",SCENARIO="Name")
allmodel <- CGEload1  
allmodel$Y <- as.numeric(levels(allmodel$Y))[allmodel$Y]

nalist <- c(varlist,"TPES","POWER","Landuse")
allplot <- as.list(nalist)
plotflag <- as.list(nalist)

#for (i in 1:nrow(varlist)){
for (i in 1:nrow(varlist)){
  if(length(allmodel[allmodel$Variable==varlist[i,1],c(1)])>0){
    plot.0 <- ggplot() + 
      geom_line(data=allmodel[allmodel$Variable==varlist[i,1] ,c(-2)],aes(x=Y, y = Value , color=SCENARIO),stat="identity") +
      geom_point(data=allmodel[allmodel$Variable==varlist[i,1] ,c(-2)],aes(x=Y, y = Value , color=SCENARIO),shape=1,size=3.0,fill="white") +
      MyThemeLine+ scale_color_manual(values=linepalette)  +
      xlab("year") + ylab(varlist[i,3])  +  ggtitle(varlist[i,2]) +
      annotate("segment",x=2005,xend=2100,y=0,yend=0,linetype="dashed",color="grey")+ 
      theme(legend.title=element_blank()) 
    outname <- paste0(outputdir,"",varlist[i,1],".png")
    ggsave(plot.0, file=outname, dpi = 150, width=7, height=4,limitsize=FALSE)
#  allplot[[i]] <- plot.1
  }
  plotflag[[i]] <- nrow(allmodel[allmodel$Variable==varlist[i,1] ,c(1)])
}

tpespalette <- c("Coal|w/o CCS"="#000000","Coal|w/ CCS"="#7f878f","Oil|w/o CCS"="#ff2800","Oil|w/ CCS"="#ffd1d1","Gas|w/o CCS"="#9a0079","Gas|w/ CCS"="#c7b2de","Hydro"="#0041ff","Nuclear"="#663300","Solar"="#b4ebfa","Wind"="#ff9900","Biomass|w/o CCS"="#35a16b","Biomass|w/ CCS"="#cbf266","Geothermal"="#edc58f","Other"="#ffff99",
                 "Solid"=pastelpal[1],"Liquid"=pastelpal[2],"Gas"=pastelpal[3],"Electricity"=pastelpal[4],"Heat"=pastelpal[5],
                 "Build-up"=pastelpal[1],"Cropland (for food)"=pastelpal[2],"Forest"=pastelpal[3],"Pasture"=pastelpal[4],"Energy Crops"=pastelpal[5],"Other Land"=pastelpal[6],"Other Arable Land"=pastelpal[7])
areamap <- read.table("../data/Areafigureorder.txt", sep="\t",header=T, stringsAsFactors=F)
areamappara <- read.table("../data/Area.map", sep="\t",header=T, stringsAsFactors=F)
area.0 <- allmodel %>% filter(Variable %in% as.vector(areamap$Variable)) %>% left_join(areamap,by="Variable") %>% ungroup()
plot.1 <- function(){
  plot <- ggplot() + geom_area(data=XX,aes(x=Y, y = Value , fill=reorder(Ind,-order)), stat="identity") + 
    ylab(ylab1) + xlab(xlab1) +labs(fill="")+ guides(fill=guide_legend(reverse=TRUE)) + MyThemeLine +
    theme(legend.position="bottom", text=element_text(size=12),  
          axis.text.x=element_text(angle=0, vjust=0.9, hjust=1, size = 12)) +
    guides(fill=guide_legend(ncol=5))
  plot2 <- plot +facet_wrap(~ SCENARIO) + scale_fill_manual(values=colorpal) + 
    annotate("segment",x=2005,xend=2100,y=0,yend=0,linetype="solid",color="grey") + theme(legend.position='bottom')
  return(plot2)
}

for(j in 1:nrow(areamappara)){
  XX <- area.0 %>% filter(Class==areamappara[j,1]) %>% select(SCENARIO,Ind,Y,Value,order)  %>% arrange(order)
  na.omit(XX$Value)
  unit_name <-areamappara[j,3]
  ylab1 <- paste0(areamappara[j,2], " (", unit_name, ")")
  xlab1 <- areamappara[j,2]
  colorpal <- tpespalette
  plot_TPES.1 <- plot.1()
  outname <- paste0(outputdir,areamappara[j,1],".png")
  ggsave(plot_TPES.1, file=outname, dpi = 450, width=9, height=6,limitsize=FALSE)
  
}

