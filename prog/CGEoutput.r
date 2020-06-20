# This program is AIM/Enduse and CGE Japan coupling analysis
#----------------------data loading and parameter settings ----------------------*
if(insflag==1){
options(CRAN="http://cran.md.tsukuba.ac.jp/")
install.packages("ggplot2", dependencies = TRUE)
install.packages("RColorBrewer", dependencies = TRUE)
#install.packages("grid", dependencies = TRUE)
#install.packages("gdxrrw", dependencies = TRUE)
install.packages("dplyr", dependencies = TRUE)
install.packages("sp", dependencies = TRUE)
install.packages("maptools", dependencies = TRUE)
install.packages("maps", dependencies = TRUE)
install.packages("ggradar", dependencies = TRUE)
install.packages("fmsb", dependencies = TRUE)
install.packages("tidyr", dependencies = TRUE)
install.packages("stringr", dependencies = TRUE)
install.packages("rJava", dependencies = TRUE)
install.packages("Rcpp", dependencies = TRUE)
install.packages("ReporteRsjars", dependencies = TRUE)
install.packages("ReporteRs", dependencies = TRUE)
install.packages("xlsx", dependencies = TRUE)
install.packages("R2PPT", dependencies = TRUE) #Rtools needs to be installed
install.packages('RDCOMClient', repos = 'http://www.omegahat.net/R/')
}

library(gdxrrw)
library(ggplot2)
library(dplyr)
library(reshape2)
library(tidyr)
library(maps)
library(grid)
library(RColorBrewer)
library(R2PPT)
library(RDCOMClient)

OrRdPal <- brewer.pal(9, "OrRd")
set2Pal <- brewer.pal(8, "Set2")
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
#file.copy("E:/sfujimori/CGE/AIMHub2.2ESIntAsia/AIMCGE/individual/IEAEB1062CGE/output/IEAEBIAMCTemplate.gdx", "../data/IEAEBIAMCTemplate.gdx",overwrite = TRUE)
filename <- c("global_17_emf.gdx")
file.copy("E:/sfujimori/CGE/AIMHub2.2ESIntAsia/anls_output/iiasa_database/gdx/global_17_emf.gdx", "../modeloutput/global_17_emf.gdx",overwrite = TRUE)
#file.copy("E:/sfujimori/CGE/AIMHub2.2ESIntAsia/anls_output/iiasa_database/gdx/JPN_emf.gdx", "../modeloutput/JPN_emf.gdx",overwrite = TRUE)
#file.copy("E:/sfujimori/CGE/AIMHub2.2ESIntAsia/anls_output/iiasa_database/gdx/CHN_emf.gdx", "../modeloutput/CHN_emf.gdx",overwrite = TRUE)
#file.copy("E:/sfujimori/CGE/AIMHub2.2ESIntAsia/anls_output/iiasa_database/gdx/IND_emf.gdx", "../modeloutput/IND_emf.gdx",overwrite = TRUE)
#file.copy("E:/sfujimori/CGE/AIMHub2.2ESIntAsia/AIMCGE/individual/AIMEnduseG2CGE/data/AIMEnduseG.gdx", "../modeloutput/AIMEnduseG.gdx",overwrite = TRUE)
#filename <- c("JPN_emf.gdx")
linepalette <- c("#4DAF4A","#FF7F00","#377EB8","#E41A1C","#984EA3","#F781BF","#8DD3C7","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#7f878f","#A65628","#FFFF33")
landusepalette <- c("#8DD3C7","#FF7F00","#377EB8","#4DAF4A","#A65628")
scenariomap <- read.table("../data/scenariomap.map", sep="\t",header=T, stringsAsFactors=F)
scenariomap2 <- read.table("../data/scenariomap2.map", sep="\t",header=T, stringsAsFactors=F)
region <- read.table("../data/region.txt", sep="\t",header=F, stringsAsFactors=F)
varlist_load <- read.table("../data/varlist.txt", sep="\t",header=F, stringsAsFactors=F)
varalllist <- read.table("../data/varalllist.txt", sep="\t",header=F, stringsAsFactors=F)
varlist <- left_join(varlist_load,varalllist,by="V1")
areapalette <- c("Coal|w/o CCS"="#000000","Coal|w/ CCS"="#7f878f","Oil|w/o CCS"="#ff2800","Oil|w/ CCS"="#ffd1d1","Gas|w/o CCS"="#9a0079","Gas|w/ CCS"="#c7b2de","Hydro"="#0041ff","Nuclear"="#663300","Solar"="#b4ebfa","Wind"="#ff9900","Biomass|w/o CCS"="#35a16b","Biomass|w/ CCS"="#cbf266","Geothermal"="#edc58f","Other"="#ffff99",
                 "Solid"=pastelpal[1],"Liquid"=pastelpal[2],"Gas"=pastelpal[3],"Electricity"=pastelpal[4],"Heat"=pastelpal[5],"Hydrogen"=pastelpal[6],
                 "Industry"=set2Pal[1],"Transport"=set2Pal[2],"Commercial"=set2Pal[3],"Residential"=set2Pal[4],
                 "Build-up"=pastelpal[1],"Cropland (for food)"=pastelpal[2],"Forest"=pastelpal[3],"Pasture"=pastelpal[4],"Energy Crops"=pastelpal[5],"Other Land"=pastelpal[6],"Other Arable Land"=pastelpal[7])
areamap <- read.table("../data/Areafigureorder.txt", sep="\t",header=T, stringsAsFactors=F)
areamappara <- read.table("../data/Area.map", sep="\t",header=T, stringsAsFactors=F)

#---IAMC tempalte loading and data merge
CGEload0 <- rgdx.param(paste0('../modeloutput/',filename),'EMFtemp1') 
#CGEload0 <- rbind(rgdx.param(paste0('../modeloutput/JPN_emf.gdx'),'EMFtemp1'), rgdx.param(paste0('../modeloutput/IND_emf.gdx'),'EMFtemp1'), rgdx.param(paste0('../modeloutput/CHN_emf.gdx'),'EMFtemp1'))
CGEload1 <- CGEload0 %>% rename("Value"=EMFtemp1,"Variable"=VEMF) %>% 
  left_join(scenariomap,by="SCENARIO") %>% filter(SCENARIO %in% as.vector(scenariomap[,1]) & REMF %in% (region$V1)) %>% 
  select(-SCENARIO) %>% rename(Region="REMF",SCENARIO="Name")

#Enduse loading
enduseflag <-0
if(enduseflag==1){
  EnduseJload0 <- rgdx.param(paste0('../modeloutput/AIMEnduse.gdx'),'EMFtemp1') %>% rename("SCENARIO"=i1,"Region"=i2,"Variable"=i3,"Y"=i4,"Value"=value) %>% mutate(Model="AIM/Enduse[Japan]")
  EnduseJload1 <- EnduseJload0 %>% left_join(scenariomap2,by="SCENARIO") %>% filter(SCENARIO %in% as.vector(scenariomap2[,1]) & Region %in% region$V1) %>% 
    select(-SCENARIO) %>% rename(SCENARIO="Name")

  EnduseGload0 <- rgdx.param(paste0('../modeloutput/AIMEnduseG.gdx'),'IAMC_template') %>% select(-i1) %>% rename("SCENARIO"=i2,"Region"=i3,"Variable"=i4,"Unit"=i5,"Y"=i6,"Value"=value)  %>% mutate(Model="AIM/Enduse[Global]")
  EnduseGload1 <- EnduseGload0 %>% left_join(scenariomap2,by="SCENARIO") %>% filter(SCENARIO %in% as.vector(scenariomap2[,1]) & Region %in% region$V1) %>% 
    select(-SCENARIO,-Unit) %>% rename(SCENARIO="Name")
}

IEAEB0 <- rgdx.param('../data/IEAEBIAMCTemplate.gdx','IAMCtemp17') %>% rename("Value"=IAMCtemp17,"Variable"=VEMF,"Y"=St,"Region"=Sr17,"SCENARIO"=SceEneMod) %>%
  select(Region,Variable,Y,Value,SCENARIO) %>% filter(Region %in% region$V1) %>% mutate(Model="Reference")
IEAEB0$Y <- as.numeric(levels(IEAEB0$Y))[IEAEB0$Y]
IEAEB1 <- filter(IEAEB0,Y<=2015 & Y>=1990)

#allmodel0 <- rbind(CGEload1,EnduseGload1,EnduseJload1)  
allmodel0 <- rbind(CGEload1)  
allmodel0$Y <- as.numeric(levels(allmodel0$Y))[allmodel0$Y]

allmodel <- rbind(allmodel0,IEAEB1)  

#---function
plot.1 <- function(XX){
  plot <- ggplot() + 
    geom_area(data=XX,aes(x=Y, y = Value , fill=reorder(Ind,-order)), stat="identity") + 
    ylab(ylab1) + xlab(xlab1) +labs(fill="")+ guides(fill=guide_legend(reverse=TRUE)) + MyThemeLine +
    theme(legend.position="bottom", text=element_text(size=12),  
          axis.text.x=element_text(angle=45, vjust=0.9, hjust=1, size = 12)) +
    guides(fill=guide_legend(ncol=5)) + scale_x_continuous(breaks=seq(miny,maxy,10)) +  ggtitle(paste(rr,areamappara$Class[j],sep=" "))
  
  plot2 <- plot +facet_wrap(Model ~ SCENARIO,ncol=4) + scale_fill_manual(values=colorpal) + 
#  plot2 <- plot +facet_grid(Model~SCENARIO) + scale_fill_manual(values=colorpal) + 
    annotate("segment",x=miny,xend=maxy,y=0,yend=0,linetype="solid",color="grey") + theme(legend.position='bottom')
  if(nrow(XX2)>=1){
    plot3 <- plot2 +    geom_area(data=XX2,aes(x=Y, y = Value , fill=reorder(Ind,-order)), stat="identity")
  }else{
    plot3 <- plot2
  }
  return(plot3)
}

#---IAMC tempalte loading and data mergeEnd
for(rr in region$V1){
  dir.create(paste0("../output/",rr))
  dir.create(paste0("../output/",rr,"/png"))
  dir.create(paste0("../output/",rr,"/ppt"))

nalist <- c(as.vector(varlist$V1),"TPES","POWER","Landuse","TFC_fuel","TFC_Sector","TFC_Ind","TFC_Tra","TFC_Res","TFC_Com")
allplot <- as.list(nalist)
plotflag <- as.list(nalist)
names(allplot) <- nalist
names(plotflag) <- nalist

maxy <- max(allmodel$Y)

#---Line figures
for (i in 1:nrow(varlist)){
  if(nrow(filter(allmodel,Variable==varlist[i,1] & Region==rr))>0){
    miny <- min(filter(allmodel,Variable==varlist[i,1] & Region==rr)$Y) 
    plot.0 <- ggplot() + 
      geom_line(data=filter(allmodel,Variable==varlist[i,1] & Model!="Reference"& Region==rr),aes(x=Y, y = Value , color=interaction(SCENARIO,Model),group=interaction(SCENARIO,Model)),stat="identity") +
      geom_point(data=filter(allmodel,Variable==varlist[i,1] & Model!="Reference"& Region==rr),aes(x=Y, y = Value , color=interaction(SCENARIO,Model),shape=Model),size=3.0,fill="white") +
      MyThemeLine + scale_color_manual(values=linepalette) + scale_x_continuous(breaks=seq(miny,maxy,10)) +
      xlab("year") + ylab(varlist[i,3])  +  ggtitle(paste(rr,varlist[i,2],sep=" ")) +
      annotate("segment",x=2005,xend=maxy,y=0,yend=0,linetype="dashed",color="grey")+ 
      theme(legend.title=element_blank()) 
    if(length(scenariomap$SCENARIO)<20){
      plot.0 <- plot.0 +
      geom_point(data=filter(allmodel,Variable==varlist[i,1] & Model=="Reference"& Region==rr),aes(x=Y, y = Value) , color="black",shape=6,size=2.0,fill="grey") 
    }
    outname <- paste0(outputdir,rr,"/png/",varlist[i,1],".png")
    ggsave(plot.0, file=outname, dpi = 150, width=10, height=6,limitsize=FALSE)
    allplot[[nalist[i]]] <- plot.0
  }
  plotflag[[nalist[i]]] <- nrow(filter(allmodel,Variable==varlist[i,1] & Model!="Reference"& Region==rr))
}
#---Area figures

for(j in 1:nrow(areamappara)){
  XX <- allmodel %>% filter(Variable %in% as.vector(areamap$Variable)) %>% left_join(areamap,by="Variable") %>% ungroup() %>% 
    filter(Class==areamappara[j,1] & Model!="Reference"& Region==rr) %>% select(Model,SCENARIO,Ind,Y,Value,order)  %>% arrange(order)
  XX2 <- allmodel %>% filter(Variable %in% as.vector(areamap$Variable)) %>% left_join(areamap,by="Variable") %>% ungroup() %>% 
    filter(Class==areamappara[j,1] & Model=="Reference"& Region==rr) %>% select(-SCENARIO,-Model,Ind,Y,Value,order)  %>% arrange(order)
  miny <- min(XX$Y,XX2$Y) 
  na.omit(XX$Value)
  unit_name <-areamappara[j,3]
  ylab1 <- paste0(areamappara[j,2], " (", unit_name, ")")
  xlab1 <- areamappara[j,2]
  colorpal <- areapalette
  plot_TPES.1 <- plot.1(XX)
  allplot[[areamappara$Class[j]]] <- plot_TPES.1 
  outname <- paste0(outputdir,rr,"/png/",areamappara[j,1],".png")
  ggsave(plot_TPES.1, file=outname, dpi = 450, width=9, height=floor(length(unique(XX$SCENARIO))/4)*3,limitsize=FALSE)
  plotflag[[areamappara$Class[j]]] <- nrow(XX)  
}

#----r2ppt
#The figure should be prearranged before going this ppt process since emf file type does not accept size changes. 
#If you really needs ppt slide, you first ouptput png and then paste it.

myPPT<-PPT.Init(method="RDCOMClient")
for (i in 1:length(nalist)){
    if(plotflag[[nalist[i]]]>0){
#      win.graph(width=1860, height=1450,pointsize = 1)
#      print(allplot[[i]])
      myPPT<-PPT.AddTitleOnlySlide(myPPT,title="Title Only",title.fontsize=40,title.font="Arial")
#      myPPT<-PPT.AddGraphicstoSlide(myPPT,size= c(10,10,700,350), dev.out.type ='emf' )
      myPPT<-PPT.AddGraphicstoSlide(myPPT,file=paste0(outputdir,rr,"/png/",nalist[i],".png"),size=c(10,10,700,500))
#      dev.off()
  }
}
myPPT<-PPT.SaveAs(myPPT,file=paste0("../output/",rr,"/ppt/",rr,"comparison.pptx"))
myPPT<-PPT.Close(myPPT)
rm(myPPT)
#      savePlot("test.emf",type="emf", device = dev.cur())
#      myPPT<-PPT.AddGraphicstoSlide(myPPT,file="test.emf",dev.out.type="emf",size=c(10,10,500,350))

}

