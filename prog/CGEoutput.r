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
#install.packages("R2PPT", dependencies = TRUE) #Rtools needs to be installed
install.packages("officer", dependencies = TRUE) #Rtools needs to be installed
install.packages('RDCOMClient', repos = 'http://www.omegahat.net/R/')
library(devtools)
devtools::install_github("tomwenseleers/export")
install.packages("furrr", dependencies = TRUE)
install.packages("purrr", dependencies = TRUE)
install.packages("progressr", dependencies = TRUE)
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
library(cowplot)
library(officer)
library(export)
library(purrr)
library(furrr)
library(progressr)

#---------------switches to specify the run condition -----
filename <- "global_17" # filename should be "global_17","CHN","JPN"....
enduseflag <- 0   # If you would like to display AIM/Enduse outputs, make this parameter 1 otherwise 0.
dirCGEoutput <-"../../anls_output/iiasa_database/gdx/"  # directory where the CGE output is located 
parallelmode <- 1 #Switch for parallel process. if you would like to use multi-processors assign 1 otherwise 0.
#parallelmode <- 0 #Switch for parallel process. if you would like to use multi-processors assign 1 otherwise 0.
threadsnum <- min(floor(availableCores()/2),24)
r2ppt <- 0 #Switch for ppt export. if you would like to export as ppt then assign 1 otherwise 0.

#---------------End of switches to specify the run condition -----

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
dir.create("../output/ppt")
outputdir <- c("../output/")
#Linefigure cross country
dir.create(paste0("../output/","merge"))
dir.create(paste0("../output/","merge","/png"))
dir.create(paste0("../output/","merge","/pngdet"))

file.copy(paste0(dirCGEoutput,filename,"_IAMC.gdx"), paste0("../modeloutput/",filename,"_IAMC.gdx"),overwrite = TRUE)
file.copy(paste0(dirCGEoutput,"../../../AIMCGE/individual/AIMEnduseG2CGE/data/merged_output.gdx"), paste0("../modeloutput/AIMEnduseG.gdx"),overwrite = TRUE)
file.copy(paste0(dirCGEoutput,"../../../AIMCGE/individual/IEAEB1062CGE/output/IEAEBIAMCTemplate.gdx"), paste0("../data/IEAEBIAMCTemplate.gdx"),overwrite = TRUE)

linepalette <- c("#4DAF4A","#FF7F00","#377EB8","#E41A1C","#984EA3","#F781BF","#8DD3C7","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#7f878f","#A65628","#FFFF33","black")
#linepalette <- c("Baseline"="#4DAF4A","GlobalOptimalZero"="#FF7F00","NDC+Zero"="#377EB8","#E41A1C","#984EA3","#F781BF","#8DD3C7","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#7f878f","#A65628","#FFFF33")

landusepalette <- c("#8DD3C7","#FF7F00","#377EB8","#4DAF4A","#A65628")
scenariomap <- read.table("../data/scenariomap.map", sep="\t",header=T, stringsAsFactors=F)
scenariomap2 <- read.table("../data/scenariomap2.map", sep="\t",header=T, stringsAsFactors=F)
region_load <- as.vector(read.table("../data/region.txt", sep="\t",header=F, stringsAsFactors=F)$V1)
region <- region_load
varlist_load <- read.table("../data/varlist.txt", sep="\t",header=F, stringsAsFactors=F)
varalllist <- read.table("../data/varalllist.txt", sep="\t",header=F, stringsAsFactors=F)
varlist <- left_join(varlist_load,varalllist,by="V1")
areapaletteload <- read.table("../data/color.txt", sep="\t",header=F, stringsAsFactors=F)
areamap <- read.table("../data/Areafigureorder.txt", sep="\t",header=T, stringsAsFactors=F)
areamappara <- read.table("../data/Area.map", sep="\t",header=T, stringsAsFactors=F)

#---IAMC tempalte loading and data merge
CGEload0 <- rgdx.param(paste0('../modeloutput/',filename,"_IAMC.gdx"),'IAMC_Template') 
Getregion <- as.vector(unique(CGEload0$REMF))
if(length(Getregion)==1){region <- Getregion}
CGEload1 <- CGEload0 %>% rename("Value"=IAMC_Template,"Variable"=VEMF) %>% 
  left_join(scenariomap,by="SCENARIO") %>% filter(SCENARIO %in% as.vector(scenariomap[,1]) & REMF %in% region) %>% 
  select(-SCENARIO) %>% rename(Region="REMF",SCENARIO="Name",Y="YEMF")

#Enduse loading
if(enduseflag==1){
  EnduseGload0 <- rgdx.param(paste0('../modeloutput/AIMEnduseG.gdx'),'data_all')  %>% rename("SCENARIO"=Sc,"Region"=Sr,"Variable"=Sv,"Y"=Sy,"Value"=data_all)  %>% mutate(Model="AIM/Enduse[Global]")
  EnduseGload1 <- EnduseGload0 %>% left_join(scenariomap2,by="SCENARIO") %>% filter(SCENARIO %in% as.vector(scenariomap2[,1]) & Region %in% region) %>% 
    select(-SCENARIO) %>% rename(SCENARIO="Name") %>% select(Region,Variable,Y,Value,SCENARIO,Model)
}

#IEA energy balance information
IEAEB0 <- rgdx.param('../data/IEAEBIAMCTemplate.gdx','IAMCtemp17') %>% rename("Value"=IAMCtemp17,"Variable"=VEMF,"Y"=St,"Region"=Sr17,"SCENARIO"=SceEneMod) %>%
  select(Region,Variable,Y,Value,SCENARIO) %>% filter(Region %in% region) %>% mutate(Model="Reference")
IEAEB0$Y <- as.numeric(levels(IEAEB0$Y))[IEAEB0$Y]
IEAEB1 <- filter(IEAEB0,Y<=2015 & Y>=1990)

#allmodel0 <- rbind(CGEload1,EnduseGload1,EnduseJload1)  
if(enduseflag==1){
  allmodel0 <- rbind(CGEload1,EnduseGload1)  
}else{
  allmodel0 <- rbind(CGEload1)  
}
allmodel0$Y <- as.numeric(levels(allmodel0$Y))[allmodel0$Y]
allmodel <- rbind(allmodel0,IEAEB1)  
maxy <- max(allmodel$Y)
linepalettewName <- linepalette
names(linepalettewName) <- unique(allmodel$SCENARIO)

#---End of IAMC tempalte loading and data merge

#---functions
#function for regional figure generation
funcplotgen <- function(rr,progr){
  progr(message='region figures')
#  for(rr in region){
  dir.create(paste0("../output/",rr))
  dir.create(paste0("../output/",rr,"/png"))
  dir.create(paste0("../output/",rr,"/pngdet"))
  dir.create(paste0("../output/",rr,"/ppt"))
  dir.create(paste0("../output/",rr,"/merge"))
  
#---Line figures
  for (i in 1:nrow(varlist)){
    if(nrow(filter(allmodel,Variable==varlist[i,1] & Region==rr & Model!="Reference"))>0){
      miny <- min(filter(allmodel,Variable==varlist[i,1] & Region==rr)$Y) 
      plot.0 <- ggplot() + 
        geom_line(data=filter(allmodel,Variable==varlist[i,1] & Model!="Reference"& Region==rr),aes(x=Y, y = Value , color=SCENARIO,group=interaction(SCENARIO,Model)),stat="identity") +
        geom_point(data=filter(allmodel,Variable==varlist[i,1] & Model!="Reference"& Region==rr),aes(x=Y, y = Value , color=SCENARIO,shape=Model),size=3.0,fill="white") +
        MyThemeLine + scale_color_manual(values=linepalettewName) + scale_x_continuous(breaks=seq(miny,maxy,10)) +
        xlab("year") + ylab(varlist[i,4])  +  ggtitle(paste(rr,varlist[i,3],sep=" ")) +
        annotate("segment",x=miny,xend=maxy,y=0,yend=0,linetype="dashed",color="grey")+ 
        theme(legend.title=element_blank()) 
      if(length(scenariomap$SCENARIO)<20){
        plot.0 <- plot.0 +
        geom_point(data=filter(allmodel,Variable==varlist[i,1] & Model=="Reference"& Region==rr),aes(x=Y, y = Value) , color="black",shape=6,size=2.0,fill="grey") 
      }
      if(varlist[i,2]==1){
        outname <- paste0(outputdir,rr,"/png/",varlist[i,1],".png")
      }else{
        outname <- paste0(outputdir,rr,"/pngdet/",varlist[i,1],".png")
      }
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
      filter(Class==areamappara[j,1] & Model=="Reference"& Region==rr) %>% select(-SCENARIO,-Model,Ind,Y,Value,order)  %>% arrange(order)%>%
      filter(Y>=2015)
    miny <- min(XX$Y,XX2$Y) 
    na.omit(XX$Value)
    
    unit_name <-areamappara[j,3] 
    ylab1 <- paste0(areamappara[j,2], " (", unit_name, ")")
    xlab1 <- areamappara[j,2]
    
    areapaletteArea <- filter(areapaletteload,V1 %in% unique(XX$Ind))$V2
    names(areapaletteArea) <- filter(areapaletteload,V1 %in% unique(XX$Ind))$V1
    colorpal <- areapaletteArea 
    
    plot2 <- ggplot() + 
      geom_area(data=XX,aes(x=Y, y = Value , fill=reorder(Ind,-order)), stat="identity") + 
      ylab(ylab1) + xlab(xlab1) +labs(fill="")+ guides(fill=guide_legend(reverse=TRUE)) + MyThemeLine +
      theme(legend.position="bottom", text=element_text(size=12),  
            axis.text.x=element_text(angle=45, vjust=0.9, hjust=1, size = 12)) +
      guides(fill=guide_legend(ncol=5)) + scale_x_continuous(breaks=seq(miny,maxy,10)) +  ggtitle(paste(rr,areamappara$Class[j],sep=" "))+
      facet_wrap(Model ~ SCENARIO,ncol=4) + scale_fill_manual(values=colorpal) + 
      annotate("segment",x=miny,xend=maxy,y=0,yend=0,linetype="solid",color="grey") + theme(legend.position='bottom')
    if(nrow(XX2)>=1){
      plot3 <- plot2 +    geom_area(data=XX2,aes(x=Y, y = Value , fill=reorder(Ind,-order)), stat="identity")
    }else{
      plot3 <- plot2
    }
    allplot[[areamappara$Class[j]]] <- plot3 
    outname <- paste0(outputdir,rr,"/png/",areamappara[j,1],".png")
    ggsave(plot3, file=outname, dpi = 450, width=9, height=floor(length(unique(XX$SCENARIO))/4+1)*4+2,limitsize=FALSE)
    plotflag[[areamappara$Class[j]]] <- nrow(XX)  
  }

#---merged figures
  pp_tfc <- plot_grid(allplot[["TFC_Ind"]],allplot[["TFC_Tra"]],allplot[["TFC_Res"]],allplot[["TFC_Com"]],ncol=2,align = "hv")
  ggsave(pp_tfc, file=paste0(outputdir,rr,"/merge/tfc.png"), width=9*2, height=(floor(length(unique(allmodel$SCENARIO))/4+1)*3+2)*2,limitsize=FALSE)
  p_legend1 <- gtable::gtable_filter(ggplotGrob(allplot[["Fin_Ene"]]), pattern = "guide-box")
  pp_tfcind <- plot_grid(allplot[["Fin_Ene"]] + theme(legend.position="none"),allplot[["Fin_Ene_Ind"]] + theme(legend.position="none"),allplot[["Fin_Ene_Tra"]] + theme(legend.position="none"),allplot[["Fin_Ene_Res"]] + theme(legend.position="none"),allplot[["Fin_Ene_Com"]] + theme(legend.position="none"),
                         allplot[["Fin_Ene_Ele_Heat"]] + theme(legend.position="none"),allplot[["Fin_Ene_Gas"]] + theme(legend.position="none"),allplot[["Fin_Ene_Liq"]] + theme(legend.position="none"),allplot[["Fin_Ene_SolidsCoa"]] + theme(legend.position="none"),allplot[["Fin_Ene_SolidsBio"]] + theme(legend.position="none"),
                         allplot[["Fin_Ene_Ind_Ele_Heat"]] + theme(legend.position="none"),allplot[["Fin_Ene_Ind_Gas"]] + theme(legend.position="none"),allplot[["Fin_Ene_Ind_Liq"]] + theme(legend.position="none"),allplot[["Fin_Ene_Ind_SolidsCoa"]] + theme(legend.position="none"),allplot[["Fin_Ene_Ind_SolidsBio"]] + theme(legend.position="none"),
                         allplot[["Fin_Ene_Com_Ele_Heat"]] + theme(legend.position="none"),allplot[["Fin_Ene_Com_Gas"]] + theme(legend.position="none"),allplot[["Fin_Ene_Com_Liq"]] + theme(legend.position="none"),allplot[["Fin_Ene_Com_SolidsCoa"]] + theme(legend.position="none"),allplot[["Fin_Ene_Com_SolidsBio"]] + theme(legend.position="none"),
                         allplot[["Fin_Ene_Res_Ele_Heat"]] + theme(legend.position="none"),allplot[["Fin_Ene_Res_Gas"]] + theme(legend.position="none"),allplot[["Fin_Ene_Res_Liq"]] + theme(legend.position="none"),allplot[["Fin_Ene_Res_SolidsCoa"]] + theme(legend.position="none"),allplot[["Fin_Ene_Res_SolidsBio"]] + theme(legend.position="none"),
                         allplot[["Fin_Ene_Tra_Ele"]] + theme(legend.position="none"),allplot[["Fin_Ene_Tra_Liq_Bio"]] + theme(legend.position="none"),allplot[["Fin_Ene_Tra_Liq_Oil"]] + theme(legend.position="none"),NULL,p_legend1,
                         nrow=6,rel_widths =c(1,1,1,1,1),align = "hv")
  ggsave(pp_tfcind, file=paste0(outputdir,rr,"/merge/tfcind.png"), width=25, height=(floor(length(unique(allmodel$SCENARIO))/4+1)*4+2)*2,limitsize=FALSE)
  pp_area <- plot_grid(allplot[["TPES"]],allplot[["Power_heat"]],allplot[["Landuse"]],ncol=1,align = "hv")
  ggsave(pp_area, file=paste0(outputdir,rr,"/merge/majorArea.png"), width=15, height=(floor(length(unique(allmodel$SCENARIO))/4+1)*4+2)*4,limitsize=FALSE)
  pp_main <- plot_grid(allplot[["GDP_MER"]] + theme(legend.position="none"),allplot[["POP"]] + theme(legend.position="none"),allplot[["Tem_Glo_Mea"]],
                       allplot[["Emi_CO2_Ene_and_Ind_Pro"]] + theme(legend.position="none"),allplot[["Emi_CO2"]] + theme(legend.position="none"),allplot[["Emi_Kyo_Gas"]],
                      allplot[["Pol_Cos_GDP_Los_rat"]] + theme(legend.position="none"),allplot[["Pol_Cos_Cns_Los_rat"]] + theme(legend.position="none"),allplot[["Prc_Car"]],
                      allplot[["Pop_Ris_of_Hun"]] + theme(legend.position="none"),allplot[["Prc_Prm_Ene_Oil"]] + theme(legend.position="none"),allplot[["Prc_Sec_Ene_Ele"]],
                         nrow=4,rel_widths =c(1,1,1.5),align = "hv")
  ggsave(pp_main, file=paste0(outputdir,rr,"/merge/main.png"), width=15, height=15,limitsize=FALSE)

#----r2ppt
#The figure should be prearranged before going this ppt process since emf file type does not accept size changes. 
#If you really needs ppt slide, you first ouptput png and then paste it.
  pptlist <- c("Fin_Ene","Fin_Ene_Ele_Heat","Fin_Ene_Gas","Fin_Ene_Liq","Fin_Ene_Solids","Fin_Ene_Res","Fin_Ene_Com","Fin_Ene_Tra","Fin_Ene_Ind","Emi_CO2_Ene_and_Ind_Pro","Pol_Cos_GDP_Los_rat","Prc_Car","TPES","Power_heat")
  TorF <- 0
  if (r2ppt==1){
    for (i in pptlist){
        if(plotflag[[i]]>0){
          TorF = TorF + 1
          if(TorF>1){
            graph2ppt(allplot[[i]], file = paste0("../output/",rr,"/ppt/",rr,"comparison.pptx"),width = 10, height = 10, append = TRUE)
          }else{
            graph2ppt(allplot[[i]], file = paste0("../output/",rr,"/ppt/",rr,"comparison.pptx"),width = 10, height = 10, append = FALSE)
          }
        }
    }
  }
}

# making cross regional figure
mergefigGen <- function(ii,progr){
  progr(message='merge figures')
#    for(ii in lst$varlist){
  if(nrow(filter(allmodel,Variable==ii  & Model!="Reference"))>0){
    miny <- 2010 
    plot.0 <- ggplot() + 
      geom_line(data=filter(allmodel,Variable==ii & Model!="Reference"),aes(x=Y, y = Value , color=SCENARIO,group=interaction(SCENARIO,Model)),stat="identity") +
      geom_point(data=filter(allmodel,Variable==ii & Model!="Reference"),aes(x=Y, y = Value , color=SCENARIO,shape=Model),size=3.0,fill="white") +
      MyThemeLine + scale_color_manual(values=linepalettewName) + scale_x_continuous(breaks=seq(miny,maxy,10)) +
      xlab("year") + ylab(varlist$V3[varlist$V1==ii])  +  ggtitle(paste("Multi-regions",varlist$V2.y[varlist$V1==ii],sep=" ")) +
      annotate("segment",x=miny,xend=maxy,y=0,yend=0,linetype="dashed",color="grey")+ 
      theme(legend.title=element_blank()) +facet_wrap(~Region,scales="free")
    if(length(scenariomap$SCENARIO)<20){
      plot.0 <- plot.0 +
        geom_point(data=filter(allmodel,Variable==ii & Model=="Reference"),aes(x=Y, y = Value) , color="black",shape=6,size=2.0,fill="grey") 
    }
    if(varlist$V2.x[varlist$V1==ii]==1){
      outname <- paste0(outputdir,"merge","/png/",ii,".png")
    }else{
      outname <- paste0(outputdir,"merge","/pngdet/",ii,".png")
    }
    ggsave(plot.0, file=outname, dpi = 150, width=15, height=12,limitsize=FALSE)
  }
}

#execute making regional figures
exe_fig_make <- function(ListIte,Xfunc){
  print(Sys.time())
  if(parallelmode==1){
    on.exit(plan(oplan), add = TRUE)
    oplan <- plan(multisession,workers=threadsnum)
    oplan
    handlers('progress')
    with_progress({
      progr <- progressor(along=ListIte)
      ListIte %>% future_map(Xfunc,progr=progr)
    })
  }else{
    progr <- progressor(along=ListIte)
    lapply(ListIte,Xfunc)  
  }
  print(Sys.time())
}

#-----------------------
nalist <- c(as.vector(varlist$V1),"TPES","POWER","Power_heat","Landuse","TFC_fuel","TFC_Sector","TFC_Ind","TFC_Tra","TFC_Res","TFC_Com")
allplot <- as.list(nalist)
plotflag <- as.list(nalist)
names(allplot) <- nalist
names(plotflag) <- nalist
allplotmerge <- as.list(nalist)
plotflagmerge <- as.list(nalist)
lst <- list()
lst$region <- region_load
lst$varlist <- as.list(as.vector(varlist$V1))
#lst$region <- "World"

#regional figure generation execution
exe_fig_make(lst$region,funcplotgen)
#cross-regional figure generation execution
exe_fig_make(lst$varlist,mergefigGen)



