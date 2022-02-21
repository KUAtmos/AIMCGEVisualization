# This program intends to visualize AIMHub results
#----------------------package installation and load ----------------------*
insflag <- 0
if(insflag==1){
  options(CRAN="http://cran.md.tsukuba.ac.jp/")
  library(devtools)
  devtools::install_github("tomwenseleers/export")
  install.packages('RDCOMClient', repos = 'http://www.omegahat.net/R/')
  liblist <- c("ggplot2","RColorBrewer","dplyr","sp","maptools","maps","ggradar","fmsb","tidyr","stringr","rJava","Rcpp","ReporteRsjars","ReporteRs","xlsx","officer","furrr","purrr","progressr")
  for(j in liblist){
    install.packages(j, dependencies = TRUE)
  }
  #install.packages("gdxrrw", dependencies = TRUE)
}

libloadlist <- c("gdxrrw","ggplot2","dplyr","reshape2","tidyr","maps","grid","RColorBrewer","R2PPT","RDCOMClient","cowplot","officer","export","purrr","furrr","progressr")
for(j in libloadlist){
  eval(parse(text=paste0("library(",j,")")))
}


#---------------switches to specify the run condition -----
filename <- "global_17" # filename should be "global_17","CHN","JPN"....
enduseflag <- 4   # If you would like to display AIM/Enduse outputs, make this parameter 1 otherwise 0.
enduseEneCost <- 0 # if you would like to display additional, energy system cost per GDP in the figure of GDP loss rate, make parameter 1 and otherwise 0.
dirCGEoutput <-"../../output/iiasa_database/gdx/"  # directory where the CGE output is located 
CGEgdxcopy <- 0 # if you would like to copy and store the CGE IAMC template file make this parameter 1, otherwise 0.
dirEnduseoutput <-"../../../Enduse/output/globalCGEInt/cons/main/"  # directory where the CGE output is located 
parallelmode <- 1 #Switch for parallel process. if you would like to use multi-processors assign 1 otherwise 0.
#parallelmode <- 0 #Switch for parallel process. if you would like to use multi-processors assign 1 otherwise 0.
threadsnum <- min(floor(availableCores()/2),24)
r2ppt <- 0 #Switch for ppt export. if you would like to export as ppt then assign 1 otherwise 0.
mergecolnum <- 4 #merge figure facet number of columns

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
outdir <- "../output/"
dirlist <- c(outdir,paste0(outdir,"ppt"),paste0(outdir,"merge"),paste0(outdir,"merge/png"),paste0(outdir,"merge/pngdet"))
for(dd in dirlist){
  if(file.exists(dd)){}else{dir.create(dd)}
}


file.copy(paste0(dirCGEoutput,"../../../AIMCGE/data/AIMHubData/main/IEAEBIAMCTemplate.gdx"), paste0("../data/IEAEBIAMCTemplate.gdx"),overwrite = TRUE)

linepalette <- c("#4DAF4A","#FF7F00","#377EB8","#E41A1C","#984EA3","#F781BF","#8DD3C7","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#7f878f","#A65628","#FFFF33","black")
landusepalette <- c("#8DD3C7","#FF7F00","#377EB8","#4DAF4A","#A65628")
#linepalette <- c("Baseline"="#4DAF4A","GlobalOptimalZero"="#FF7F00","NDC+Zero"="#377EB8","#E41A1C","#984EA3","#F781BF","#8DD3C7","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#7f878f","#A65628","#FFFF33")

#File loading and parameter configuration
fileloadlist <- read.table("../data/loadfilelist.txt", sep="\t",header=T, stringsAsFactors=F)
for(i in 1:nrow(fileloadlist)){
  eval(parse(text=paste0(fileloadlist[i,]$paraname," <- read.table('",fileloadlist[i,]$filename,"', sep='\t',header=T, stringsAsFactors=F)")))
}
region_load <- as.vector(read.table("../data/region.txt", sep="\t",header=F, stringsAsFactors=F)$V1)
region <- region_load
varlist <- left_join(varlist_load,varalllist,by="V1")

#---IAMC tempalte loading and data merge
if(CGEgdxcopy==1){
  file.copy(paste0(dirCGEoutput,filename,"_IAMC.gdx"), paste0("../modeloutput/",filename,"_IAMC.gdx"),overwrite = TRUE)
  CGEload0 <- rgdx.param(paste0('../modeloutput/',filename,"_IAMC.gdx"),'IAMC_Template') 
}else{  
  CGEload0 <- rgdx.param(paste0(dirCGEoutput,filename,"_IAMC.gdx"),'IAMC_Template') 
}
Getregion <- as.vector(unique(CGEload0$REMF))
if(length(Getregion)==1){region <- Getregion}
CGEload1 <- CGEload0 %>% rename("Value"=IAMC_Template,"Variable"=VEMF) %>% 
  left_join(scenariomap,by="SCENARIO") %>% filter(SCENARIO %in% as.vector(scenariomap[,1]) & REMF %in% region) %>% 
  select(-SCENARIO) %>% rename(Region="REMF",SCENARIO="Name",Y="YEMF")

#Enduse loading

if(enduseflag>=1){
  for(ii in 1:enduseflag){
    fileid <- ii
    if(file.exists(paste0(dirEnduseoutput,"../../../globalCGEInt",fileid,"/cons/main/merged_output.gdx"))){
      file.copy(paste0(dirEnduseoutput,"../../../globalCGEInt",fileid,"/cons/main/merged_output.gdx"), paste0("../modeloutput/AIMEnduseG",ii,".gdx"),overwrite = TRUE)
      eval(parse(text=paste0("EnduseGloadX0_",ii," <- rgdx.param(paste0('../modeloutput/AIMEnduseG",ii,".gdx'),'data_all')  %>% rename('SCENARIO'=Sc,'Region'=Sr,'Variable'=Sv,'Y'=Sy,'Value'=data_all)  %>% mutate(Model=paste0('Enduse[Global]-',",ii,"))%>% left_join(scenariomap2,by='SCENARIO')")))
      eval(parse(text=paste0("EnduseGloadX1_",ii," <- EnduseGloadX0_",ii,"  %>% filter(SCENARIO %in% as.vector(scenariomap2[,1]) & Region %in% region) %>% select(-SCENARIO) %>% rename(SCENARIO='Name') %>% select(Region,Variable,Y,Value,SCENARIO,Model)")))
      if(enduseEneCost==1){
        eval(parse(text=paste0("EnduseGload_cost <- filter(EnduseGloadX0_",ii,", Variable %in% c('Pol_Cos_Add_Tot_Ene_Sys_Cos','GDP_MER') & SCENARIO %in% as.vector(scenariomap2[,1]) & Region %in% region) %>% spread(key=Variable,value=Value,fill=0)")))
        EnduseGload_cost$Pol_Cos_GDP_los_rat <- EnduseGload_cost$Pol_Cos_Add_Tot_Ene_Sys_Cos/EnduseGload_cost$GDP_MER*100
        eval(parse(text=paste0("EnduseGloadX2_",ii," <- rbind(EnduseGloadX1_",ii,", select(EnduseGload_cost,-GDP_MER,-Pol_Cos_Add_Tot_Ene_Sys_Cos,-SCENARIO) %>% rename(Value=Pol_Cos_GDP_los_rat) %>% mutate(Variable='Pol_Cos_GDP_Los_rat')%>% rename(SCENARIO='Name'))")))
      }else{
        eval(parse(text=paste0("EnduseGloadX2_",ii," <- EnduseGloadX1_",ii)))
      }
      eval(parse(text=paste0("EnduseGloadX1_",ii," <- 0")))
    }
  }
}
allmodel0 <- CGEload1  
if(enduseflag>=1){
  for(ii in 1:enduseflag){
    eval(parse(text=paste0("allmodeltmp <- rbind(allmodel0,EnduseGloadX2_",ii,")")))
    allmodel0 <- allmodeltmp
  }
}
allmodel0$Y <- as.numeric(levels(allmodel0$Y))[allmodel0$Y]

#IEA energy balance information
IEAEB0 <- rgdx.param('../data/IEAEBIAMCTemplate.gdx','IAMCtemp17') %>% rename("Value"=IAMCtemp17,"Variable"=VEMF,"Y"=St,"Region"=Sr17,"SCENARIO"=SceEneMod) %>%
  select(Region,Variable,Y,Value,SCENARIO) %>% filter(Region %in% region) %>% mutate(Model="Reference")
IEAEB0$Y <- as.numeric(levels(IEAEB0$Y))[IEAEB0$Y]
IEAEB1 <- filter(IEAEB0,Y<=2015 & Y>=1990)

allmodel <- rbind(allmodel0,IEAEB1)  
#%>% filter(Model!="Enduse[Global]-4" & Model!="CGEEnd5" & Model!="CGEEnd4")
maxy <- max(allmodel$Y)
linepalettewName <- linepalette
names(linepalettewName) <- unique(allmodel$SCENARIO)

#---End of IAMC tempalte loading and data merge

#Decomposition analysis data load
Decom1 <- rgdx.param(paste0(dirCGEoutput,"../../global/",filename,"/gdx/analysis.gdx"),'Loss_dcp_gdp' ) %>% rename("value"=Loss_dcp_gdp,"Sector"=SCO2_S,"Element"=decele) 
flabel <- c("change in % of GDP","sectors")


#---functions
#function for regional figure generation
funcplotgen <- function(rr,progr){
  progr(message='region figures')

#---Line figures
  for (i in 1:nrow(varlist)){
    if(nrow(filter(allmodel,Variable==varlist[i,1] & Region==rr & Model!="Reference"))>0){
      miny <- min(filter(allmodel,Variable==varlist[i,1] & Region==rr)$Y) 
      plot.0 <- ggplot() + 
        geom_line(data=filter(allmodel,Variable==varlist[i,1] & Model!="Reference"& Region==rr),aes(x=Y, y = Value , color=SCENARIO,group=interaction(SCENARIO,Model)),stat="identity") +
        geom_point(data=filter(allmodel,Variable==varlist[i,1] & Model!="Reference"& Region==rr),aes(x=Y, y = Value , color=SCENARIO,shape=Model),size=3.0,fill="white") +
        MyThemeLine + scale_color_manual(values=linepalettewName) + scale_x_continuous(breaks=seq(miny,maxy,10)) +
        scale_shape_manual(values = 1:length(unique(allmodel$Model))) +
        xlab("year") + ylab(varlist[i,4])  +  ggtitle(paste(rr,varlist[i,3],sep=" ")) +
        annotate("segment",x=miny,xend=maxy,y=0,yend=0,linetype="dashed",color="grey")+ 
        theme(legend.title=element_blank()) 
      if(length(scenariomap$SCENARIO)<20){
        plot.0 <- plot.0 +
        geom_point(data=filter(allmodel,Variable==varlist[i,1] & Model=="Reference"& Region==rr),aes(x=Y, y = Value) , color="black",shape=0,size=2.0,fill="grey") 
      }
      if(varlist[i,2]==1){
        outname <- paste0(outdir,rr,"/png/",varlist[i,1],".png")
      }else{
        outname <- paste0(outdir,rr,"/pngdet/",varlist[i,1],".png")
      }
      ggsave(plot.0, file=outname, dpi = 150, width=10, height=6,limitsize=FALSE)
      allplot[[nalist[i]]] <- plot.0
    }
    plotflag[[nalist[i]]] <- nrow(filter(allmodel,Variable==varlist[i,1] & Model!="Reference"& Region==rr))
  }
  #---merged figures
  #Final energy consumption area
  pp_tfc <- plot_grid(allplot[["TFC_Ind"]],allplot[["TFC_Tra"]],allplot[["TFC_Res"]],allplot[["TFC_Com"]],ncol=2,align = "hv")
  ggsave(pp_tfc, file=paste0(outdir,rr,"/merge/tfc.png"), width=9*2, height=(floor(length(unique(allmodel$SCENARIO))/4+1)*3+2)*3,limitsize=FALSE)
  p_legend1 <- gtable::gtable_filter(ggplotGrob(allplot[["Fin_Ene"]]), pattern = "guide-box")
  #Final energy consumption by sectors and fuels
  pp_tfcind <- plot_grid(allplot[["Fin_Ene"]] + theme(legend.position="none"),allplot[["Fin_Ene_Ind"]] + theme(legend.position="none"),allplot[["Fin_Ene_Tra"]] + theme(legend.position="none"),allplot[["Fin_Ene_Res"]] + theme(legend.position="none"),allplot[["Fin_Ene_Com"]] + theme(legend.position="none"),
                         allplot[["Fin_Ene_Ele_Heat"]] + theme(legend.position="none"),allplot[["Fin_Ene_Gas"]] + theme(legend.position="none"),allplot[["Fin_Ene_Liq"]] + theme(legend.position="none"),allplot[["Fin_Ene_SolidsCoa"]] + theme(legend.position="none"),allplot[["Fin_Ene_SolidsBio"]] + theme(legend.position="none"),
                         allplot[["Fin_Ene_Ind_Ele_Heat"]] + theme(legend.position="none"),allplot[["Fin_Ene_Ind_Gas"]] + theme(legend.position="none"),allplot[["Fin_Ene_Ind_Liq"]] + theme(legend.position="none"),allplot[["Fin_Ene_Ind_SolidsCoa"]] + theme(legend.position="none"),allplot[["Fin_Ene_Ind_SolidsBio"]] + theme(legend.position="none"),
                         allplot[["Fin_Ene_Com_Ele_Heat"]] + theme(legend.position="none"),allplot[["Fin_Ene_Com_Gas"]] + theme(legend.position="none"),allplot[["Fin_Ene_Com_Liq"]] + theme(legend.position="none"),allplot[["Fin_Ene_Com_SolidsCoa"]] + theme(legend.position="none"),allplot[["Fin_Ene_Com_SolidsBio"]] + theme(legend.position="none"),
                         allplot[["Fin_Ene_Res_Ele_Heat"]] + theme(legend.position="none"),allplot[["Fin_Ene_Res_Gas"]] + theme(legend.position="none"),allplot[["Fin_Ene_Res_Liq"]] + theme(legend.position="none"),allplot[["Fin_Ene_Res_SolidsCoa"]] + theme(legend.position="none"),allplot[["Fin_Ene_Res_SolidsBio"]] + theme(legend.position="none"),
                         allplot[["Fin_Ene_Tra_Ele"]] + theme(legend.position="none"),allplot[["Fin_Ene_Tra_Liq_Bio"]] + theme(legend.position="none"),allplot[["Fin_Ene_Tra_Liq_Oil"]] + theme(legend.position="none"),NULL,p_legend1,
                         nrow=6,rel_widths =c(1,1,1,1,1),align = "hv")
  ggsave(pp_tfcind, file=paste0(outdir,rr,"/merge/tfcind.png"), width=25, height=20,limitsize=FALSE)
  #Main indicators
  pp_main <- plot_grid(allplot[["GDP_MER"]] + theme(legend.position="none"),allplot[["POP"]] + theme(legend.position="none"),allplot[["Tem_Glo_Mea"]],
                       allplot[["Emi_CO2_Ene_and_Ind_Pro"]] + theme(legend.position="none"),allplot[["Emi_CO2"]] + theme(legend.position="none"),allplot[["Emi_Kyo_Gas"]],
                       allplot[["Pol_Cos_GDP_Los_rat"]] + theme(legend.position="none"),allplot[["Pol_Cos_Cns_Los_rat"]] + theme(legend.position="none"),allplot[["Prc_Car"]],
                       allplot[["Pop_Ris_of_Hun"]] + theme(legend.position="none"),allplot[["Prc_Prm_Ene_Oil"]] + theme(legend.position="none"),allplot[["Prc_Sec_Ene_Ele"]],
                       nrow=4,rel_widths =c(1,1,1.5),align = "hv")
  ggsave(pp_main, file=paste0(outdir,rr,"/merge/main.png"), width=15, height=15,limitsize=FALSE)
  
  #----r2ppt
  #The figure should be prearranged before going this ppt process since emf file type does not accept size changes. 
  #If you really needs ppt slide, you first output png and then paste it.
  pptlist <- c("Fin_Ene","Fin_Ene_Ele_Heat","Fin_Ene_Gas","Fin_Ene_Liq","Fin_Ene_Solids","Fin_Ene_Res","Fin_Ene_Com","Fin_Ene_Tra","Fin_Ene_Ind","Emi_CO2_Ene_and_Ind_Pro","Pol_Cos_GDP_Los_rat","Prc_Car")
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
#Function for Decomposition
funcDecGen <- function(rr,progr){
  progr(message='region figures')
  Decom2 <- Decom1 %>% filter(SCENARIO %in% scenariomap$SCENARIO & Y %in% c(2030,2050) & Sector %in% c("IND","SER","PWR","OEN","TRS","AGR") & R==rr) 
  plotdec <- ggplot() + geom_bar(data=filter(Decom2,Element %in% c("fd_output","output_va","va","residual1")),aes(x=Sector, y = value*100 , fill=Element), stat="identity") +
    geom_point(data=filter(Decom2,Element %in% c("fd")),aes(x=Sector, y = value*100 ),color="black", stat="identity") +
    ylab(flabel[1]) + xlab(flabel[2]) +labs(fill="") +
    guides(fill=guide_legend(reverse=TRUE)) + 
    MyThemeLine + theme(legend.position="bottom", text=element_text(size=12))+
    guides(fill=guide_legend(ncol=5))+ggtitle(paste0(rr," decomposition"))+
    facet_grid(Y~SCENARIO,scales="free_x") + annotate("segment",x=0,xend=6,y=0,yend=0,linetype="dashed",color="grey")
  outname <- paste0(outdir,rr,"/merge/","decomp.png")
  ggsave(plotdec, file=outname, width=floor(length(unique(Decom2$SCENARIO))/2+1)*4, height=10,limitsize=FALSE)    
}
#function for regional area figure generation
funcAreaPlotGen <- function(rr,progr){
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
      facet_wrap(Model ~ SCENARIO,ncol=mergecolnum) + scale_fill_manual(values=colorpal) + 
      annotate("segment",x=miny,xend=maxy,y=0,yend=0,linetype="solid",color="grey") + theme(legend.position='bottom')+
      ggtitle(paste0(rr,areamappara[j,]$Class))
    if(nrow(XX2)>=1){
      plot3 <- plot2 +    geom_area(data=XX2,aes(x=Y, y = Value , fill=reorder(Ind,-order)), stat="identity")
    }else{
      plot3 <- plot2
    }
    allplot[[areamappara$Class[j]]] <- plot3 
    outname <- paste0(outdir,rr,"/merge/",areamappara[j,1],".png")
    ggsave(plot3, file=outname, dpi = 450, width=9, height=floor(length(unique(XX$SCENARIO))/4+1)*8+2,limitsize=FALSE)
    plotflag[[areamappara$Class[j]]] <- nrow(XX)  
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
      scale_shape_manual(values = 1:length(unique(allmodel$Model))) +
      xlab("year") + ylab(varlist$V3[varlist$V1==ii])  +  ggtitle(paste("Multi-regions",varlist$V2.y[varlist$V1==ii],sep=" ")) +
      annotate("segment",x=miny,xend=maxy,y=0,yend=0,linetype="dashed",color="grey")+ 
      theme(legend.title=element_blank()) +facet_wrap(~Region,scales="free")
    if(length(scenariomap$SCENARIO)<20){
      plot.0 <- plot.0 +
        geom_point(data=filter(allmodel,Variable==ii & Model=="Reference"),aes(x=Y, y = Value) , color="black",shape=0,size=2.0,fill="grey") 
    }
    if(length(varlist$V2.x[varlist$V1==ii])==1){
      outname <- paste0(outdir,"merge","/png/",ii,".png")
    }else{
      outname <- paste0(outdir,"merge","/pngdet/",ii,".png")
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
#Parameter configuration for iterations
nalist <- c(as.vector(varlist$V1),"TPES","POWER","Power_heat","Landuse","TFC_fuel","TFC_Sector","TFC_Ind","TFC_Tra","TFC_Res","TFC_Com","Investment")
allplot <- as.list(nalist)
plotflag <- as.list(nalist)
names(allplot) <- nalist
names(plotflag) <- nalist
allplotmerge <- as.list(nalist)
plotflagmerge <- as.list(nalist)
lst <- list()
lst$region <- region_load
#lst$region <- "World"
lst$varlist <- as.list(as.vector(varlist$V1))

#Creat directories
for(rr in lst$region){
  regoutdir <- paste0("../output/",rr)
  dirlist <- c(regoutdir,paste0(regoutdir,"/png"),paste0(regoutdir,"/pngdet"),paste0(regoutdir,"/ppt"),paste0(regoutdir,"/merge"))
  for(dd in dirlist){
    if(file.exists(dd)){}else{dir.create(dd)}
  }
}

#regional figure generation execution
exe_fig_make(lst$region,funcplotgen)
#regional Decomposition figure generation execution
exe_fig_make(lst$region,funcDecGen)
#regional area figure generation execution
exe_fig_make(lst$region,funcAreaPlotGen)


#cross-regional figure generation execution
exe_fig_make(lst$varlist,mergefigGen)



