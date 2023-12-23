# This program intends to visualize AIMHub results
#----------------------package installation and load ----------------------*
insflag <- 0
if(insflag==1){
  options(CRAN="http://cran.md.tsukuba.ac.jp/")
  install.packages("devtools", dependencies = TRUE)
  library(devtools)
  devtools::install_github("tomwenseleers/export")
  install.packages('RDCOMClient', repos = 'http://www.omegahat.net/R/')
  liblist <- c("reshape2","cowplot","ggplot2","RColorBrewer","dplyr","sp","maptools","maps","ggradar","fmsb","tidyr","stringr","rJava","Rcpp","ReporteRsjars","ReporteRs","xlsx","officer","furrr","purrr","progressr")
  for(j in liblist){
    install.packages(j, dependencies = TRUE)
  }
  #install.packages("gdxrrw", dependencies = TRUE)
}

libloadlist <- c("gdxrrw","ggplot2","dplyr","reshape2","tidyr","maps","grid","RColorBrewer","cowplot","hms","purrr","furrr","progressr")
for(j in libloadlist){
  eval(parse(text=paste0("library(",j,")")))
}
args <- commandArgs(trailingOnly = TRUE)
#arguments are: 1:gams sys,2:number of CPU, 3:visualizaton scenario name specification auto or not, 4:file location, 5:submodule switch, 6:enduse iteration switch, 7: GDX file name, 8: region code/global
default_args <- c("/opt/gams/gams37.1_linux_x64_64_sfx", min(floor(availableCores()/2),24), "on", "global/global_17","1","0","global_17_IAMC","global")   # Default value but gams path should be modified if GUI based R is used
#default_args <- c("/opt/gams/gams37.1_linux_x64_64_sfx", min(floor(availableCores()/2),24), "on", "country/CHN","2","1","IAMCTemplate_Iteon_CHN","CHN")   # Default value but gams path should be modified if GUI based R is used

default_flg <- is.na(args[1:8])
args[default_flg] <- default_args[default_flg]
gams_sys_dir <- as.character(args[1])
AscenarionameAuto <- as.character(args[3])
igdx(gams_sys_dir)
Iterationflag <- as.numeric(args[6])   # If you would like to display AIM/Enduse outputs, make this parameter 1 otherwise 0.
decompositionflag <- 0  #if you would like to run decomposition analysis turn on 1, otherwise 0.
threadsnum <-  as.numeric(args[2])
sizememory <- 1000*1024^2 
options(future.globals.maxSize= sizememory)


#---------------switches to specify the run condition -----
submodule <- as.numeric(args[5]) #0 if this repository cloned to the AIMHub dir
if(submodule==1){
	maindirloc <- "../../"
	outdir <- paste0("../../../../output/fig_",args[8],"/") #Output directory 
	AIMHubdir <- "../../../" 
	VarListPath <- paste0(AIMHubdir,"tools/iiasa_data_submission/data/all_list.txt")
}else if(submodule==0){
	maindirloc <- ""
	outdir <- "../output/" 
	AIMHubdir <- "../../" 
	VarListPath <- paste0(AIMHubdir,"tools/iiasa_data_submission/data/all_list.txt")
}else if(submodule==2){
  maindirloc <- "../../"
  if(Iterationflag==0){
    Itename <- "Iteoff"
    VisualizationScenarioFile <- "VisualizationScenariomap"
  }else{
    Itename <- "Iteon"
    VisualizationScenarioFile <- "VisualizationScenariomapIte"
  }
  outdir <- paste0("../../../../../../IntTool/output/fig_",args[8],"/") #Output directory 
  AIMHubdir <- "../../../" 
	VarListPath <- paste0(outdir,"../../define/iamctemp/VariableFullList.txt")
}
outdirmd <- paste0(outdir,"modeloutput/") #output direcotry to save temporary GDX file
filename <- args[7] # filename should be "global_17","CHN","JPN"....
CGEgdxcopy <- 0 # if you would like to copy and store the CGE IAMC template file make this parameter 1, otherwise 0.
parallelmode <- 1 #Switch for parallel process. if you would like to use multi-processors assign 1 otherwise 0.
print(threadsnum) 
r2ppt <- 0 #Switch for ppt export. if you would like to export as ppt then assign 1 otherwise 0.
mergecolnum <- 6 #merge figure facet number of columns
RegSpec <- 0 #Regional specification if turned into 1, the regional list is loaded for the regional plot
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
    axis.line=element_line(colour="black"),
    panel.background=element_rect(fill = "white"),
    panel.grid.major=element_blank(),
    strip.background=element_rect(fill="white", colour="white"),
    strip.text.x = element_text(size=10, colour = "black", angle = 0,face="bold"),
    axis.text.x=element_text(angle=45, vjust=0.9, hjust=1, margin = unit(c(t = 0.3, r = 0, b = 0, l = 0), "cm")),
    axis.text.y=element_text(margin = unit(c(t = 0, r = 0.3, b = 0, l = 0), "cm")),
#    axis.text.x=element_text(size = 10,angle=45, vjust=0.9, hjust=1, margin = unit(c(t = 0.3, r = 0, b = 0, l = 0), "cm")),
#    axis.text.y=element_text(size = 10,margin = unit(c(t = 0, r = 0.3, b = 0, l = 0), "cm")),
#    legend.text = element_text(size = 10),
#    legend.title = element_text(size = 10),
    axis.ticks.length=unit(-0.15,"cm")
  )

#-- data load
dirlist <- c(outdir,paste0(outdir,"data"),paste0(outdir,"byRegion"),paste0(outdirmd),paste0(outdir,"multiRegR5"),paste0(outdir,"ppt"),paste0(outdir,"misc"),paste0(outdir,"data"))
for(dd in dirlist){
  if(file.exists(dd)){}else{dir.create(dd)}
}


file.copy(paste0(AIMHubdir,"data/AIMHubData/main/IEAEBIAMCTemplate.gdx"), paste0("../data/IEAEBIAMCTemplate.gdx"),overwrite = TRUE)

linepalette <- c("#4DAF4A","#FF7F00","#377EB8","#E41A1C","#984EA3","#F781BF","#8DD3C7","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#7f878f","#A65628","#FFFF33","black")
landusepalette <- c("#8DD3C7","#FF7F00","#377EB8","#4DAF4A","#A65628")
#linepalette <- c("Baseline"="#4DAF4A","GlobalOptimalZero"="#FF7F00","NDC+Zero"="#377EB8","#E41A1C","#984EA3","#F781BF","#8DD3C7","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FFED6F","#7f878f","#A65628","#FFFF33")

#File loading and parameter configuration
varalllist <- read.table(VarListPath, sep="\t",header=F, stringsAsFactors=F)
varlist_load <- read.table('../data/varlist.txt',sep='\t',header=T)
areamap <- read.table('../data/Areafigureorder.txt',sep='\t',header=T)
areamappara <- read.table('../data/area.map',sep='\t',header=T)
R5R_load <- read.table('../data/regionR5.txt',sep='\t',header=F)
R17R_load <- read.table('../data/regionR17.txt',sep='\t',header=F)
region_load <- as.vector(read.table("../data/region.txt", sep="\t",header=F, stringsAsFactors=F)$V1)
region <- region_load
R5R <- as.vector(R5R_load$V1)
R17R <- as.vector(R17R_load$V1)
R17p5R <- c(R5R,R17R[-18])
varlist <- left_join(varlist_load,varalllist,by=c("V1"))
Ylist <- seq(2010,2100,by=5)

areapaletteload <- select(areamap,Class,Ind,color) %>% rename(V0=Class,V1=Ind,V2=color) 

#---IAMC tempalte loading and data merge
dirCGEoutput <- paste0(maindirloc,"../../output/iiasa_database/gdx/")  # directory where the CGE output is located 
if(submodule!=2){
  if(AscenarionameAuto=="on"){  #scenario mapping specification
    scenariomap_load <- read.table(paste0(dirCGEoutput,'../../',args[4],'/txt/scenario_list.txt'), header=F,stringsAsFactors=F)
    scenariomap <- cbind(scenariomap_load,scenariomap_load,"CGE")
    names(scenariomap) <- c("SCENARIO","Name","ModName")
  }else{
    scenariomap <- read.table('../data/scenariomap.map',sep='\t',header=T)
  }
  
  if(CGEgdxcopy==1){ # Data file loading
    file.copy(paste0(dirCGEoutput,filename,".gdx"), paste0(outdirmd,filename,".gdx"),overwrite = TRUE)
    CGEload0 <- rgdx.param(paste0(outdirmd ,filename,".gdx"),'IAMC_Template') 
  }else{  
    CGEload0 <- rgdx.param(paste0(dirCGEoutput,filename,".gdx"),'IAMC_Template') 
  }
  CGEload1 <- CGEload0 %>% rename("Value"=IAMC_Template,"Var"=VEMF) %>% 
    left_join(scenariomap,by="SCENARIO") %>% filter(SCENARIO %in% as.vector(scenariomap[,1])) %>% 
    select(-SCENARIO) %>% rename(Region="REMF",SCENARIO="Name",Y="YEMF")
}else{
  if(AscenarionameAuto=="on"){  #scenario mapping specification
    scenariomap_load <- rgdx.set(paste0(outdir,'/../iamc/',filename,'.gdx'),'SCENARIO') 
    scenariomap <- cbind(scenariomap_load,scenariomap_load)
    names(scenariomap) <- c("SCENARIO","Name")
  }else{
    scenariomap <- read.table(paste0(outdir,'../../define/iamctemp/VisualizationScenariomap.map'),sep='\t',header=T)
  }
  CGEload0 <- rgdx.param(paste0(outdir,'/../iamc/',filename,'.gdx'),'MergedIAMC') 
  CGEload1 <- CGEload0 %>% rename("Value"=mergedIAMC,"Var"=VIAMC,Region="RIAMC",ModName="Modelset") %>% left_join(scenariomap,by="SCENARIO") %>% select(Region,Var,Y,Value,SCENARIO,ModName)
}
Getregion <- as.vector(unique(CGEload1$Region))
if(length(Getregion)==1){region <- Getregion}else{region <- R17p5R}
CGEload1$Y <- as.numeric(levels(CGEload1$Y))[CGEload1$Y]
CGEload0 <- 0
allmodel0 <- CGEload1


#IEA energy balance information
IEAEB0 <- rgdx.param('../data/IEAEBIAMCTemplate.gdx','IAMCtemp17') %>% rename("Value"=IAMCtemp17,"Var"=VEMF,"Y"=St,"Region"=Sr17,"SCENARIO"=SceEneMod) %>%
  select(Region,Var,Y,Value,SCENARIO) %>% filter(Region %in% c(R5R,R17R)) %>% mutate(ModName="Reference")
IEAEB0$Y <- as.numeric(levels(IEAEB0$Y))[IEAEB0$Y]
IEAEB1 <- filter(IEAEB0,Y<=2020 & Y>=1990)

allmodel <- rbind(allmodel0,IEAEB1) %>% select(ModName,Region,Var,SCENARIO,Y,Value) 
maxy <- max(allmodel$Y)
#maxy <- 2050
linepalettewName <- linepalette[1:length(unique(allmodel$SCENARIO))]
names(linepalettewName) <- unique(allmodel$SCENARIO)
allmodel <- filter(allmodel,Y <= maxy)

#Extract data
ExtData <- filter(allmodel0,Var %in% varlist$V1) %>% left_join(unique(varlist %>% rename(Var=V1,Variable=V2.y,Unit=V3) %>% select(Var,Variable,Unit))) %>% 
  select(-Var) %>% rename(Model=ModName,Year=Y) %>%
  select(Model,SCENARIO,Region,Variable,Unit,Year,Value) %>% 
  spread(value=Value,key=Year)
write.csv(x = ExtData, row.names = FALSE,file = paste0(outdir,"/data/exportdata.csv"))
ExtData <- 0
#---End of IAMC tempalte loading and data merge


#---functions
#function for regional figure generation
funcplotgen <- function(rr,progr){
  progr(message='region figures')
  Data4plot0 <- filter(allmodel,Region==rr)
  #---Line figures
#  for (rr in region_spe){ Data4plot0 <- filter(allmodel,Region==rr) #For debug
  for (i in 1:nrow(varlist)){
    Data4plot <- filter(Data4plot0,Var==varlist$V1[i])
    if(nrow(filter(Data4plot,ModName!="Reference"))>0){
      miny <- min(Data4plot$Y) 
      linepalettewName1 <- linepalette[1:length(unique(Data4plot$SCENARIO))]
      names(linepalettewName1) <- unique(Data4plot$SCENARIO)
      numitem <- length(as.vector(unique(Data4plot$SCENARIO)))+length(as.vector(unique(Data4plot$ModName))) #Get number of items
      plot.0 <- ggplot() + 
        geom_line(data=filter(Data4plot,ModName!="Reference"),aes(x=Y, y = Value , color=SCENARIO,group=interaction(SCENARIO,ModName)),stat="identity") +
        geom_point(data=filter(Data4plot,ModName!="Reference"),aes(x=Y, y = Value , color=SCENARIO,shape=ModName),size=1.5,fill="white") +
        MyThemeLine + scale_color_manual(values=linepalettewName1) + scale_x_continuous(breaks=seq(miny,maxy,10)) +
        scale_shape_manual(values = 1:length(unique(allmodel$ModName))) +
        xlab("year") + ylab(paste0(varlist$V2.y[i],"(",varlist$V3[i],")") ) +  ggtitle(paste0(rr,expression("\n"),varlist$V2.y[i])) +
        annotate("segment",x=miny,xend=maxy,y=0,yend=0,linetype="dashed",color="grey")+ 
        theme(legend.title=element_blank()) 
      if(length(scenariomap$SCENARIO)<40){
        plot.0 <- plot.0 +
        geom_point(data=filter(Data4plot,ModName=="Reference"),aes(x=Y, y = Value) , color="black",shape=0,size=1.5,fill="grey") 
      }
      outname <- paste0(outdir,"byRegion/",rr,"/png/",varlist$V1[i],"_",rr,".png")
      ggsave(plot.0, file=outname, dpi = 150, width=max(1,numitem/10)*5, height=max(1,numitem/10)*3.5,limitsize=FALSE)
      allplot[[nalist[i]]] <- plot.0
      allplot_nonleg[[nalist[i]]] <- plot.0+ theme(legend.position="none")
    }
    plotflag[[nalist[i]]] <- nrow(filter(Data4plot,ModName!="Reference"))
  }
  #---merged figures
  #Final energy consumption by sectors and fuels
  p_legend1 <- gtable::gtable_filter(ggplotGrob(allplot[["Fin_Ene"]]), pattern = "guide-box")
  pp_tfcind <- plot_grid(
                         allplot_nonleg[["Fin_Ene"]],    allplot_nonleg[["Fin_Ene_Ele_Heat"]],    allplot_nonleg[["Fin_Ene_Liq_and_Gas"]],
                           allplot_nonleg[["Fin_Ene_Gas"]],allplot_nonleg[["Fin_Ene_Liq"]],allplot_nonleg[["Fin_Ene_SolidsCoa"]],allplot_nonleg[["Fin_Ene_SolidsBio"]],allplot_nonleg[["Fin_Ene_Hyd"]],
                         allplot_nonleg[["Fin_Ene_Ind"]],allplot_nonleg[["Fin_Ene_Ind_Ele_Heat"]],allplot_nonleg[["Fin_Ene_Ind_Liq_and_Gas"]],
                           allplot_nonleg[["Fin_Ene_Ind_Gas"]],allplot_nonleg[["Fin_Ene_Ind_Liq"]],allplot_nonleg[["Fin_Ene_Ind_SolidsCoa"]],allplot_nonleg[["Fin_Ene_Ind_SolidsBio"]],allplot_nonleg[["Fin_Ene_Ind_Hyd"]],
                         allplot_nonleg[["Fin_Ene_Com"]],allplot_nonleg[["Fin_Ene_Com_Ele_Heat"]],allplot_nonleg[["Fin_Ene_Com_Liq_and_Gas"]],
                           allplot_nonleg[["Fin_Ene_Com_Gas"]],allplot_nonleg[["Fin_Ene_Com_Liq"]],allplot_nonleg[["Fin_Ene_Com_SolidsCoa"]],allplot_nonleg[["Fin_Ene_Com_SolidsBio"]],allplot_nonleg[["Fin_Ene_Com_Hyd"]],
                         allplot_nonleg[["Fin_Ene_Res"]],allplot_nonleg[["Fin_Ene_Res_Ele_Heat"]],allplot_nonleg[["Fin_Ene_Res_Liq_and_Gas"]],
                           allplot_nonleg[["Fin_Ene_Res_Gas"]],allplot_nonleg[["Fin_Ene_Res_Liq"]],allplot_nonleg[["Fin_Ene_Res_SolidsCoa"]],allplot_nonleg[["Fin_Ene_Res_SolidsBio"]],allplot_nonleg[["Fin_Ene_Res_Hyd"]],
                         allplot_nonleg[["Fin_Ene_Tra"]],allplot_nonleg[["Fin_Ene_Tra_Ele"]],     allplot_nonleg[["Fin_Ene_Tra_Liq_and_Gas"]],
                           allplot_nonleg[["Fin_Ene_Tra_Gas"]],allplot_nonleg[["Fin_Ene_Tra_Liq_Bio"]],allplot_nonleg[["Fin_Ene_Tra_Liq_Oil"]],allplot_nonleg[["Fin_Ene_Tra_Hyd"]],p_legend1,
                         nrow=5,ncol=8,rel_widths =c(1,1,1,1,1,1,1,1),align = "hv")
  ggsave(pp_tfcind, file=paste0(outdir,"byRegion/",rr,"/merge/tfcind_",rr,".png"), width=30, height=20,limitsize=FALSE)
  #Main indicators
  p_legend1 <- gtable::gtable_filter(ggplotGrob(allplot[["GDP_MER"]]), pattern = "guide-box")
  if(nrow(filter(Data4plot0,Var=="Pol_Cos_GDP_Los_rat"))>0){
    p_legend2 <- gtable::gtable_filter(ggplotGrob(allplot[["Pol_Cos_GDP_Los_rat"]]), pattern = "guide-box")
    pp_main <- plot_grid(allplot_nonleg[["GDP_MER"]],allplot_nonleg[["Pop"]],allplot_nonleg[["Tem_Glo_Mea"]],p_legend1,
                       allplot_nonleg[["Emi_CO2_Ene_and_Ind_Pro"]],allplot_nonleg[["Emi_CO2"]],allplot_nonleg[["Emi_Kyo_Gas"]],p_legend1,
                       allplot_nonleg[["Pol_Cos_GDP_Los_rat"]],allplot_nonleg[["Pol_Cos_Cns_Los_rat"]],allplot_nonleg[["Prc_Car"]],p_legend2,
                       allplot_nonleg[["Pop_Ris_of_Hun"]],allplot_nonleg[["Prc_Prm_Ene_Oil"]],allplot_nonleg[["Prc_Sec_Ene_Ele"]],p_legend1,
                       nrow=4,rel_widths =c(1,1,1,0.7),align = "hv")
  }else{
    pp_main <- plot_grid(allplot_nonleg[["GDP_MER"]],allplot_nonleg[["Pop"]],allplot_nonleg[["Tem_Glo_Mea"]],p_legend1,
                         allplot_nonleg[["Emi_CO2_Ene_and_Ind_Pro"]],allplot_nonleg[["Emi_CO2"]],allplot_nonleg[["Emi_Kyo_Gas"]],p_legend1,
                         allplot_nonleg[["Pop_Ris_of_Hun"]],allplot_nonleg[["Prc_Prm_Ene_Oil"]],allplot_nonleg[["Prc_Sec_Ene_Ele"]],p_legend1,
                         nrow=3,rel_widths =c(1,1,1,0.7),align = "hv")
  }
  ggsave(pp_main, file=paste0(outdir,"byRegion/",rr,"/merge/main_",rr,".png"), width=15, height=15,limitsize=FALSE)

#Emissions
  p_legend1 <- gtable::gtable_filter(ggplotGrob(allplot[["Emi_CO2"]]), pattern = "guide-box")
  pp_main <- plot_grid(allplot_nonleg[["Emi_CO2"]],allplot_nonleg[["Emi_CH4"]],allplot_nonleg[["Emi_N2O"]],allplot_nonleg[["Emi_F_G"]],
                       allplot_nonleg[["Emi_Sul"]],allplot_nonleg[["Emi_NOx"]],allplot_nonleg[["Emi_BC"]],allplot_nonleg[["Emi_OC"]],
                       allplot_nonleg[["Emi_VOC"]],allplot_nonleg[["Emi_NH3"]],allplot_nonleg[["Emi_CO"]],allplot_nonleg[["Emi_Kyo_Gas"]],
                       allplot_nonleg[["Tem_Glo_Mea"]],allplot_nonleg[["Frc"]],p_legend1,
                       nrow=4,rel_widths =c(1,1,1,1),align = "hv")
  ggsave(pp_main, file=paste0(outdir,"byRegion/",rr,"/merge/Emissions_",rr,".png"), width=15, height=15,limitsize=FALSE)

#Land use
  p_legend1 <- gtable::gtable_filter(ggplotGrob(allplot[["Lan_Cov_Pst"]]), pattern = "guide-box")
  pp_main <- plot_grid(allplot_nonleg[["Lan_Cov_Cro"]],allplot_nonleg[["Lan_Cov_Pst"]],allplot_nonleg[["Lan_Cov_Cro_Ene_Cro"]],allplot_nonleg[["Lan_Cov_Frs_Aff_and_Ref"]],
                       allplot_nonleg[["Lan_Cov_Cro_Non_Ene_Cro"]],allplot_nonleg[["Lan_Cov_Frs"]],allplot_nonleg[["Lan_Cov_Frs_Nat_Frs"]],
                       allplot_nonleg[["Lan_Cov_Oth_Nat_Lan"]],allplot_nonleg[["Lan_Cov_Oth_Lan"]],allplot_nonleg[["Lan_Cov_Bui_Are"]],p_legend1,
                       nrow=3,rel_widths =c(1,1,1,1),align = "hv")
  ggsave(pp_main, file=paste0(outdir,"byRegion/",rr,"/merge/LanduseLine_",rr,".png"), width=15, height=15,limitsize=FALSE)

#Price
  p_legend1 <- gtable::gtable_filter(ggplotGrob(allplot[["Prc_Prm_Ene_Oil"]]), pattern = "guide-box")
  pp_main <- plot_grid(allplot_nonleg[["Prc_Prm_Ene_Oil"]],allplot_nonleg[["Prc_Prm_Ene_Gas"]],allplot_nonleg[["Prc_Prm_Ene_Bio"]],
                       allplot_nonleg[["Prc_Sec_Ene_Liq_Bio"]],allplot_nonleg[["Prc_Sec_Ene_Ele_MWh"]],allplot_nonleg[["Prc_Sec_Ene_Hyd"]],
                       allplot_nonleg[["Prc_Fin_Ene_Res_and_Com_Res"]],allplot_nonleg[["Prc_Fin_Ene_Res_and_Com_Res_Ele"]],allplot_nonleg[["Prc_Fin_Ene_Res_and_Com_Res_wo_car_pri"]],
                       allplot_nonleg[["Prc_Agr_NonEneCro_Ind"]],allplot_nonleg[["Prc_Agr_Liv_Ind"]],p_legend1,
                       nrow=4,rel_widths =c(1,1,1),align = "hv")
  ggsave(pp_main, file=paste0(outdir,"byRegion/",rr,"/merge/Price_",rr,".png"), width=15, height=15,limitsize=FALSE)

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
#function for regional area figure generation
funcAreaPlotSpe <- function(ZZ,ZZ2,ZZ3,AreaItem){
  miny <- min(ZZ$Y,ZZ2$Y) 
  na.omit(ZZ$Value)
  ylab1 <- paste0(areamappara$Var[areamappara$Class==AreaItem], " (", areamappara$Unit[areamappara$Class==AreaItem], ")")
  xlab1 <- areamappara$Var[areamappara$Class==AreaItem]

  areapaletteArea <- filter(areapaletteload,V0==AreaItem & V1 %in% unique(ZZ$Ind))$V2
  names(areapaletteArea) <- filter(areapaletteload,V0==AreaItem & V1 %in% unique(ZZ$Ind))$V1
  colorpal <- areapaletteArea 
  
  plotX <- ggplot() + 
    geom_area(data=filter(ZZ,Y<=maxy),aes(x=Y, y = Value , fill=reorder(Ind,-order)), stat="identity") + 
    ylab(ylab1) + xlab(xlab1) +labs(fill="")+ guides(fill=guide_legend(reverse=TRUE)) + MyThemeLine +
    theme(legend.position="bottom", text=element_text(size=12),  
          axis.text.x=element_text(angle=45, vjust=0.9, hjust=1, size = 12)) +
    guides(fill=guide_legend(ncol=5)) + scale_x_continuous(breaks=seq(miny,maxy,10)) + scale_fill_manual(values=colorpal) +
    annotate("segment",x=miny,xend=maxy,y=0,yend=0,linetype="solid",color="grey") + theme(legend.position='bottom')+
    geom_line(data=filter(ZZ3,Y<=maxy),aes(x=Y, y = Value ), color="black",linetype="dashed",size=1.2)
  if(nrow(ZZ2)>=1){
    plotY <- plotX +    geom_area(data=ZZ2,aes(x=Y, y = Value , fill=reorder(Ind,-order)), stat="identity", alpha=0.7)
  }else{
    plotY <- plotX
  }
  
  return(plotY)
}

funcAreaPlotGen <- function(rr,progr){
  Data4plot <- filter(allmodel,Region==rr)
#  for( rr in as.vector(region_load)){
  for(AreaItem in as.vector(areamappara$Class)){
    XX <- Data4plot %>% filter(Var %in% as.vector(areamap$Var)) %>% left_join(areamap,by="Var") %>% ungroup() %>% 
      filter(Class==AreaItem & ModName!="Reference") %>% select(ModName,SCENARIO,Ind,Y,Value,order)  %>% arrange(order)
    if(nrow(XX)>0){
      XX2 <- Data4plot %>% filter(Var %in% as.vector(areamap$Var)) %>% left_join(areamap,by="Var") %>% ungroup() %>% 
        filter(Class==AreaItem & ModName=="Reference") %>% select(-SCENARIO,-ModName,Ind,Y,Value,order)  %>% arrange(order)%>%
        filter(Y<=2015)
      XX3 <- Data4plot %>% filter(Var==areamappara$lineVar[areamappara$Class==AreaItem] & ModName!="Reference") %>% select(ModName,SCENARIO,Var,Y,Value)
      numitem <- length(as.vector(unique(Data4plot$ModName))) #Get number of items
      plot1 <- funcAreaPlotSpe(XX,XX2,XX3,AreaItem)
      plot3 <- plot1 + ggtitle(paste(rr,AreaItem,sep=" "))+facet_wrap(ModName ~ SCENARIO)
      allplot[[AreaItem]] <- plot3 
      outname <- paste0(outdir,"byRegion/",rr,"/merge/",AreaItem,"_",rr,".png")
      ggsave(plot3, file=outname, width=mergecolnum*2, height=max(1,floor(length(unique(XX$SCENARIO))/mergecolnum))*10+2,limitsize=FALSE)
      plotflag[[AreaItem]] <- nrow(XX)  
    }
    #Final energy consumption area
    pp_tfc <- plot_grid(allplot[["TFC_Ind"]],allplot[["TFC_Tra"]],allplot[["TFC_Res"]],allplot[["TFC_Com"]],ncol=2,align = "hv")
    ggsave(pp_tfc, file=paste0(outdir,"byRegion/",rr,"/merge/tfc_",rr,".png"), width=9*2, height=(floor(length(unique(allmodel$SCENARIO))/4+1)*3+2)*3,limitsize=FALSE)
  }
}

# making cross regional figure
plotXregion <-function(InputX,ii,rr){
#  for(ii in lst$varlist){
  InputX <- filter(allmodel,Region %in% rr)
  linepalettewName1 <- linepalette[1:length(unique(filter(InputX,Var==ii)$SCENARIO))]
  names(linepalettewName1) <- unique(filter(InputX,Var==ii)$SCENARIO)
  Data4Plot <- filter(InputX,Var==ii)
  miny <- min(Data4Plot$Y,2010) 
  plot.0 <- ggplot() + 
    geom_line(data=filter(Data4Plot, ModName!="Reference" & Y<=maxy),aes(x=Y, y = Value , color=SCENARIO,group=interaction(SCENARIO,ModName)),stat="identity") +
    geom_point(data=filter(Data4Plot, ModName!="Reference" & Y<=maxy),aes(x=Y, y = Value , color=SCENARIO,shape=ModName),size=1.5,fill="white") +
    MyThemeLine + scale_color_manual(values=linepalettewName1) + scale_x_continuous(breaks=seq(miny,maxy,10)) +
    scale_shape_manual(values = 1:length(unique(InputX$ModName))) +
    xlab("year") + ylab(paste0(varlist$V2.y[varlist$V1==ii],"(",varlist$V3[varlist$V1==ii],")"))  +  ggtitle(paste("Multi-regions",expression("\n"),varlist$V2.y[varlist$V1==ii],sep=" ")) +
    annotate("segment",x=miny,xend=maxy,y=0,yend=0,linetype="dashed",color="grey")+ 
    theme(legend.title=element_blank()) +facet_wrap(~Region,scales="free")
  if(length(scenariomap$SCENARIO)<20){
    plot.0 <- plot.0 +
      geom_point(data=filter(Data4Plot, ModName=="Reference"),aes(x=Y, y = Value) , color="black",shape=0,size=1.5,fill="grey") 
  }
  return(plot.0)
}
mergefigGen <- function(ii,progr){
  progr(message='merge figures')
#  for(ii in lst$varlist){
  if(nrow(filter(allmodel,Var==ii  & ModName!="Reference"))>0){
    plot.reg <- plotXregion(filter(allmodel,Region %in% R17R),ii,R17R)
    outname <- paste0(outdir,"multiReg/png/",ii,"_R17.png")
    ggsave(plot.reg, file=outname, dpi = 150, width=15, height=12,limitsize=FALSE)
  }
  if(nrow(filter(allmodel,Var==ii & Region %in% R5R & ModName!="Reference"))>0){
    plot.reg <- plotXregion(filter(allmodel,Region %in% R5R),ii,R5R)
    outname <- paste0(outdir,"multiRegR5/png/",ii,"_R5.png")
    ggsave(plot.reg, file=outname, dpi = 150, width=12, height=7.5,limitsize=FALSE)
  }
}

#function for area figure cross region
funcAreaXregionPlotGen <- function(AreaItem,progr){
  Data4plot <- allmodel_area %>% left_join(areamap,by="Var") %>% ungroup() %>% 
      filter(Class==AreaItem) %>% select(ModName,SCENARIO,Region,Ind,Y,Value,order)  %>% arrange(order)
  ModelList <- unique(as.vector(Data4plot$ModName))
  ScenarioList <- unique(as.vector(scenariomap$SCENARIO))
  #scenario model loop
  for(SC in ScenarioList){
    for(MD in ModelList){
      if(MD!="Reference"){
        XX <- Data4plot %>% filter(ModName==MD & SCENARIO==SC) %>% select(Region,Ind,Y,Value,order)
        XX2 <- Data4plot %>% filter(ModName=="Reference") %>% select(Region,Ind,Y,Value,order)  %>% arrange(order)%>% filter(Y<=2015)
        XX3 <- allmodel %>% filter(ModName==MD & SCENARIO==SC & Var %in% as.vector(areamappara$lineVar[areamappara$Class==AreaItem]) & ModName!="Reference") %>% select(Region,Var,Y,Value)
        if(nrow(XX)>0){      
          numitem <- length(as.vector(unique(Data4plot$Region))) #Get number of items
          plot1 <- funcAreaPlotSpe(XX,XX2,XX3,AreaItem)
          plot3 <- plot1 + facet_wrap( ~ Region,scales="free_y",ncol=mergecolnum) + ggtitle(paste(AreaItem,sep=" "))
          outname <- paste0(outdir,"multiReg/merge/",SC,"_",MD,"_",AreaItem,".png")
          ggsave(plot3, file=outname, width=mergecolnum*3, height=max(1,floor(numitem/mergecolnum))*5+2,limitsize=FALSE)
        }
      }
    }
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
allplot_nonleg <- as.list(nalist)
plotflag <- as.list(nalist)
names(allplot) <- nalist
names(allplot_nonleg) <- nalist
names(plotflag) <- nalist
allplotmerge <- as.list(nalist)
plotflagmerge <- as.list(nalist)
lst <- list()
if(args[8]=="global"){
  if(RegSpec==0){
    lst$region <- as.list(R17p5R)
  }else{
    lst$region <- as.list(region_load)
  }
}else{
  lst$region <-as.list(args[8])
}
lst$varlist <- as.list(as.vector(varlist$V1))
lst$R5R <- as.list(R5R)
lst$Area <- as.list(as.vector(unique(areamappara$Class)))

#Creat directories
for(rr in lst$region){
  regoutdir <- paste0(outdir,"byRegion/",rr)
  dirlist <- c(regoutdir,paste0(regoutdir,"/png"),paste0(regoutdir,"/pngdet"),paste0(regoutdir,"/ppt"),paste0(regoutdir,"/merge"))
  for(dd in dirlist){
    if(file.exists(dd)){}else{dir.create(dd)}
  }
}
ffff <- 1
if(ffff==1){
#regional figure generation execution
  exe_fig_make(lst$region,funcplotgen)
#regional area figure generation execution
  exe_fig_make(lst$region,funcAreaPlotGen)

#cross-regional figure generation execution
  if(length(Getregion)!=1){
    exe_fig_make(lst$varlist,mergefigGen)
#X regional for area figure
    print(as.vector(scenariomap$SCENARIO))
    allmodel_area <- filter(allmodel, Var %in% as.vector(areamap$Var)) 
    exe_fig_make(lst$Area,funcAreaXregionPlotGen)
  }
}

#regional Decomposition figure generation execution
#Decomposition analysis data load

if(Iterationflag>0){
  symDim <- 6
  attr(allmodel, "symName") <- "allmodel"
  lst3 <- wgdx.reshape(allmodel,symDim)
  wgdx.lst(gdxName = paste0(outdir,"data/allcombine.gdx"),lst3)
  write.csv(x = as.vector(unique(allmodel$SCENARIO)), row.names = FALSE,file = paste0(outdir,"/data/scenario.csv"))
  write.csv(x = as.vector(unique(allmodel$Region)), row.names = FALSE,file = paste0(outdir,"/data/region.csv"))

  system(paste0("gams analysis.gms --outdir=",outdir," --Region=",args[8]))
  source("Discrepancy.R")
}

decompositionflag <- 1
if(decompositionflag>0){
  Decom1 <- rgdx.param(paste0(dirCGEoutput,'../../',args[4],'/gdx/analysis.gdx'),'Loss_dcp_gdp' ) %>% rename("value"=Loss_dcp_gdp,"Sector"=SCO2_S,"Element"=decele) 
  scenariomap_load2 <- read.table(paste0(dirCGEoutput,'../../',args[4],'/txt/scenario_list.txt'), header=F,stringsAsFactors=F)
  scenariomap2 <- cbind(scenariomap_load2,scenariomap_load2,"CGE")
  names(scenariomap2) <- c("SCENARIO","Name","ModName")
  flabel <- c("change in % of GDP","sectors")
  #Function for Decomposition
  funcDecGen <- function(rr,progr){
    progr(message='region figures')
#  for(rr in as.vector(lst$region)){
    Decom2 <- Decom1 %>% filter(SCENARIO %in% scenariomap2$SCENARIO & Y %in% c(2030,2050,2100) & Sector %in% c("IND","SER","PWR","OEN","TRS","AGR") & R==rr) 
    if(nrow(filter(Decom2,Element %in% c("fd_output","output_va","va","residual1")))>0){

    plotdec <- ggplot() + geom_bar(data=filter(Decom2,Element %in% c("fd_output","output_va","va","residual1")),aes(x=Sector, y = value*100 , fill=Element), stat="identity") +
      geom_point(data=filter(Decom2,Element %in% c("fd")),aes(x=Sector, y = value*100 ),color="black", stat="identity") +
      ylab(flabel[1]) + xlab(flabel[2]) +labs(fill="") +
      guides(fill=guide_legend(reverse=TRUE)) + 
      MyThemeLine + theme(legend.position="bottom", text=element_text(size=12))+
      guides(fill=guide_legend(ncol=5))+ggtitle(paste0(rr,expression("\n")," decomposition"))+
      facet_grid(Y~SCENARIO,scales="free_x") + annotate("segment",x=0,xend=6,y=0,yend=0,linetype="dashed",color="grey")
    outname <- paste0(outdir,"byRegion/",rr,"/merge/",rr,"_decomp.png")
    ggsave(plotdec, file=outname, width=floor(length(unique(Decom2$SCENARIO))/2+1)*4, height=10,limitsize=FALSE)    
    }
  }
  exe_fig_make(lst$region,funcDecGen)
}
