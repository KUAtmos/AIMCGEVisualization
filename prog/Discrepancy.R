

Varass <- c("Fin_Ene", "Fin_Ene_Ele_Heat", "Fin_Ene_Liq_and_Gas", "Fin_Ene_SolidsCoa", "Fin_Ene_SolidsBio", "Fin_Ene_Hyd", "Fin_Ene_Ind", "Fin_Ene_Com", "Fin_Ene_Res", "Fin_Ene_Tra",
            "Sec_Ene_Ele_Heat","Sec_Ene_Ele_Heat_Coa","Sec_Ene_Ele_Heat_Gas","Sec_Ene_Ele_Heat_Nuc","Sec_Ene_Ele_Hyp","Sec_Ene_Ele_Solar","Sec_Ene_Ele_Win","Sec_Ene_Ele_Heat_Bio","Sec_Ene_Ele_Heat_Geo",
            "Sec_Ene_Hyd","Sec_Ene_Hyd_Gas","Sec_Ene_Hyd_Bio","Sec_Ene_Hyd_Ele"
#            ,"Fin_Ene_Ind_Ele_Heat", "Fin_Ene_Ind_Liq_and_Gas", "Fin_Ene_Ind_SolidsCoa", "Fin_Ene_Ind_SolidsBio", "Fin_Ene_Ind_Hyd", "Fin_Ene_Com_Ele_Heat", "Fin_Ene_Com_Liq_and_Gas", "Fin_Ene_Com_SolidsCoa", "Fin_Ene_Com_SolidsBio", "Fin_Ene_Com_Hyd", "Fin_Ene_Res_Ele_Heat", "Fin_Ene_Res_Liq_and_Gas", "Fin_Ene_Res_SolidsCoa", "Fin_Ene_Res_SolidsBio", "Fin_Ene_Res_Hyd", "Fin_Ene_Tra_Ele", "Fin_Ene_Tra_Liq_and_Gas", "Fin_Ene_Tra_Liq_Bio", "Fin_Ene_Tra_Liq_Oil", "Fin_Ene_Tra_Hyd"
)
ScenarioNameMap <- data.frame(c("Baseline","Mitigation"),c("Baseline","Mitigation"))
names(ScenarioNameMap) <- c("target","origin")

# The statistics in analysis.gdx compare the iteration results (AIMHub0, AIMHub1, ..., AIMTech0, ...).
# When the merged data has no iteration results (for example the merged IAMC data without the iteration models) the statistics are empty.
# An empty plot with facet_wrap stops R (Faceting variables must have at least one value), so the figures are skipped in that case.

#COnvergence check figures and list file
Conv_load <- rgdx.param(paste0(outdir,'data/analysis.gdx'),'ConvStat') %>% rename("Value"=ConvStat) %>% filter(Var %in% Varass) %>%
  left_join(varalllist %>% rename(Var=V1) %>% select(-V3,-V4)) %>% rename (VarCode=Var) %>% rename(Var=V2) %>% left_join(ScenarioNameMap %>% rename(SCENARIO=origin))
Conv_plot <- filter(Conv_load, Indi=="MAEStand" & ModSpe=="AIMHub")
if(nrow(Conv_plot)>0){
  plot.0 <- ggplot() +
    geom_line(data=Conv_plot,aes(x=Ite, y = Value * 100, color=Var,group=Var),stat="identity") +
    MyThemeLine +  xlab("Iteration") + ylab("Discrepancy (%)") +
    facet_wrap(Region ~ target) + guides(color = guide_legend(ncol = 1, byrow = TRUE, title=""))
  outname <- paste0(outdir,"misc/R5Conv.png")
  ggsave(plot.0, file=outname, width=15, height=15,limitsize=FALSE)
}else{
  print("Discrepancy.R: no iteration results for the convergence figure (R5Conv.png skipped)")
}

Exceedlist <- filter(Conv_load, Indi=="MAEStand" & ModSpe=="AIMHub" & Ite=="i3" & Value >= 0.05)
write.csv(x = Exceedlist, file = paste0(outdir,"misc/Exceedlist_conv.csv"))

#Difference figures
Diff_load <- rgdx.param(paste0(outdir,'data/analysis.gdx'),'Stat') %>% rename("Value"=Stat) %>% filter(Var %in% Varass) %>%
  left_join(varalllist %>% rename(Var=V1) %>% select(-V3,-V4)) %>% rename (VarCode=Var) %>% rename(Var=V2) %>% left_join(ScenarioNameMap %>% rename(SCENARIO=origin))
Diff_plot <- filter(Diff_load, Indi=="MAEStand")
if(nrow(Diff_plot)>0){
  plot.0 <- ggplot() +
    geom_line(data=Diff_plot,aes(x=Ite, y = Value * 100, color=Var,group=Var),stat="identity") +
    MyThemeLine +  xlab("Iteration") + ylab("Discrepancy (%)") +
    facet_wrap(Region ~ target) + guides(color = guide_legend(ncol = 1, byrow = TRUE, title=""))
  outname <- paste0(outdir,"misc/R5ModelDif.png")
  ggsave(plot.0, file=outname, width=15, height=15,limitsize=FALSE)
}else{
  print("Discrepancy.R: no iteration results for the model difference figure (R5ModelDif.png skipped)")
}

Exceedlist <- filter(Diff_load, Indi=="MAEStand" & Ite=="i3" & Value >= 0.1)
write.csv(x = Exceedlist, file = paste0(outdir,"misc/Exceedlist_dif.csv"))
