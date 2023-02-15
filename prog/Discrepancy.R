

Varass <- c("Fin_Ene", "Fin_Ene_Ele_Heat", "Fin_Ene_Liq_and_Gas", "Fin_Ene_SolidsCoa", "Fin_Ene_SolidsBio", "Fin_Ene_Hyd", "Fin_Ene_Ind", "Fin_Ene_Com", "Fin_Ene_Res", "Fin_Ene_Tra",
            "Sec_Ene_Ele_Heat","Sec_Ene_Ele_Heat_Coa","Sec_Ene_Ele_Heat_Gas","Sec_Ene_Ele_Heat_Nuc","Sec_Ene_Ele_Hyd","Sec_Ene_Ele_Solar","Sec_Ene_Ele_Win","Sec_Ene_Ele_Heat_Bio","Sec_Ene_Ele_Heat_Geo",
            "Sec_Ene_Hyd","Sec_Ene_Hyd_Gas","Sec_Ene_Hyd_Bio","Sec_Ene_Hyd_Ele",
            "Fin_Ene_Ind_Ele_Heat", "Fin_Ene_Ind_Liq_and_Gas", "Fin_Ene_Ind_SolidsCoa", "Fin_Ene_Ind_SolidsBio", "Fin_Ene_Ind_Hyd", "Fin_Ene_Com_Ele_Heat", "Fin_Ene_Com_Liq_and_Gas", "Fin_Ene_Com_SolidsCoa", "Fin_Ene_Com_SolidsBio", "Fin_Ene_Com_Hyd", "Fin_Ene_Res_Ele_Heat", "Fin_Ene_Res_Liq_and_Gas", "Fin_Ene_Res_SolidsCoa", "Fin_Ene_Res_SolidsBio", "Fin_Ene_Res_Hyd", "Fin_Ene_Tra_Ele", "Fin_Ene_Tra_Liq_and_Gas", "Fin_Ene_Tra_Liq_Bio", "Fin_Ene_Tra_Liq_Oil", "Fin_Ene_Tra_Hyd"
)
ScenarioNameMap <- data.frame(c("Baseline","Mitigation"),c("Baseline","0500CDAC"))
names(ScenarioNameMap) <- c("target","origin")

#COnvergence check figures and list file
Conv_load <- rgdx.param('../output/data/analysis.gdx','ConvStat') %>% rename("Value"=ConvStat) %>% filter(Var %in% Varass) %>% 
  left_join(varalllist %>% rename(Var=V1) %>% select(-V3,-V4)) %>% rename (VarCode=Var) %>% rename(Var=V2) %>% left_join(ScenarioNameMap %>% rename(SCENARIO=origin))
plot.0 <- ggplot() + 
  geom_line(data=filter(Conv_load, Indi=="MAEStand" & ModSpe=="CGE"),aes(x=Ite, y = Value * 100, color=Var,group=Var),stat="identity") +
  MyThemeLine +  xlab("Iteration") + ylab("Discrepancy (%)") +
  facet_wrap(Region ~ target) + guides(color = guide_legend(ncol = 1, byrow = TRUE, title=""))
outname <- paste0(outdir,"misc/R5Conv.png")
ggsave(plot.0, file=outname, width=15, height=15,limitsize=FALSE)

Exceedlist <- filter(Conv_load, Indi=="MAEStand" & ModSpe=="CGE" & Ite=="i3" & Value >= 0.05)
write.csv(x = Exceedlist, file = "../output/misc/Exceedlist_conv.csv")

#Difference figures
Diff_load <- rgdx.param('../output/data/analysis.gdx','Stat') %>% rename("Value"=Stat) %>% filter(Var %in% Varass) %>% 
  left_join(varalllist %>% rename(Var=V1) %>% select(-V3,-V4)) %>% rename (VarCode=Var) %>% rename(Var=V2) %>% left_join(ScenarioNameMap %>% rename(SCENARIO=origin))
plot.0 <- ggplot() + 
  geom_line(data=filter(Diff_load, Indi=="MAEStand"),aes(x=Ite, y = Value * 100, color=Var,group=Var),stat="identity") +
  MyThemeLine +  xlab("Iteration") + ylab("Discrepancy (%)") +
  facet_wrap(Region ~ target) + guides(color = guide_legend(ncol = 1, byrow = TRUE, title=""))
outname <- paste0(outdir,"misc/R5ModelDif.png")
ggsave(plot.0, file=outname, width=15, height=15,limitsize=FALSE)

Exceedlist <- filter(Diff_load, Indi=="MAEStand" & Ite=="i3" & Value >= 0.05)
write.csv(x = Exceedlist, file = "../output/misc/Exceedlist_dif.csv")
