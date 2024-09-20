library(readxl)
library(readr)

#Two AR6 files should be located in data directory if you would like to process AR6 original data.
#The original data is too large for git management.
AR6_data_file <- '../data/AR6_Scenarios_Database_World_v1.1.csv'
AR6reg_data_file <- '../data/AR6_Scenarios_Database_R5_regions_v1.1.csv'
AR6_meta_file <- '../data/AR6_Scenarios_Database_metadata_indicators_v1.1.xlsx'
varalllist <- read.table(VarListPath, sep="\t",header=F, stringsAsFactors=F)
cmapping <- data.frame(Category=c("C1","C2","C3","C4","C5","C6","C7","C8"),Category4=c("C1-C2","C1-C2","C3-C4","C3-C4","C5-C6","C5-C6","C7-C8","C7-C8"))
CategorySet <- c("C1","C2","C3","C4","C5","C6","C7","C8")

#Loading data
AR6_global_meta <- read_xlsx(path=AR6_meta_file,sheet='meta_Ch3vetted_withclimate') %>% 
  select(Model,Scenario,Category,IMP_marker)

AR6_global_load <- read_csv(file=AR6_data_file) %>% 
  pivot_longer(cols=-c(Model,Scenario,Region,Variable,Unit),names_to='Y',values_to='Value') %>% 
  filter(!is.na(Value)) %>% 
  mutate(Y=as.numeric(Y)) %>% 
  inner_join(AR6_global_meta,by=c('Model','Scenario')) %>%
  rename(SCENARIO=Scenario) %>% 
  inner_join(varalllist %>% rename(Var=V1,Variable=V2)) %>% select(-V3,-V4,-IMP_marker)

AR6_regional_load <- read_csv(file=AR6reg_data_file) %>% 
  pivot_longer(cols=-c(Model,Scenario,Region,Variable,Unit),names_to='Y',values_to='Value') %>% 
  filter(!is.na(Value)) %>% 
  mutate(Y=as.numeric(Y)) %>% 
  inner_join(AR6_global_meta,by=c('Model','Scenario')) %>%
  rename(SCENARIO=Scenario) %>% 
  inner_join(varalllist %>% rename(Var=V1,Variable=V2)) %>% select(-V3,-V4,-IMP_marker)

#Cleaning the data and extracting each classification 
AR6DB.nat <- rbind(AR6_global_load,filter(AR6_regional_load,Region!="World"))  %>% select(Model,SCENARIO,Region,Variable,Category,Var,Unit,Y,Value)
AR6DB.rcate <- left_join(AR6DB.nat,cmapping,by="Category") %>% select(-Category) %>%rename(Category=Category4) %>% select(Model,SCENARIO,Region,Variable,Category,Var,Unit,Y,Value)
AR6DB <- rbind(AR6DB.nat,AR6DB.rcate) %>% filter(Var %in% varlist$V1)
AR6_global_load <- 0
AR6_regional_load <- 0

AR6DB_spread <- AR6DB %>% spread(key=Y,value=Value)
AR6Reg <- as.vector(unique(AR6DB$Region))
AR6Var <- as.vector(unique(AR6DB$Var))
AR6Cat <- as.vector(unique(AR6DB$Category))
AR6Varariable <- as.vector(unique(AR6DB.nat$Variable))

#Extract median, min and max for each category, region and variables
w <- 0
for(j in AR6Reg){
    for(k in AR6Var){
    for(i in AR6Cat){ 
      if(nrow(filter(AR6DB,Region==j,Var==k,Category==i))>=1){
        w <- w +1
        AR6DBtmp <- AR6DB_spread %>% filter(Region==j,Var==k,Category==i) %>% select(-Model,-SCENARIO,-Unit,-Region,-Var,-Category,-Variable)
        AR6DB.ind <- rbind(AR6DBtmp %>% summarise(across(everything(), ~ max(ifelse(is.infinite(.x), NA, .x), na.rm = TRUE)))%>% mutate(ind="max",Region=j,Var=k,Category=i),
                           AR6DBtmp %>% summarise(across(everything(), ~ min(ifelse(is.infinite(.x), NA, .x), na.rm = TRUE)))%>% mutate(ind="min",Region=j,Var=k,Category=i),
                           AR6DBtmp %>% summarise(across(everything(), ~ median(ifelse(is.infinite(.x), NA, .x), na.rm = TRUE)))%>% mutate(ind="med",Region=j,Var=k,Category=i),
                           AR6DBtmp %>% summarise(across(everything(), ~ quantile(ifelse(is.infinite(.x), NA, .x), na.rm = TRUE, probs = 0.9))) %>% mutate(ind="90p",Region=j,Var=k,Category=i),
                           AR6DBtmp %>% summarise(across(everything(), ~ quantile(ifelse(is.infinite(.x), NA, .x), na.rm = TRUE, probs = 0.1))) %>% mutate(ind="10p",Region=j,Var=k,Category=i)
        )
        if(w==1){
          AR6DBall <- AR6DB.ind
        } else{
          AR6DBall <- rbind(AR6DBall,AR6DB.ind)
        }
      }
    }
  }
}
AR6DBall0 <- AR6DBall  
AR6DBall_inf <- AR6DBall 
AR6DBall_inf[] <- lapply(AR6DBall_inf, function(x) ifelse(is.infinite(x), NA, x))
AR6DBall_inf2 <- pivot_longer(AR6DBall_inf,cols=c("2010","2020","2030","2040","2050","2060","2070","2080","2090","2100"),names_to = "Y", values_to = "Value") %>%
   select(Region,Var,Category,Y,ind,Value) %>% pivot_wider(names_from = ind, values_from = Value) %>% mutate(Y=as.numeric(Y))
write.csv(x = AR6DBall_inf2, row.names = FALSE,file = paste0("../data/AR6Slected.csv"))
#AR6DBIndload <- AR6DBall_inf2 %>% filter(Category %in% CategorySub)


