libloadlist <- c("gdxrrw","ggplot2","dplyr","reshape2","tidyr","maps","grid","RColorBrewer","cowplot","hms","purrr","furrr","progressr","readr","readxl")
for(j in libloadlist){
  eval(parse(text=paste0("library(",j,")")))
}

#Two AR6 files should be located in data directory if you would like to process AR6 original data.
#The original data is too large for git management.
outdir <- "../../../../output/global/" 
AR6_data_file <- '../data/AR6_Scenarios_Database_World_v1.1.csv'
AR6reg_data_file <- '../data/AR6_Scenarios_Database_R5_regions_v1.1.csv'
AR6_meta_file <- '../data/AR6_Scenarios_Database_metadata_indicators_v1.1.xlsx'
varalllist <- read.table(paste0(outdir,"../include_code/IAMCTemp/all_list.txt"), sep="\t",header=F, stringsAsFactors=F)
cmapping <- data.frame(Category=c("C1","C2","C3","C4","C5","C6","C7","C8"),Category4=c("C1-C2","C1-C2","C3-C4","C3-C4","C5-C6","C5-C6","C7-C8","C7-C8"))
CategorySet <- c("C1","C2","C3","C4","C5","C6","C7","C8")
args <- commandArgs(trailingOnly = TRUE)
default_args <- c("/opt/gams/gams37.1_linux_x64_64_sfx", min(floor(availableCores()/2),24),"off","global/global_17","2","0","IAMCTemplate_Iteoff_global","global","on","scenarioMIP")
default_flg <- is.na(args[1:10])
args[default_flg] <- default_args[default_flg]
gams_sys_dir <- as.character(args[1])
igdx(gams_sys_dir)
sizememory <- 10000*1024^2 
options(future.globals.maxSize= sizememory)
threadsnum <-  as.numeric(args[2])

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
AR6DB <- rbind(AR6DB.nat,AR6DB.rcate) %>% filter(Var %in% varalllist$V1)

AR6DB_spread <- AR6DB %>% spread(key=Y,value=Value)
AR6Reg <- as.vector(unique(AR6DB$Region))
AR6Var <- as.vector(unique(AR6DB$Var))
AR6Cat <- as.vector(unique(AR6DB$Category))
AR6Varariable <- as.vector(unique(AR6DB.nat$Variable))

AR6_global_load <- 0
AR6_regional_load <- 0
AR6DB.nat <- 0
AR6DB.rcate <- 0
AR6DB <- 0

lst <- list()
lst$var <- as.list(AR6Var)
#lst$var <- head(lst$var, 10)
parallelmode <- 1
#Extract median, min and max for each category, region and variables
funcplotgen <- function(k,progr){
w <- 0
AR6DBall <- NULL
  progr(message='region figures')
#  for(k in AR6Var){
    for(j in AR6Reg){
    for(i in AR6Cat){ 
      if(nrow(filter(AR6DB_spread,Region==j,Var==k,Category==i))>=1){
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
return(AR6DBall)
}
# progress bar parallell computation
plan(multisession, workers = threadsnum)
with_progress({
  progr <- progressor(along = lst$var)  # 最初の10個に限定
  AR6DBall <- future_map(lst$var, 
                         ~ funcplotgen(.x, progr = progr),
                         .options = furrr_options(seed = TRUE))
})

AR6DBall0 <- AR6DBall  
AR6DBall_inf <- AR6DBall 
names(AR6DBall) <- as.vector(lst$var)
AR6DBall_inf <- bind_rows(AR6DBall)
AR6DBall_inf[] <- lapply(AR6DBall_inf, function(x) ifelse(is.infinite(x), NA, x))
AR6DBall_inf2 <- pivot_longer(AR6DBall_inf,cols=c("2010","2020","2030","2040","2050","2060","2070","2080","2090","2100"),names_to = "Y", values_to = "Value") %>%
   select(Region,Var,Category,Y,ind,Value) %>% pivot_wider(names_from = ind, values_from = Value) %>% mutate(Y=as.numeric(Y))
write.csv(x = AR6DBall_inf2, row.names = FALSE,file = paste0("../data/AR6Slected.csv"))
#AR6DBIndload <- AR6DBall_inf2 %>% filter(Category %in% CategorySub)

#ここからは試作
function(){
variable_sets <- read.table("../data/AR6_calc_set.txt", sep=",",header=F, stringsAsFactors=F)
colnames(variable_sets) <- c("Variable1", "Variable2", "Ratios")

# 比率を計算する関数の定義
calculate_ratios <- function(df, variable1, variable2, ratio_name) {
  #  df %>%
  a <-  AR6DB %>%
    filter(Var %in% c(variable1, variable2)) %>%
    pivot_wider(names_from = Var, values_from = Value) %>%
    mutate(!!ratio_name := !!sym(variable1) / !!sym(variable2)) %>%
    select(Model,SCENARIO,Region,Variable,Category,Var,Unit,Y, !!ratio_name)
}

# すべての変数セットに対して処理を行い、結果をリストに格納
results_list <- lapply(1:nrow(variable_sets), function(i) {
  variable_set <- variable_sets[i, ]
  calculate_ratios(AR6DB, variable_set$Variable1, variable_set$Variable2, variable_set$Ratios)
})

# 結果を1つのデータフレームに統合
AR6DB <- bind_rows(results_list)
}
