library(openxlsx)
args <- commandArgs(trailingOnly = TRUE)
#arguments are: 5:submodule switch, 7: GDX file name, 8: Model name
default_args <- c("1", "global_17_IAMC")   # Default value but gams path should be modified if GUI based R is used
args <- commandArgs(trailingOnly = TRUE)
default_flg <- is.na(args[1:10])
args[default_flg] <- default_args[default_flg]

Hub_IntTool_Flag <- as.numeric(args[1])
filename <- args[2] # filename should be "global_17","CHN","JPN"....

if(Hub_IntTool_Flag==1){
  diroutput <- paste0("../output/iiasa_database/txt/")  # directory where the CGE output is located 
}else if(Hub_IntTool_Flag==2){
  diroutput <- paste0("../../../../../../output/iamc/") #Output directory 
}

# Read CSV file and unload as xlsx
csv_data <- read.csv(paste0(diroutput,filename,".csv"), check.names = FALSE)
current_date <- format(Sys.Date(), "%Y%m%d")  # date for YYYYMMDD format
output_file <- paste0(diroutput,"IAMC", current_date, ".xlsx")
write.xlsx(csv_data, file=output_file, sheetName="data")
