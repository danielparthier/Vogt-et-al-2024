setwd(dir = projectDir)
source("R/LoadDataFunctions.R")

n_samples_per_chain <- ceiling(n_samples/n_chains)

AssayList <- LoadData("fullATR")
stanData <- ConstructStanData(AssayList)

dirPath <- "Output/20240611_ATR/"
if(!dir.exists(paths = dirPath)) {
  dir.create(path = dirPath)
}

ModelRun(stanData = stanData,
         n_chains = n_chains,
         n_samples = n_samples_per_chain,
         output_path = dirPath,
         prior_only = F)

stanDataPrior <- stanData
stanDataPrior$priorOnly <- 1

ModelRun(stanData = stanData,
         n_chains = n_chains,
         n_samples = n_samples_per_chain,
         output_path = dirPath,
         prior_only = T)

saveRDS(object = AssayList, file = paste0(dirPath, "AssayList.rds"))
