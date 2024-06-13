### functions to load data and preprocess it for Stan model ###

source("R/HelperFunctions.R")
LoadData <- function(DataSet, illumination_range = c(0, Inf)) {
    DataContainer <- list()
    illumination_range[2] <- illumination_range[2] + 1
    suppressMessages(
        {
    if(DataSet=="includeRetinals") {
        DataContainer <- list(
            AssayData = list(
                readxl::read_xlsx(path = "Data/summary JellyRho new.xlsx", skip = 10, col_names = F, n_max = 12),
                readxl::read_xlsx(path = "Data/summary bovRho 12.02.2024.xlsx", skip = 11, col_names = F, n_max = 12),
                readxl::read_xlsx(path = "Data/summary bovRho retinal isomers 10s.xlsx", skip = 10, col_names = F, n_max = 12),
                readxl::read_xlsx(path = "Data/summary JellyRho retinal isomers 10s.xlsx", skip = 10, col_names = F, n_max = 12)
            ),
            LoadMetaData = list(
                readxl::read_xlsx(path = "Data/summary JellyRho new.xlsx", col_names = F, n_max = 9),
                readxl::read_xlsx(path = "Data/summary bovRho 12.02.2024.xlsx", col_names = F, n_max = 9, skip = 1),
                readxl::read_xlsx(path = "Data/summary bovRho retinal isomers 10s.xlsx", col_names = F, n_max = 9),
                readxl::read_xlsx(path = "Data/summary JellyRho retinal isomers 10s.xlsx", col_names = F, n_max = 9)
            )
        )
    } else if (DataSet=="fullATR" | DataSet=="fullATR10s" | DataSet=="fullATRTrimmed" | DataSet=="fullATRBov10sOnly" | DataSet=="fullATRcleaned_bov") {
        DataContainer <- list(
            AssayData = list(
                readxl::read_xlsx(path = "Data/summary JellyRho new.xlsx", skip = 10, col_names = F, n_max = 12),
                readxl::read_xlsx(path = "Data/summary bovRho 12.02.2024.xlsx", skip = 11, col_names = F, n_max = 12)
            ),
            LoadMetaData = list(
                readxl::read_xlsx(path = "Data/summary JellyRho new.xlsx", col_names = F, n_max = 9),
                readxl::read_xlsx(path = "Data/summary bovRho 12.02.2024.xlsx", col_names = F, n_max = 9, skip = 1)
            )
        )
    } else if(DataSet=="jellyOnly") {
        DataContainer <- list(
            AssayData = list(
                readxl::read_xlsx(path = "Data/summary JellyRho new.xlsx", skip = 10, col_names = F, n_max = 12),
                readxl::read_xlsx(path = "Data/summary JellyRho retinal isomers 10s.xlsx", skip = 10, col_names = F, n_max = 12)
            ),
            LoadMetaData = list(
                readxl::read_xlsx(path = "Data/summary JellyRho new.xlsx", col_names = F, n_max = 9),
                readxl::read_xlsx(path = "Data/summary JellyRho retinal isomers 10s.xlsx", col_names = F, n_max = 9)
            )
        )
        
    } else if(DataSet=="jellyRetinalsOnly") {
        DataContainer <- list(
            AssayData = list(
                readxl::read_xlsx(path = "Data/summary JellyRho retinal isomers 10s.xlsx", skip = 10, col_names = F, n_max = 12)
            ),
            LoadMetaData = list(
                readxl::read_xlsx(path = "Data/summary JellyRho retinal isomers 10s.xlsx", col_names = F, n_max = 9)
            )
        )
        
    } else {
        stop("DataSet not recognized")
        }

suppressWarnings({
    DataContainer <- lapply(1:length(DataContainer$LoadMetaData), function(index) {
    AssayDataConstruction(DataContainer$AssayData[[index]], DataContainer$LoadMetaData[[index]])})})

AssayDataExpInfo <- rbindlist(lapply(DataContainer, function(x) {
  x$AssayDataExpInfo
}), use.names = T, fill = T)

AssayData <- rbindlist(lapply(DataContainer, function(x) {
  x$AssayData
}), use.names = T, fill = T)
rm(DataContainer)
}
)
AssayData <- AssayData[illumination_time %between% illumination_range,]
if(DataSet=="fullATR10s") {
    AssayData <- AssayData[illumination_time==10000,]
} else if(DataSet=="fullATRTrimmed") {
    AssayData <- AssayData[illumination_time %between% illumination_range,]
} else if(DataSet=="fullATRBov10sOnly") {
    AssayData <- AssayData[(construct=="JellyRho" & illumination_time <60000)|(construct=="bovRho" & illumination_time==10000),]
} else if (DataSet == "includeRetinals") {
    AssayData <- AssayData[illumination_time %between% illumination_range,]
} else if (DataSet == "fullATRcleaned_bov") {
    AssayData <- AssayData[!(!(illumination_time %between% c(1000,30000)) & construct=="bovRho"),]
} else if (DataSet == "jellyOnly" | DataSet == "jellyRetinalsOnly") {
} else {
    stop("DataSet not recognized")
}

AssayData[,Plate:=.GRP,by=.(date, plate)]
### reconstruct AssayDataExpInfo after subsetting ###

LightActivity <- AssayData[Wavelength==1,Activity,]
AssayDataControl <- AssayData[Wavelength==0,]
AssayDataControl[,Light_Activity:=LightActivity,]
rm(LightActivity)
AssayData <- AssayData[Wavelength>1,]
AssayDataExpInfo <- AssayData[
    ,.(count=.N,
       ExpID=paste(sep = "_",
                          construct,
                          date,
                          retinal,
                          culture,
                          plate)),
    by=.(construct, date, retinal, culture, plate, illumination_time, Plate)]
setnames(x = AssayDataControl, old = "Activity", new = "Dark_Activity")

AssayData[,ConstructRetinal:=factor(interaction(construct, retinal)),]
AssayData[,colnames(AssayData) := lapply(.SD, FUN = function(x) if(is.character(x)) as.factor(x) else x),]

AssayDataExpInfo[,ConstructRetinal:=factor(interaction(construct, retinal)),]
AssayDataExpInfo[,colnames(AssayDataExpInfo) := lapply(.SD, FUN = function(x) if(is.character(x)) as.factor(x) else x),]
AssayDataControl[,ConstructRetinal:=factor(interaction(construct, retinal)),]
AssayDataControl[,colnames(AssayDataControl) := lapply(.SD, FUN = function(x) if(is.character(x)) as.factor(x) else x),]
AssayDataControl[Dark_Activity<=0,Dark_Activity:=1e-12,][Light_Activity<=0,Light_Activity:=1e-12,]

AssayData <- merge(AssayData, AssayDataControl[,.(Dark_Activity=mean(Dark_Activity)),by=.(Plate,ConstructRetinal)],
                   by=c("Plate", "ConstructRetinal"), all.x=T)
AssayData[,CorrectedActivity:=Activity-Dark_Activity,]

list(AssayData = AssayData,
     AssayDataExpInfo = AssayDataExpInfo,
     AssayDataControl = AssayDataControl)
}
#### Load and preprocess data ####

ConstructStanData <- function(AssayList) {
### set up data for stan model
ModelMatrix <- as.matrix(
  model.matrix(data = AssayList$AssayData, object = Activity ~ 1 + ConstructRetinal))

DarkModelMatrix <- as.matrix(
  model.matrix(data = AssayList$AssayDataControl, object = Dark_Activity ~ 1 + ConstructRetinal))

list(N=AssayList$AssayData[,.N,],
                 priorOnly=0,
                 GroupN = AssayList$AssayData[,.N,by=ConstructRetinal][,.N,],
                 Group = as.matrix(model.matrix(object = ~ 1 + ConstructRetinal, data = AssayList$AssayData)[,-1]),
                 ExperimentN=AssayList$AssayData[,.N,by=ExpID][,.N,],
                 Experiment=as.integer(factor(AssayList$AssayData$ExpID)),
                 PlateN=AssayList$AssayData[,.N,by=Plate][,.N,],
                 Plate=AssayList$AssayData$Plate,
                 ExperimentConstruct=AssayList$AssayData[,.(ConstructRetinal=as.integer(factor(ConstructRetinal)), ExpID),][,.(ConstructRetinal=unique(ConstructRetinal)),by=ExpID][,ConstructRetinal,],
                 Dark_Group=as.matrix(model.matrix(object = ~ 1 + ConstructRetinal, data = AssayList$AssayDataControl)[,-1]),
                 DarkGroupN=AssayList$AssayDataControl[,.N,by=ConstructRetinal][,.N,],
                 DarkPlateID=AssayList$AssayDataControl$Plate,
                 Activity=AssayList$AssayData$Activity,
                 lambda=AssayList$AssayData$Wavelength,
                 illumination_time=AssayList$AssayData$illumination_time,
                 DarkN=AssayList$AssayDataControl[,.N,],
                 DarkExperimentID=as.integer(factor(AssayList$AssayDataControl$ExpID)),
                 DarkActivity=AssayList$AssayDataControl$Dark_Activity)
}

