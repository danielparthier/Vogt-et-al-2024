# check and install required packages
requiredPackages <- c("data.table", "ggplot2", "ggdist", "readxl", "cmdstanr", "bayestestR", "posterior", "doFuture", "rmarkdown", "knitr")
invisible(lapply(requiredPackages,
                 FUN = function(package) {
                   if(!require(package, character.only = T)) {
                     if(package == "cmdstanr") {
                       message("Downloading and installing cmdstanr")
                       install.packages("cmdstanr", repos = c('https://stan-dev.r-universe.dev', getOption("repos")))
                       message("Please visit https://mc-stan.org/cmdstanr/articles/cmdstanr.html check for how to finish installing cmdstanr")
                     } else {
                       message(paste("Downloading and installing", package))
                       install.packages(package)
                     }
                   }
                   })
          )
source("R/settings.R")

LoadFiles <- function(FolderPath, FileString, ConstructString) {
  FileDirectory <- paste(sep = "/", FolderPath, FileString)
  AssayData <- suppressWarnings(suppressMessages(rbindlist(lapply(X = seq_along(FileDirectory), FUN = function(x) {
    AssayData <- data.table(read_excel(FileDirectory[x], skip = 8, col_names = F, n_max = 10))
    DatesPlate <- data.table(read_excel(FileDirectory[x], skip = 2, col_names = F, n_max = 6))
    setnames(x = DatesPlate, new = c("first", 1:(ncol(DatesPlate)-1)))
    DatesPlate$first <- c("date", "plate", "row", "type", "timing")
    setnames(x = AssayData, new = c("Wavelength", 1:(ncol(DatesPlate)-1)))
    
    AssayData <- melt.data.table(data = AssayData, id.vars = "Wavelength", variable.name = "Col", value.name = "Activity",na.rm = F)
    DatesPlate <- melt.data.table(data = DatesPlate, id.vars = "first", variable.name = "Col", na.rm = F)
    setkey(x = AssayData, "Col")
    setkey(x = DatesPlate, "Col")
    
    AssayData <- merge.data.table(x = AssayData, y = DatesPlate[first=="date", .(Col, Date=value)], by = "Col")
    AssayData <- merge.data.table(x = AssayData, y = DatesPlate[first=="plate", .(Col, plate=value)], by = "Col")
    AssayData <- merge.data.table(x = AssayData, y = DatesPlate[first=="timing", .(Col, timing=strtoi(value))], by = "Col")
    AssayData <- merge.data.table(x = AssayData, y = DatesPlate[first=="type", .(Col, type=value)], by = "Col")
    AssayData <- AssayData[!is.na(Activity),]
    AssayData[,`:=`(ID=as.integer(factor(Col)),
                    Date=as.Date(x = strtoi(Date), origin = "1899-12-30", tryFormats = c("%Y-%m-%d", "%Y/%m/%d")),
                    Construct=ConstructString[x],
                    Col=NULL),]
  }))))
  AssayData[,`:=`(Experiment=paste(Date,plate, sep = "_"), ExperimentSub=paste(Construct,Date,plate, sep = "_")),]
  AssayData[is.na(type),type:="ATR",]
  AssayData[,`:=`(Date=factor(Date),
                  Construct=factor(Construct),
                  ExperimentSub=factor(ExperimentSub),
                  Experiment=factor(Experiment)),]
  return(AssayData)
}

AssayDataConstruction <- function(AssayData, LoadMetaData) {
  ColumnSelect <- as.vector(unlist(LoadMetaData[1,]))
  ColumnSelect[1] <- NA
  ColumnSelect <- grep(ColumnSelect, pattern = "[A-z]")
  LoadMetaData <- data.table(t(LoadMetaData))
  setnames(x = LoadMetaData, new = unlist(LoadMetaData[1,]))
  LoadMetaData <- LoadMetaData[!is.na(construct),][-1,]
  ColumnNames <- colnames(LoadMetaData)
  
  LoadMetaData <- LoadMetaData[,.(construct = ifelse(test = any(grepl(pattern = "construct", x = ColumnNames)), yes = construct, no = NA),
                                  date = ifelse(test = any(grepl(pattern = "date", x = ColumnNames)), date, no = NA),
                                  retinal = ifelse(test = any(grepl(pattern = "retinal", x = ColumnNames)), yes = retinal, no = NA),
                                  culture = ifelse(test = any(grepl(pattern = "culture", x = ColumnNames)), yes = culture, no = NA),
                                  plate = ifelse(test = any(grepl(pattern = "plate", x = ColumnNames)), yes = plate, no = NA),
                                  illumination_time = ifelse(test = any(grepl(pattern = "ill. time", x = ColumnNames, )), yes = as.numeric(`ill. time ms`), no = NA),
                                  well_n = ifelse(test = any(grepl(pattern = "well", x = ColumnNames)), `n (=well)`, no = NA)), by=rownames(LoadMetaData)]
  LoadMetaData[,plate:=gsub(pattern = "plate| ", replacement = "", x = plate),by=plate]
  LoadMetaData[,well_n:=gsub(pattern = "n| ", replacement = "", x = well_n),by=well_n]
  LoadMetaData[grepl(pattern = "jelly", x = construct, ignore.case = T),construct:="JellyOp-mScarlet"][
    grepl(pattern = "bov", x = construct, ignore.case = T),construct:="bovRho-mScarlet"]
  LoadMetaData[grepl(pattern = "atr", x = retinal, ignore.case = T),retinal:="ATR"]
  LoadMetaData[,rownames:=NULL,]
  LoadMetaData <- LoadMetaData[rep(1: .N, each = 12),]
  
  AssayData <- data.table(AssayData[,c(1,ColumnSelect-1)])
  AssayData <- melt.data.table(data = AssayData)[,c(1,3)]
  setnames(x = AssayData, new = c("Wavelength", "Activity"))
  AssayData <- data.table(LoadMetaData, AssayData)
  
  # encode black current with 0 and white current with 1
  AssayData[,Wavelength:=fifelse(Wavelength=="dark",0,fifelse(Wavelength=="white",1, strtoi(Wavelength))),]
  AssayData <- AssayData[!is.na(Activity),]
  AssayData[,`:=`(ExpID=factor(paste(sep = "_",
                                     construct, date, retinal, culture, plate))),]
  
    AssayDataExpInfo <- AssayData[
    ,.(count=.N,
       ExpID=paste(sep = "_",
                          construct,
                          date,
                          retinal,
                          culture,
                          plate)),
    by=.(construct, date, retinal, culture, plate, illumination_time)]
  return(list(AssayData = AssayData, AssayDataExpInfo = AssayDataExpInfo))
}

ExtractDensity <- function(ParameterDF) {
  rbindlist(lapply(X = 1:(dim(ParameterDF)[2]-3), FUN = function(x) {
    tmpDens <- density(ParameterDF[[x]])
    return(data.table(x=tmpDens$x, y=tmpDens$y, Parameter=names(ParameterDF)[x]))
  }))
}

InitChains <- function(Chains, stanData) {
  set.seed(seed = 20240301)
  lapply(1:Chains, FUN = function (x) {
    list(shape=rgamma(n = 1, shape = 200, rate = 200/50),
         betaBand_a=rgamma(n = 1, shape = 2000, rate = 2000/28.2),
         betaBand_b=rgamma(n = 1, shape = 2000, rate = 2000/80),
         A=rgamma(n = 1, shape = 2000, rate = 2000/70),
         B=rgamma(n = 1, shape = 2000, rate = 2000/28),
         C_raw=rnorm(n = 1, mean = 0, sd = 0.1),
         D_raw=rnorm(n = 1, mean = 0, sd = 0.1),
         c_raw=rnorm(n = 1, mean = 0, sd = 0.1),
         b_raw=rnorm(n = 1, mean = 2.5, sd = 0.1),
         a_1_raw=rnorm(n = 1, mean = 2, sd = 0.1),
         a_2_raw=rnorm(n = 1, mean = 0, sd = 0.1),
         lambda_max_raw = rnorm(n = 1, mean = 0, sd = .01),
         betaBandFactor = rgamma(n = 1, shape = 2000, rate = 2000/0.5),
         betaBandwidth = rgamma(n = 1, shape = 2000, rate = 2000/0.26),
         tau_Plate = rgamma(n = stanData$GroupN+3, shape = 1000, rate = 1000/0.1),
         L_Omega_Plate = matrix(data = rnorm(n = (stanData$GroupN+3)^2, mean = 0, sd = 0.01), nrow = stanData$GroupN+3, ncol = stanData$GroupN+3),
         z_Plate = matrix(data = rgamma(n = (stanData$GroupN+3)*stanData$PlateN, shape = 100, rate = 100/0.1), nrow = stanData$GroupN+3, ncol = stanData$PlateN),
         tau_Group = rgamma(n = 1, shape = 1000, rate = 1000/0.1),
         z_Group = rgamma(n = stanData$GroupN, shape = 100, rate = 100/0.1),
         betaBandMean = rbeta(n = 1, shape1 = 10, shape2 = 10),
         betaBandShape = rgamma(n = 1, shape = 200, rate = 200/20),
         betaMax_1 = rgamma(n = 1, shape = 4000, rate = 4000/189),
         betaMax_2 = rgamma(n = 1, shape = 3000, rate = 3000/0.315),
         bBetaBand_1_raw = rnorm(n = 1, mean = 0, sd = 0.1),
         beta_AlphaGroup_raw = rnorm(n = stanData$GroupN-1, mean = 0, sd = 0.001),
         beta_I_50_Group = rnorm(n = stanData$GroupN-1, mean = 0, sd = 0.001),
         beta_activity_max_Group = rnorm(n = stanData$GroupN-1, mean = 0, sd = 0.001),
         activity_max = rgamma(n = 1, shape = 1000, rate = 1000/0.4),
         I_50 = rgamma(n = 1, shape = 20, rate = 20/10000),
         Dark_shape = rgamma(n = 1, shape = 3, rate = 3/15),
         Dark_beta = rnorm(n = 1, mean = log(0.05), sd = 1),
         Dark_Group_beta = rnorm(n = stanData$DarkGroupN-1, mean = 0, sd = 0.25))
  })
}

ActivityOptim <- function(LambdaMax=500, x, i=1) {
    lambda_beta_max <-  BaseVariables$betaMax_1[i] + BaseVariables$betaMax_2[i] * LambdaMax

    b_betaBand <- BaseVariables$bBetaBand_1[i] + BaseVariables$betaBandwidth[i] * LambdaMax
    beta_band <- BaseVariables$betaBandFactor[i] * exp( -((x-lambda_beta_max) / b_betaBand)^2);
    max_div_lambda <- LambdaMax / x;
    
    mu_vec <- 1/(exp(BaseVariables$A[i] * (BaseVariables$a_1[i]+BaseVariables$a_2[i]*exp( - ((LambdaMax-300)^2)/11940.0) - max_div_lambda)) +
                 exp(BaseVariables$B[i] * (BaseVariables$b[i] - max_div_lambda)) +
                 exp(BaseVariables$C[i] * (BaseVariables$c[i] - max_div_lambda)) +
                 BaseVariables$D[i])
    mu_vec <- mu_vec * (exp(BaseVariables$A[i] * (BaseVariables$a_1[i]+BaseVariables$a_2[i]*exp(- ((LambdaMax-300)^2)/11940)-1)) +
                        exp(BaseVariables$B[i] * (BaseVariables$b[i]-1)) +
                        exp(BaseVariables$C[i] * (BaseVariables$c[i]-1)) +
                        BaseVariables$D[i])
    mu_vec <- mu_vec + beta_band
    mu_vec <- mu_vec / (1+BaseVariables$betaBandFactor[i] * exp( - ((LambdaMax-lambda_beta_max) / (b_betaBand))^2))
    return(mu_vec)
}

ActivitySimVec <- function(WaveVec, LambdaExp) {
  IterationN <- nrow(LambdaExp)
  rbindlist(foreach(i=lapply(X=1:(ncol(LambdaExp)-3), FUN = function(x) {
    list(lambda_alpha = as.numeric(LambdaExp[[x]]),
         ExpID=x)
  })) %dopar% {
    Exp <- sapply(X = WaveVec, FUN = function(w) {
      ActivitySim(LambdaMax = i$lambda_alpha, Wavelength = w)
    })
    colnames(Exp) <- WaveVec
    rownames(Exp) <- 1:IterationN
    Exp <- melt(data.table(Exp))
    setnames(x = Exp, old=c("variable", "value"), new=c("Wavelength", "Activity"))
    Exp[,`:=`(Iteration=rowid(Wavelength), ExpID=i$ExpID, Wavelength=rep(WaveVec,each=IterationN)),]
    message("Experiment: ", i$ExpID,"\r",appendLF=FALSE)
    flush.console()
    return(Exp)
  })
}

ActivitySim <- function(LambdaMax, Wavelength) {
  lambda_beta_max_exp <-  BaseVariables$betaMax_1
  + BaseVariables$betaMax_2 * LambdaMax
  b_betaBand_exp <- BaseVariables$bBetaBand_1
  + BaseVariables$betaBandwidth * LambdaMax
  beta_band_exp <- BaseVariables$betaBandFactor * exp( -((Wavelength-lambda_beta_max_exp) / b_betaBand_exp)^2)
  tmp_a <- BaseVariables$a_1 + BaseVariables$a_2 * exp( - (LambdaMax-300)^2.0/11940.0)
  Activity <- 1/(exp(BaseVariables$A * (tmp_a - LambdaMax / Wavelength))
                 + exp(BaseVariables$B * (BaseVariables$b - LambdaMax / Wavelength))
                 + exp(BaseVariables$C * (BaseVariables$c - LambdaMax / Wavelength))
                 + BaseVariables$D)
  Activity <- Activity * (exp(BaseVariables$A * (tmp_a - 1))
                          + exp(BaseVariables$B * (BaseVariables$b - 1))
                          + exp(BaseVariables$C * (BaseVariables$c - 1))
                          + BaseVariables$D)
  Activity <- Activity + beta_band_exp
  Activity / (1+BaseVariables$betaBandFactor * exp( - ((LambdaMax-lambda_beta_max_exp) / (b_betaBand_exp))^2))
}

PeakOutput <- function(Lambda, x_i) {
  suppressWarnings(optimize(ActivityOptim,
                            interval=c(1, 800),
                            i=x_i,
                            LambdaMax=Lambda,
                            maximum=TRUE, tol = 1e-5)$maximum)
}

calculate_HDI <- function(x) {
  hdi95 <- unlist(bayestestR::hdi(x = x, ci = 0.95))
  hdi80 <- unlist(bayestestR::hdi(x = x, ci = 0.80))
  return(list(N=length(x),Mean=mean(x),Median=median(x),SD=sd(x),HDI95.CI_low = hdi95[2], HDI95.CI_high = hdi95[3], HDI80.CI_low = hdi80[2], HDI80.CI_high = hdi80[3]))
}

calculate_HDI_95 <- function(x) {
  hdi95 <- unlist(bayestestR::hdi(x = x, ci = 0.95))
  return(list(Mean=mean(x),Median=median(x),CIlow = hdi95[2], CIhigh = hdi95[3]))
}

lambda_max_Plate_Group <- function(Group, Plate, Sample=NA) {
  if(anyNA(Sample)) {
    Sample <- seq_along(BaseVariables$lambda_max)
  }
    if(Group == 1) {
      as.vector(exp(log(BaseVariables$lambda_max[Sample]) + betas_Plate_list[[3]][,Plate][Sample]))  
    } else {
      as.vector(exp(log(BaseVariables$lambda_max[Sample]) + beta_AlphaGroup[,Group-1][Sample] + betas_Plate_list[[3]][,Plate][Sample]))  
    }
}

activity_max_Plate_Group <- function(Group, Plate, Sample=NA) {
  if(anyNA(Sample)) {
    Sample <- seq_along(BaseVariables$activity_max)
  }
    if(Group == 1) {
        as.vector(exp(log(BaseVariables$activity_max[Sample]) + betas_Plate_list[[2]][,Plate][Sample]))
    } else {
        as.vector(exp(log(BaseVariables$activity_max[Sample]) + betas_Plate_list[[2]][,Plate][Sample] + beta_activity_max_Group[,Group-1][Sample]))
    }
}

shape_Plate_Group <- function(Group, Plate, Sample=NA) {
  if(anyNA(Sample)) {
    Sample <- seq_along(BaseVariables$shape)
  }
  as.vector(exp(log(BaseVariables$shape[Sample]) + betas_Plate_list[[1]][,Plate][Sample] + betas_Group[,Group][Sample]))
}

I_50_Group_Plate_Wavelength <- function(Group, Plate, Wavelength, Sample=NA) {
  if(anyNA(Sample)) {
    lambda_max_tmp <- lambda_max_Plate_Group(Group, Plate)
    Sample <- seq_along(lambda_max_tmp)
  } else {
    lambda_max_tmp <- lambda_max_Plate_Group(Group, Plate)[Sample]
  }
  spectrum <- ActivityOptim(LambdaMax=lambda_max_tmp, x=Wavelength, i=Sample)
  if(Group==1) {
    as.vector(exp(log(BaseVariables$I_50[Sample])) / spectrum)
  } else {
    as.vector(exp(log(BaseVariables$I_50[Sample]) + beta_I_50_Group[,Group-1][Sample]) / spectrum)
  }
}

dark_mu_Plate_Group <- function(Group, Plate, Sample=NA) {
  if(anyNA(Sample)) {
    Sample <- seq_along(BaseVariables$Dark_beta)
  }
    tmpOut <- BaseVariables$Dark_beta[Sample] + betas_Plate_list[[4]][,Plate][Sample]# + betas_Group_list[[3]][,Group]
    if(Group != 1) {
        return(exp(as.vector(tmpOut + Dark_Group_beta[,Group-1][Sample] + betas_Plate_list[[3+Group]][,Plate][Sample])))
        #        return(exp(as.vector(tmpOut + Dark_Group_beta[,Group-1] + betas_Plate_list[[4+Group]][,Plate])))
    } else {
        return(exp(as.vector(tmpOut)))
    }
}

PeakDetectionGroup <- function(Sample, GroupSelection, interval=c(200,800)) {
  tmp <- optimize(ActivityOptim,
                  interval=interval,
                  i=Sample,
                  LambdaMax=lambda_max_Group[Sample,GroupSelection],
                  maximum=TRUE, tol = 1e-5)
  return(list(PeakLocation=tmp$maximum, PeakAmplitude=tmp$objective))
  }

PeakDetection <- function(Sample, lambda_max, interval=c(400,700)) {
  tmp <- optimize(ActivityOptim,
                  interval=interval,
                  i=Sample,
                  LambdaMax=lambda_max,
                  maximum=TRUE, tol = 1e-5)
  return(list(PeakLocation=tmp$maximum, PeakAmplitude=tmp$objective))
  }

Spectrum <- function(x, Wavelength=381:661, Group, Plate, illumination_time, Sample=NA, HDI=F, dark_activity=T) {
  
  lambda_max_Plate_Group <- function(Group, Plate, Sample=NA) {
    if(anyNA(Sample)) {
      Sample <- seq_along(BaseVariables$lambda_max)
    }
    if(Group == 1) {
      as.vector(exp(log(BaseVariables$lambda_max[Sample]) + betas_Plate_list[[3]][,Plate][Sample]))  
    } else {
      as.vector(exp(log(BaseVariables$lambda_max[Sample]) + beta_AlphaGroup[,Group-1][Sample] + betas_Plate_list[[3]][,Plate][Sample]))  
    }
  }
  
  activity_max_Plate_Group <- function(Group, Plate, Sample=NA) {
    if(anyNA(Sample)) {
      Sample <- seq_along(BaseVariables$activity_max)
    }
    if(Group == 1) {
      as.vector(exp(log(BaseVariables$activity_max[Sample]) + betas_Plate_list[[2]][,Plate][Sample]))
    } else {
      as.vector(exp(log(BaseVariables$activity_max[Sample]) + betas_Plate_list[[2]][,Plate][Sample] + beta_activity_max_Group[,Group-1][Sample]))
    }
  }
  
  shape_Plate_Group <- function(Group, Plate, Sample=NA) {
    if(anyNA(Sample)) {
      Sample <- seq_along(BaseVariables$shape)
    }
    as.vector(exp(log(BaseVariables$shape[Sample]) + betas_Plate_list[[1]][,Plate][Sample] + betas_Group[,Group][Sample]))
  }
  
  I_50_Group_Plate_Wavelength <- function(Group, Plate, Wavelength, Sample=NA) {
    lambda_max_tmp <- lambda_max_Plate_Group(Group, Plate)[Sample]
    if(anyNA(Sample)) {
      Sample <- seq_along(lambda_max_tmp)
    }
    spectrum <- ActivityOptim(LambdaMax=lambda_max_tmp, x=Wavelength, i=Sample)
    if(Group==1) {
      as.vector(exp(log(BaseVariables$I_50[Sample])) / spectrum)
    } else {
      as.vector(exp(log(BaseVariables$I_50[Sample]) + beta_I_50_Group[,Group-1][Sample]) / spectrum)
    }
  }
  
  dark_mu_Plate_Group <- function(Group, Plate, Sample=NA) {
    if(anyNA(Sample)) {
      Sample <- seq_along(BaseVariables$Dark_beta)
    }
    tmpOut <- BaseVariables$Dark_beta[Sample] + betas_Plate_list[[4]][,Plate][Sample]# + betas_Group_list[[3]][,Group]
    if(Group != 1) {
      return(exp(as.vector(tmpOut + Dark_Group_beta[,Group-1][Sample] + betas_Plate_list[[3+Group]][,Plate][Sample])))
    } else {
      return(exp(as.vector(tmpOut)))
    }
  }
  
  ActivityOptim <- function(LambdaMax=500, x, i=1) {
    lambda_beta_max <-  BaseVariables$betaMax_1[i] + BaseVariables$betaMax_2[i] * LambdaMax
    
    b_betaBand <- BaseVariables$bBetaBand_1[i] + BaseVariables$betaBandwidth[i] * LambdaMax
    beta_band <- BaseVariables$betaBandFactor[i] * exp( -((x-lambda_beta_max) / b_betaBand)^2);
    max_div_lambda <- LambdaMax / x;
    
    mu_vec <- 1/(exp(BaseVariables$A[i] * (BaseVariables$a_1[i]+BaseVariables$a_2[i]*exp( - ((LambdaMax-300)^2)/11940.0) - max_div_lambda)) +
                   exp(BaseVariables$B[i] * (BaseVariables$b[i] - max_div_lambda)) +
                   exp(BaseVariables$C[i] * (BaseVariables$c[i] - max_div_lambda)) +
                   BaseVariables$D[i])
    mu_vec <- mu_vec * (exp(BaseVariables$A[i] * (BaseVariables$a_1[i]+BaseVariables$a_2[i]*exp(- ((LambdaMax-300)^2)/11940)-1)) +
                          exp(BaseVariables$B[i] * (BaseVariables$b[i]-1)) +
                          exp(BaseVariables$C[i] * (BaseVariables$c[i]-1)) +
                          BaseVariables$D[i])
    mu_vec <- mu_vec + beta_band
    mu_vec <- mu_vec / (1+BaseVariables$betaBandFactor[i] * exp( - ((LambdaMax-lambda_beta_max) / (b_betaBand))^2))
    return(mu_vec)
  }
  
  BaseVariables <- posterior::subset_draws(posterior::as_draws_df(x), variable = c("A", "B", "C", "D",
                                                                                   "a_1", "a_2", "b", "c",
                                                                                   "betaMax_1", "betaMax_2",
                                                                                   "betaBandFactor", "betaBandwidth", "bBetaBand_1",
                                                                                   "shape", "lambda_max", "Dark_beta", "activity_max", "I_50"))
  if(anyNA(Sample) | HDI) {
    if(HDI & (length(Sample)<length(BaseVariables$lambda_max))) {
      if(!anyNA(Sample)) {
        message("Overwrite Sample: Use full length of Sample for HDI")
      }
    }
    Sample <- seq_along(BaseVariables$lambda_max)
  }
  posterior::subset_draws(posterior::as_draws_matrix(x), variable = "beta_AlphaGroup")
  ModelDraws <- posterior::as_draws_matrix(x)
  beta_AlphaGroup <- posterior::subset_draws(ModelDraws, variable = c("beta_AlphaGroup"))
  beta_activity_max_Group <- posterior::subset_draws(ModelDraws, variable = c("beta_activity_max_Group"))
  Dark_Group_beta <-  posterior::subset_draws(ModelDraws, variable = c("Dark_Group_beta"))
  beta_I_50_Group <- posterior::subset_draws(ModelDraws, variable = c("beta_I_50_Group"))
  betas_Group <- posterior::subset_draws(ModelDraws, variable = c("betas_Group"))
  GroupN <- ncol(betas_Group)
  betas_Plate <- posterior::subset_draws(ModelDraws, variable = c("betas_Plate"))
  betas_Plate_list <- lapply(X = 1:(3+GroupN), FUN = function(plate) {
    betas_Plate[,grepl(pattern = paste0(",",plate,"\\]"), x = colnames(betas_Plate))]
  })
  
  GroupPlateIndex <- data.table(Group, Plate,illumination_time)
  LambdaMaxPlateGroup <- GroupPlateIndex[,.(LambdaMax=lambda_max_Plate_Group(Group, Plate, Sample),
                                            DarkMu=dark_mu_Plate_Group(Group, Plate,Sample)),
                                         by=.(Group,Plate,illumination_time)][,Sample:=Sample,by=.(Group,Plate,illumination_time)]
  if(HDI==T) {
      I50_activityMax_Group_Plate_DT <- GroupPlateIndex[,.(Wavelength=Wavelength),by=.(Group,Plate,illumination_time)][
        ,.(DarkActivity=dark_mu_Plate_Group(Group,Plate,Sample),
           activityMax=activity_max_Plate_Group(Group,Plate,Sample),
           I50=I_50_Group_Plate_Wavelength(Group, Plate, Wavelength = Wavelength,Sample),
           shape=shape_Plate_Group(Group, Plate,Sample),
           Sample=Sample),by=.(Group,Plate,Wavelength, illumination_time)]
      if(dark_activity) {
        I50_activityMax_Group_Plate_DT[,Activity_mu:=(activityMax*illumination_time)/(I50+illumination_time)+DarkActivity,][,Activity:=rgamma(n = .N, shape = shape, rate = shape/Activity_mu),]
      } else {
        I50_activityMax_Group_Plate_DT[,Activity_mu:=(activityMax*illumination_time)/(I50+illumination_time),][,Activity:=rgamma(n = .N, shape = shape, rate = shape/Activity_mu),]
      }
      rbindlist(l = list(I50_activityMax_Group_Plate_DT[,c(calculate_HDI(x = Activity_mu), "Type"="Activity_mu"),by=.(Group, Plate,Wavelength)],I50_activityMax_Group_Plate_DT[,c(calculate_HDI(x = Activity), "Type"="Activity"),by=.(Group, Plate,Wavelength)]))
  } else {
    I50_activityMax_Group_Plate_DT <- GroupPlateIndex[,.(Wavelength=Wavelength),by=.(Group,Plate,illumination_time)][
      ,.(DarkActivity=dark_mu_Plate_Group(Group,Plate,Sample),
         activityMax=activity_max_Plate_Group(Group,Plate,Sample),
         I50=I_50_Group_Plate_Wavelength(Group, Plate, Wavelength = Wavelength,Sample),
         shape=shape_Plate_Group(Group, Plate,Sample),
         Sample=Sample),by=.(Group,Plate,Wavelength, illumination_time)]
    if(dark_activity) {
      I50_activityMax_Group_Plate_DT[,Activity_mu:=(activityMax*illumination_time)/(I50+illumination_time)+DarkActivity,][,Activity:=rgamma(n = .N, shape = shape, rate = shape/Activity_mu),]
    } else {
      I50_activityMax_Group_Plate_DT[,Activity_mu:=(activityMax*illumination_time)/(I50+illumination_time),][,Activity:=rgamma(n = .N, shape = shape, rate = shape/Activity_mu),]
    }
    return(I50_activityMax_Group_Plate_DT)
  }
}


DiffPaste <- function(x, y) {
  paste(x, y, sep = " - ")
}

ModelRun <- function(stanData,
                     n_chains,
                     n_samples,
                     output_path,
                     adapt_delta = 0.99,
                     iter_warmup = 2000,
                     save = T,
                     prior_only = F) {
  Model <- cmdstanr::cmdstan_model("stan/SpectralEstimation.stan")
  if(prior_only) {
    stanDataTmp <- stanData
    stanDataTmp$priorOnly <- 1
  } else {
    stanDataTmp <- stanData
    stanDataTmp$priorOnly <- 0
  }

  ModelFit <- Model$sample(data = stanDataTmp,
                           init = InitChains(Chains = n_chains, stanData = stanDataTmp),
                           parallel_chains = n_chains,
                           chains = n_chains,
                           iter_warmup = iter_warmup,
                           iter_sampling = n_samples,
                           seed=123,
                           refresh = 200,
                           adapt_delta = adapt_delta,
                           save_latent_dynamics = F)
  if(save) {
    if(prior_only) {
      dirPath <- paste0(output_path, "Modelout_prior.rds")
    } else {
      dirPath <- paste0(output_path, "Modelout.rds")
    }
    ModelFit$save_object(file = dirPath)
  }
  return(ModelFit)
}
