setwd(dir = projectDir)
source("R/LoadDataFunctions.R")

plan(multisession, workers = n_cores)
AxisTickLength <- 10

ModelDirectory <- "Output/20240611_10sATR/"

Model <- readRDS(file = paste0(ModelDirectory,"Modelout.rds"))
ConstructColours <- c("bovRho-mScarlet.ATR"="green","JellyOp-mScarlet.ATR"="orange")
dir.create(path = paste0(ModelDirectory,"Plots"), showWarnings = F)
dir.create(path = paste0(ModelDirectory,"Tables"), showWarnings = F)

AssayList <- LoadData("fullATR10s", illumination_range = c(10000, 10000))
AssayData <- AssayList$AssayData
AssayData <- merge(AssayData,AssayList$AssayDataControl[,.(DarkActivity=mean(Dark_Activity)),by=ExpID], by="ExpID", allow.cartesian=T, all.x=T)
AssayData[,Activity_corrected:=Activity-DarkActivity,]
stanData <- ConstructStanData(AssayList)

WavelengthVec <- AssayData[,sort(unique(Wavelength)),]

GroupNames <- AssayData[,levels(ConstructRetinal),]
GroupN <- length(GroupNames)
BaseVariables <- Model$draws(variables = c("A", "B", "C", "D",
                                           "a_1", "a_2", "b", "c",
                                           "betaMax_1", "betaMax_2",
                                           "betaBandFactor", "betaBandwidth", "bBetaBand_1",
                                           "shape", "lambda_max", "Dark_beta", "activity_max", "I_50"), format = "df")
beta_AlphaGroup <- Model$draws(variables = c("beta_AlphaGroup"), format = "matrix")
beta_activity_max_Group <- Model$draws(variables = c("beta_activity_max_Group"), format = "matrix")
Dark_Group_beta <- Model$draws(variables = c("Dark_Group_beta"), format = "matrix")
beta_I_50_Group <- Model$draws(variables = c("beta_I_50_Group"), format = "matrix")

betas_Plate <- Model$draws(variables = c("betas_Plate"), format = "matrix")
betas_Plate_list <- lapply(X = 1:(3+GroupN), FUN = function(x) {
  betas_Plate[,grepl(pattern = paste0(",",x,"\\]"), x = colnames(betas_Plate))]
})

betas_Group <- Model$draws(variables = c("betas_Group"), format = "matrix")
betas_Group_list <- lapply(X = 1, FUN = function(x) {
  betas_Group[,grepl(pattern = paste0(",",x,"\\]"), x = colnames(betas_Group))]
})

GroupPlateIndex <- AssayData[,.N,by=.(Group=as.integer(ConstructRetinal), Plate)]


### compute scaling for different groups and illumination times
IlluminationTimes <- AssayData[,sort(unique(illumination_time)),]
AssayData[,MeasurementIndex:=.I]

### Peak detection and quantification
LambdaMaxPlateGroup <- GroupPlateIndex[,.(LambdaMax=lambda_max_Plate_Group(Group, Plate)),by=.(Group,Plate)][,Sample:=seq_len(.N),by=.(Group,Plate)]
LambdaMaxPlateGroup <- rbindlist(foreach(x = split(x = LambdaMaxPlateGroup, by=c("Group", "Plate"))) %dofuture% {
    tmp <- copy(x)
    tmp[,PeakDetection(Sample = Sample, lambda_max = LambdaMax),by=.(Sample, Group, Plate)]
})
for(GroupSelection in 1:GroupN) {
  LambdaMaxPlateGroup[Group==GroupSelection,GroupName:=GroupNames[GroupSelection],]
}
LambdaMaxPlateGroup[,GroupName:=factor(x=GroupName, levels = c("bovRho-mScarlet.ATR","JellyOp-mScarlet.ATR")),]

CompMat <- outer(levels(LambdaMaxPlateGroup[,GroupName,]), levels(LambdaMaxPlateGroup[,GroupName,]), "DiffPaste")
CompVec <- CompMat[lower.tri(CompMat)]
LambdaMaxPlateGroupDiff <- rbindlist(lapply(CompVec, function(x) {
  GroupComparison <- unlist(strsplit(x, " - "))
  tmpDiff <- merge(LambdaMaxPlateGroup[GroupName == GroupComparison[1],.(PeakLocation1=mean(PeakLocation)),by=.(GroupName, Sample)], LambdaMaxPlateGroup[GroupName == GroupComparison[2],.(PeakLocation2=mean(PeakLocation)),by=.(GroupName, Sample)], by=c("Sample"))
  return(tmpDiff[,.(PeakDiff=PeakLocation1-PeakLocation2),by=Sample][,GroupComparison:=x,])
}))


LambdaMaxPlateGroup[,.(PeakLocation=mean(PeakLocation)),by=.(Sample, GroupName)][,calculate_HDI(PeakLocation),by=GroupName]
GroupComparisonProb <- LambdaMaxPlateGroupDiff[,.(prob=min(sum(PeakDiff>0)/.N,sum(PeakDiff<0)/.N),
                           direction=fifelse(sum(PeakDiff>0)/.N<0.5, yes="smaller", no="greater")),by=GroupComparison]
GroupComparisonHDI <- LambdaMaxPlateGroupDiff[,calculate_HDI(PeakDiff),by=GroupComparison]

(PeakDifferencePlot <- ggplot(LambdaMaxPlateGroupDiff, aes(x=PeakDiff, group=GroupComparison))+#, group=Group, colour=as.factor(Group)))+
  stat_halfeye(aes(fill = after_stat(x > 0)), .width=c(0.8,0.9))+
  scale_fill_manual(values=c("grey", "lightblue"))+
  coord_cartesian(expand = F, ylim = c(0,1), clip="off")+
  scale_y_continuous(name = "Density")+
  geom_vline(xintercept = 0, linetype="dashed")+
  facet_wrap(~GroupComparison)+
  theme_classic()+
  theme(strip.background = element_blank(), strip.text = element_text(size=12), axis.ticks.length = unit(x=AxisTickLength, units="pt"))+
  guides(x = guide_axis(minor.ticks = TRUE)))
ggsave(paste0(ModelDirectory,"Plots/PeakDifferencePlot.pdf"), PeakDifferencePlot, width = 10, height = 10)
saveRDS(PeakDifferencePlot, file = paste0(ModelDirectory,"Plots/PeakDifferencePlot.rds"))

(ConstructPeaksPlatePlot <- ggplot(LambdaMaxPlateGroup, aes(x=PeakLocation, y=Plate, group=interaction(GroupName,Plate), fill=GroupName, colour=GroupName))+
  stat_slab(height=1.5, colour="black", alpha=0.6)+
  geom_vline(data = LambdaMaxPlateGroup[,.(Mean=mean(PeakLocation)),by=GroupName], linetype="dashed", aes(xintercept = Mean, colour=GroupName))+
  coord_cartesian(expand = F)+
  scale_y_continuous(name = "Plate",  breaks = GroupPlateIndex$Plate)+
  scale_x_continuous(name = "wavelength in nm",  breaks = seq(0,800,10))+
  scale_colour_manual(name = "Construct", values=ConstructColours, aesthetics=c("colour", "fill"))+
  theme_classic()+
  theme(strip.background = element_blank(), strip.text = element_text(size=12), axis.ticks.length = unit(x=AxisTickLength, units="pt"))+
  guides(x = guide_axis(minor.ticks = TRUE)))
ggsave(paste0(ModelDirectory,"Plots/ConstructPeaksPlatePlot.pdf"), ConstructPeaksPlatePlot, width = 10, height = 10)
saveRDS(ConstructPeaksPlatePlot, file = paste0(ModelDirectory,"Plots/ConstructPeaksPlatePlot.rds"))

(ConstructPeaksPlot <- ggplot(LambdaMaxPlateGroup[,.(PeakLocation=mean(PeakLocation)),by=.(GroupName, Sample)], aes(x=PeakLocation, group=GroupName, fill=GroupName, colour=GroupName))+
  stat_slab(height=1.5, colour="black", alpha=0.6)+
  geom_vline(data = LambdaMaxPlateGroup[,.(PeakLocation=mean(PeakLocation)),by=.(GroupName, Sample)][,.(Mean=mean(PeakLocation)),by=GroupName], linetype="dashed", aes(xintercept = Mean, colour=GroupName))+
  coord_cartesian(expand = F)+
  scale_y_continuous(name = "Density")+
  scale_colour_manual(name = "Construct", values=ConstructColours, aesthetics=c("colour", "fill"))+
  scale_x_continuous(name = "wavelength in nm",  breaks = seq(0,800,10))+
  theme_classic()+
  theme(strip.background = element_blank(), strip.text = element_text(size=12), axis.ticks.length = unit(x=AxisTickLength, units="pt"))+
  guides(x = guide_axis(minor.ticks = TRUE)))
ggsave(paste0(ModelDirectory,"Plots/ConstructPeaksPlot.pdf"), ConstructPeaksPlot, width = 10, height = 10)
saveRDS(ConstructPeaksPlot, file = paste0(ModelDirectory,"Plots/ConstructPeaksPlot.rds"))

KableStringConstructPeaks <- gsub(pattern = "NA", replacement = "  ", x = knitr::kable(LambdaMaxPlateGroup[,.(PeakLocation=mean(PeakLocation)),by=.(GroupName, Sample)][,calculate_HDI(PeakLocation),by=GroupName][,-2,], digits = c(0,1,1,1,1,1,1,1)))
KableStringConstructPeakDiff <- gsub(pattern = "NA", replacement = "  ", x = knitr::kable(GroupComparisonHDI[,-2,], digits = c(0,1,1,1,1,1,1,1)))
KableStringConstructPeakDiffProb <- gsub(pattern = "NA", replacement = "  ", x = knitr::kable(GroupComparisonProb, digits = c(0,6)))

RmdString <- paste0(ModelDirectory,"Tables/Summary.Rmd")
cat(paste0("## Construct Comparison\n"),
    "\n### Activity Peaks (in nm)",
    KableStringConstructPeaks,
    "\n### Differences in Activity Peaks (in nm)", 
    KableStringConstructPeakDiff,
    "\n### Difference Probabilities (from 0 nm shift)", 
    KableStringConstructPeakDiffProb,
    sep="\n", file=RmdString)
render(RmdString, output_dir = paste0(ModelDirectory,"Tables/"), output_format = c("pdf_document"), clean = T)
file.remove(RmdString)

Spectrum_plate_DT <- rbindlist(foreach(x = 350:700) %dofuture% {
    GroupPlateIndex[,.(Wavelength = x, Spectrum=list(ActivityOptim(LambdaMax = lambda_max_Plate_Group(Group, Plate), x = x, i = seq_along(BaseVariables$lambda_max)))),by=.(Group,Plate)]
})


Spectrum_plate_DT <- Spectrum_plate_DT[,.(Spectrum=unlist(Spectrum), Sample = 1:length(unlist(Spectrum))),by=.(Group, Plate, Wavelength)]
Spectrum_plate_DT <- merge(Spectrum_plate_DT, LambdaMaxPlateGroup, by=c("Group", "Sample", "Plate"))

Spectrum_plate_DT[,AdjustedSpectrum:=Spectrum/PeakAmplitude,]
Spectrum_Group_HDI <- Spectrum_plate_DT[,.(AdjustedSpectrum=mean(AdjustedSpectrum)),by=.(Group, Wavelength,GroupName, Sample)][,calculate_HDI(AdjustedSpectrum),by=.(Group, Wavelength,GroupName)]


(NormalisedGroupSpectrumLog <- ggplot(data=Spectrum_Group_HDI[Wavelength %between% c(381,661)], aes(y=Mean, x=Wavelength, group=Group, colour=GroupName))+
  geom_ribbon(aes(ymin=HDI95.CI_low, ymax=HDI95.CI_high, fill=GroupName,group=Group,x=Wavelength),alpha=0.1, colour=NaN)+
  geom_line()+
  coord_cartesian(ylim=c(1e-5,1), expand = F)+
  scale_colour_manual(name = "Construct", values=ConstructColours, aesthetics=c("colour", "fill"))+
  scale_y_continuous(trans="log10", name="relative photosensitivity", breaks=c(10^(-5:0)), minor_breaks = unique(c((1:10)/1e5,(1:10)/1e4,(1:10)/1e3,(1:10)/1e2,(1:10)/1e1,1:10)), labels=expression(10^{-5}, 10^{-4}, 10^{-3}, 10^{-2}, 10^{-1}, 10^{0}),expand = c(0,0))+
  scale_x_continuous(name = "wavelength in nm",  breaks = seq(300,700,100), minor_breaks = seq(300,700,20),limits=c(360,700), expand = c(0,0))+
  theme_classic()+
  theme(axis.ticks.length = unit(x=AxisTickLength, units="pt"),legend.position = "inside", legend.position.inside =  c(.4, .1))+
  guides(y = guide_axis(minor.ticks = TRUE), x = guide_axis(minor.ticks = TRUE)))
ggsave(paste0(ModelDirectory,"Plots/NormalisedGroupSpectrumLog.pdf"), NormalisedGroupSpectrumLog, width = 10, height = 10)
saveRDS(NormalisedGroupSpectrumLog, file = paste0(ModelDirectory,"Plots/NormalisedGroupSpectrumLog.rds"))

(NormalisedGroupSpectrum <- ggplot(data=Spectrum_Group_HDI[Wavelength %between% c(381,661)], aes(y=Mean, x=Wavelength, group=Group, colour=GroupName))+
  geom_ribbon(aes(ymin=HDI95.CI_low, ymax=HDI95.CI_high, fill=GroupName,group=Group,x=Wavelength),alpha=0.1, colour=NaN)+
  geom_line()+
  coord_cartesian(expand = F, clip="off")+
  scale_colour_manual(name = "Construct", values=ConstructColours, aesthetics=c("colour", "fill"))+
  scale_y_continuous(name="relative photosensitivity", limits=c(0,1), expand=c(0,0))+
  scale_x_continuous(name = "wavelength in nm",  breaks = seq(300,700,100), minor_breaks = seq(300,700,20),limit=c(360,700, expand=c(0,0)))+
  theme_classic()+
  theme(axis.ticks.length = unit(x=AxisTickLength, units="pt"),legend.position = "inside", legend.position.inside =  c(.4, .1))+
  guides(x = guide_axis(minor.ticks = TRUE)))
ggsave(paste0(ModelDirectory,"Plots/NormalisedGroupSpectrum.pdf"), NormalisedGroupSpectrum, width = 10, height = 10)
saveRDS(NormalisedGroupSpectrum, file = paste0(ModelDirectory,"Plots/NormalisedGroupSpectrum.rds"))

### calculate actual activity from I50 and activityMax
I50_activityMax_Group_Plate_DT <- GroupPlateIndex[,.(Wavelength=seq(380,700,5)),by=.(Group,Plate)][,.(DarkActivity=dark_mu_Plate_Group(Group,Plate), activityMax=activity_max_Plate_Group(Group,Plate), I50=I_50_Group_Plate_Wavelength(Group, Plate, Wavelength = Wavelength), shape=shape_Plate_Group(Group, Plate)),by=.(Group,Plate,Wavelength)]
I50_activityMax_Group_Plate_DT[,Activity_mu:=(activityMax*10000)/(I50+10000)+DarkActivity,]
I50_activityMax_Group_Plate_DT[,Activity_corrected:=(activityMax)*10000/(I50+10000),]
for(GroupSelection in 1:GroupN) {
  I50_activityMax_Group_Plate_DT[Group==GroupSelection,GroupName:=GroupNames[GroupSelection],]
}
I50_activityMax_Group_Plate_DT[,Activity:=rgamma(.N, shape, shape/Activity_mu),]
I50_activityMax_Group_Plate_DT[,Sample:=1:.N,by=.(Group, Plate, Wavelength)]
Activity_Group_HDI <- I50_activityMax_Group_Plate_DT[,.(Activity_mu=mean(Activity_mu)),by=.(Sample,Group, Wavelength,GroupName)][,calculate_HDI(Activity_mu),by=.(Group, Wavelength,GroupName)]
Activity_corrected_Group_HDI <- I50_activityMax_Group_Plate_DT[,.(Activity_corrected=mean(Activity_corrected)),by=.(Sample,Group, Wavelength,GroupName)][,calculate_HDI(Activity_corrected),by=.(Group, Wavelength,GroupName)]
Activity_Group_HDI[,ConstructRetinal:=GroupName,][,ConstructRetinal:=factor(ConstructRetinal, levels = c("JellyOp-mScarlet.ATR", "bovRho-mScarlet.ATR"))]
Activity_corrected_Group_HDI[,ConstructRetinal:=GroupName,][,ConstructRetinal:=factor(ConstructRetinal, levels = c("JellyOp-mScarlet.ATR", "bovRho-mScarlet.ATR"))]


(GroupActivityWithDarkSpectrum <- ggplot(data=Activity_Group_HDI[Wavelength %between%c(380,665),], aes(y=Mean, x=Wavelength, group=Group, colour=GroupName))+
  geom_ribbon(aes(ymin=HDI95.CI_low, ymax=HDI95.CI_high, fill=GroupName,group=Group,x=Wavelength),alpha=0.1, colour=NaN)+
  geom_line()+
   scale_colour_manual(name = "Construct", values=ConstructColours, aesthetics=c("colour", "fill"))+
  geom_point(data=AssayData[,.(Activity=mean(Activity)),by=.(Wavelength, ConstructRetinal)], aes(x=Wavelength, y=Activity, group=ConstructRetinal, colour=ConstructRetinal), size=3, inherit.aes = F)+
  geom_errorbar(data=AssayData[,.(Mean=mean(Activity), SD=sd(Activity), SEM=sd(Activity)/sqrt(.N)),by=.(Wavelength, ConstructRetinal)][,ConstructRetinal:=factor(ConstructRetinal, levels=c("JellyOp-mScarlet.ATR", "bovRho-mScarlet.ATR")),], aes(x=Wavelength, ymin=Mean-SEM, ymax=Mean+SEM, group=ConstructRetinal, colour=ConstructRetinal), width=6, inherit.aes = F)+
  scale_y_continuous(name="luminescence (norm.)", expand=c(0,0), breaks=seq(0,0.3,0.05), limits=c(0,0.25))+
  scale_x_continuous(name = "wavelength in nm",  breaks = seq(300,700,100), minor_breaks = seq(300,700,20), expand=c(0,0), limit=c(360,700))+
  theme_classic()+theme(strip.background = element_blank(), strip.text = element_text(size=12),axis.ticks.length = unit(x=AxisTickLength, units="pt"),legend.position = "inside", legend.position.inside =  c(.95, .9))+
    guides(x = guide_axis(minor.ticks = TRUE)))
ggsave(paste0(ModelDirectory,"Plots/GroupActivityWithDarkSpectrum.pdf"), GroupActivityWithDarkSpectrum, width = 10, height = 10)
saveRDS(GroupActivityWithDarkSpectrum, file = paste0(ModelDirectory,"Plots/GroupActivityWithDarkSpectrum.rds"))

(GroupActivityCorrectedSpectrum <- ggplot(data=Activity_corrected_Group_HDI[Wavelength %between%c(380,665),], aes(y=Mean, x=Wavelength, group=Group, colour=GroupName))+
  geom_ribbon(aes(ymin=HDI95.CI_low, ymax=HDI95.CI_high, fill=GroupName,group=Group,x=Wavelength),alpha=0.1, colour=NaN)+
  geom_line()+
   scale_colour_manual(name = "Construct", values=ConstructColours, aesthetics=c("colour", "fill"))+
  geom_point(data=AssayData[,.(Activity=mean(Activity_corrected, na.rm=T)),by=.(Wavelength, ConstructRetinal)], aes(x=Wavelength, y=Activity, group=ConstructRetinal, colour=ConstructRetinal), size=3, inherit.aes = F)+
  geom_errorbar(data=AssayData[,.(Mean=mean(Activity_corrected, na.rm=T), SD=sd(Activity_corrected, na.rm=T), SEM=sd(Activity_corrected, na.rm=T)/sqrt(.N)),by=.(Wavelength, ConstructRetinal)][,ConstructRetinal:=factor(ConstructRetinal, levels=c("JellyOp-mScarlet.ATR", "bovRho-mScarlet.ATR")),], aes(x=Wavelength, ymin=Mean-SEM, ymax=Mean+SEM, group=ConstructRetinal, colour=ConstructRetinal), width=6, inherit.aes = F)+
  scale_y_continuous(name="luminescence (norm.)", expand=c(0,0), breaks=seq(0,0.3,0.05), limits=c(-0.02,0.25))+
  geom_hline(yintercept = 0, linetype="dashed")+
  scale_x_continuous(name = "wavelength in nm",  breaks = seq(300,700,100), minor_breaks = seq(300,700,20), expand=c(0,0), limit=c(360,700))+
  theme_classic()+theme(strip.background = element_blank(), strip.text = element_text(size=12),axis.ticks.length = unit(x=AxisTickLength, units="pt"),legend.position = "inside", legend.position.inside =  c(.95, .9))+
  guides(x = guide_axis(minor.ticks = TRUE)))
ggsave(paste0(ModelDirectory,"Plots/GroupActivityCorrectedSpectrum.pdf"), GroupActivityCorrectedSpectrum, width = 10, height = 10)
saveRDS(GroupActivityCorrectedSpectrum, file = paste0(ModelDirectory,"Plots/GroupActivityCorrectedSpectrum.rds"))