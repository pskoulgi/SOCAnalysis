library(dplyr)
library(tidyr)
library(reshape2)
library(broom)
library(ggplot2)

# Read-in soil data
{
  soilCN <- read.table("F:/Workspace/MSc/Project/soc-fragments/Pradeep_CoorgSoilCN.csv",
                       as.is=T, header=T, sep=",")
  
  # Clean-up col names, drop needless cols
  soilCN <- dplyr::select(soilCN, c(SAMPLE, MASS, C_PERC, N_PERC))
  soilCN <- soilCN %>% rename(POINT=SAMPLE) %>% rename(SOIL.C=C_PERC) %>%
    rename(SOIL.N=N_PERC) %>% rename(SOIL.WT=MASS)
  
  # Remove "Blank", "Soil 502-309" standard sample rows
  soilCN <- filter(soilCN, grepl("CT|FR", POINT))
  
  # Copy CT4_PT0 row to CT5_PT0 row since 
  # they are the same point
  t <- soilCN[grep("CT4_PT0",soilCN$POINT), ]
  t$POINT[1] <- "CT5_PT0"
  soilCN <- bind_rows(soilCN, t)
  
  soilCN <- separate(soilCN, col="POINT", into=c("SITE.ID", "POINT.TYPE", "DUMMY", "POINT.IN.SITE", "POINT.NO"), sep=c(3,4,6,7), remove=F)
  soilCN <- dplyr::select(soilCN, -DUMMY)
  soilCN$SITE.ID <- as.factor(soilCN$SITE.ID)
  soilCN <- soilCN %>% group_by(SITE.ID) %>%
    mutate(ADJ.SOIL.C=SOIL.C-SOIL.C[POINT.IN.SITE=="0"])

  # Group points according to sites
  soilCN <- arrange(soilCN, POINT)
  
  # Read-in soil SIR data
  soilSIR <- read.table("F:/Workspace/MSc/Project/soc-fragments/soilSIR.csv",
                        as.is=T, header=T, sep=",")
  
  concNaOH <- 2 # in N
  
  avgBlanks <- soilSIR %>% group_by(SET) %>% filter(grepl("B", SAMPLE)) %>% 
    summarise(AVG.BLANKS=mean(TITRE))
  
  soilSIRS1 <- soilSIR %>% filter(SET == "S1") %>% 
    mutate(ADJ.TITRE = -TITRE+avgBlanks$AVG.BLANKS[avgBlanks$SET=="S1"])
  soilSIRS2 <- soilSIR %>% filter(SET == "S2") %>% 
    mutate(ADJ.TITRE = -TITRE+avgBlanks$AVG.BLANKS[avgBlanks$SET=="S2"])
  soilSIR <- bind_rows(soilSIRS1, soilSIRS2)
  
  soilSIR$ADJ.TITRE[soilSIR$ADJ.TITRE<0] = NA  
  soilSIR <- filter(soilSIR, !grepl("B", SAMPLE))
  
  soilSIR <- mutate(soilSIR, C.RELEASED = 1000000*
                      (ADJ.TITRE*(HCL_N/concNaOH))*12/(1000*26*10))
  soilSIR$C.RELEASED[soilSIR$C.RELEASED>25] = NA
  soilSIR <- soilSIR %>% rename(POINT=SAMPLE) %>%
    dplyr::select(POINT, C.RELEASED)
  
  belowGround <- left_join(soilCN, soilSIR, by="POINT")
  belowGround <- rename(belowGround, SOIL.SIR=C.RELEASED)
}

# Read-in aboveground data
{
  # Litter C, N data
  litterCN <- read.table("F:/Workspace/MSc/Project/soc-fragments/Pradeep_CoorgLitterCN.csv",
                         as.is=T, header=T, sep=",")
  
  litterCN <- dplyr::select(litterCN, c(SAMPLE, C_PERC, N_PERC))
  litterCN <- litterCN %>% rename(POINT=SAMPLE) %>% rename(LITTER.C=C_PERC) %>%
    rename(LITTER.N=N_PERC)
  
  # Remove "Blank", "EDTA*" standard sample rows
  litterCN <- filter(litterCN, grepl("CT|FR", POINT))
  
  # Group points according to sites
  litterCN <- arrange(litterCN, POINT)
  
  # Calculate C/N ratio
  litterCN <- mutate(litterCN, LITTER.CN.RATIO=LITTER.C/LITTER.N)
  
  # Litter P data
  litterP <- read.table("F:/Workspace/MSc/Project/soc-fragments/ICP_calculated_P_pradeep.csv",
                        as.is=T, header=T, sep=",")
  
  litterP <- dplyr::select(litterP, c(SAMPLE, P_PERC))
  litterP <- litterP %>% rename(POINT=SAMPLE) %>% rename(LITTER.P=P_PERC) %>%
    filter(grepl("CT|FR", POINT))
  
  # Litter weights
  litterWts <- read.table("F:/Workspace/MSc/Project/soc-fragments/litterProperties.csv",
                          as.is=T, header=T, sep=",")
  litterWts$AVGWT <- rowMeans(litterWts[ ,c("DRYWT1", "DRYWT2", "DRYWT3")], na.rm=TRUE)
  litterAvgWt <- litterWts %>% dplyr::select(c(POINT, AVGWT)) %>% rename(LITTER.WT = AVGWT)
  
  litter <- left_join(litterCN, litterP, by="POINT")
  litter <- left_join(litter, litterAvgWt, by="POINT")
  
  # Tree influence potential
  treeDistrib <- read.table("F:/Workspace/MSc/Project/soc-fragments/treeInfluencePotentialData.csv",
                            as.is=T, header=T, sep=",")
  
  gbhs <- dplyr::select(treeDistrib, c(GBH1, GBH2, GBH3, GBH4))
  
  treeDistrib$DBH <- (gbhs/pi)^2 %>% rowSums(na.rm=T) %>% sqrt()
  
  treeDistrib <- treeDistrib %>% separate(col="POINT",
                                          into=c("SITE.TYPE", "T"), sep=2,
                                          remove=F) %>%
    dplyr::select(-T) %>% separate(col="POINT",
                                   into=c("SITE.ID", "PT.TYPE", "T", "PT.IN.SITE", "TT"),
                                   sep=c(3,4,6,7), remove=F) %>% dplyr::select(-T, -TT)
  
  treeDistrib <- treeDistrib %>% filter(DBH>=0.1) %>% filter(PT.TYPE!="D") %>%
    filter(PT.IN.SITE==1) %>% dplyr::select(POINT, SITE.TYPE, SITE.ID, DISTANCE, DBH)
    
  
  binBounds = c(0.1, 0.15, 0.25, 0.45, 0.85, 2)
  treeDistrib$SIZE.CLASS[(binBounds[1]<treeDistrib$DBH) &
                           (treeDistrib$DBH<=binBounds[2])] = "S1"
  treeDistrib$SIZE.CLASS[(binBounds[2]<treeDistrib$DBH) &
                           (treeDistrib$DBH<=binBounds[3])] = "S2"
  treeDistrib$SIZE.CLASS[(binBounds[3]<treeDistrib$DBH) &
                           (treeDistrib$DBH<=binBounds[4])] = "S3"
  treeDistrib$SIZE.CLASS[(binBounds[4]<treeDistrib$DBH) &
                           (treeDistrib$DBH<=binBounds[5])] = "S4"
  treeDistrib$SIZE.CLASS[(binBounds[5]<treeDistrib$DBH) &
                           (treeDistrib$DBH<=binBounds[6])] = "S5"
  treeDistrib$SIZE.CLASS <- as.ordered((treeDistrib$SIZE.CLASS))

  # Compute TIP
  cOpt = -10^2/log(1/(0.1))
  treeDistrib$t <- (treeDistrib$DBH)*exp(-(treeDistrib$DISTANCE)^2/cOpt)
  treeDistrib$POINT <- as.factor(treeDistrib$POINT)
  tip <- treeDistrib %>% group_by(POINT) %>% summarise(TIP=sum(t, na.rm=T))
  
  litter <- left_join(litter, tip, by="POINT")
  
  # Compute basal area
  treeDistrib$TreeBasArea <- (pi*(treeDistrib$DBH)^2)/4
  basArea <- treeDistrib %>% group_by(POINT) %>% 
    summarise(BAS.AREA=sum(TreeBasArea, na.rm=T)*10000/(pi*10^2))
  treeDens <- treeDistrib %>% group_by(POINT) %>% 
    summarise(TREE.DENS=n()*10000/(pi*10^2))
  largeTreeDens <- treeDistrib %>% mutate(IS.LARGE = DBH>0.7) %>%
    group_by(POINT) %>% summarise(NO.LARGE.TREES = sum(IS.LARGE)*10000/(pi*10^2), 
                                  PROP.LARGE.TREES = sum(IS.LARGE)/n())
  largeBasalAreaFrac <- treeDistrib %>% mutate(IS.LARGE = DBH>0.7) %>%
    group_by(POINT) %>% summarise(LARGE.TREE.BA = sum(TreeBasArea[DBH>0.7],
                                                      na.rm=T)*
                                                  10000/(pi*10^2),
                                  LARGE.TREE.BA.FRAC =sum(TreeBasArea[DBH>0.7])/
                                                      sum(TreeBasArea))
  largeTreeDens <- left_join(largeTreeDens, largeBasalAreaFrac, by="POINT")
  
  # Compute tree density by size class
  treeDensBySizeClass <- treeDistrib %>% group_by(POINT, SIZE.CLASS) %>% 
    summarise(TREE.DENS=n()*10000/(pi*10^2))
  treeDensBySizeClass <- treeDensBySizeClass %>% group_by(POINT) %>% 
    mutate(PROP.STEMS=TREE.DENS/sum(TREE.DENS))
  treeDensBySizeClass <- treeDensBySizeClass %>% separate(col="POINT",
                           into=c("SITE.TYPE", "T"), sep=2,
                           remove=F) %>% dplyr::select(-T) %>% 
    separate(col="POINT", into=c("SITE.ID", "T"), sep=3, remove=F) %>% 
    dplyr::select(-T)
  
  treeDensMeans <- treeDensBySizeClass %>% group_by(SITE.TYPE, SIZE.CLASS) %>%
    summarise(MEAN=mean(TREE.DENS), SD=sd(TREE.DENS), SE=sd(TREE.DENS)/sqrt(n()))
  treeStemPropsMeans <- treeDensBySizeClass %>% group_by(SITE.TYPE, SIZE.CLASS) %>%
    summarise(MEAN=mean(PROP.STEMS), SD=sd(PROP.STEMS), SE=sd(PROP.STEMS)/sqrt(n()))
    
  aboveGround <- left_join(litter, basArea, by="POINT")
  aboveGround <- left_join(aboveGround, treeDens, by="POINT")
    
}

# Anand's data
{
  anandLitter <- read.table("F:/Workspace/MSc/Project/scratchPad/Anand_leaf_litter_elements.csv",
                         as.is=T, header=T, sep=",")
  
  anandLitter <- anandLitter %>% dplyr::select(Name:Forest.type, contains("Soil")) %>%
    filter(Forest.type != "") %>% mutate(Litter.CtoN=X..C/X.N)
  
  anandLitter$Forest.type[anandLitter$Forest.type != "Fragment"] = "Contiguous forest"

  
  
  anandSoil <- read.table("F:/Workspace/MSc/Project/scratchPad/Anand_soildata.csv",
                          as.is=T, header=T, sep=",")
}

# Anand's species-wise leaf litter quality data
{
  # Replicated species-wise litter quality data
  anandSpWiseLitter <- read.table("F:/Workspace/MSc/Project/scratchPad/Anand_specieswise_litter_CN.csv",
                            as.is=T, header=T, sep=",")
  anandSpWiseLitter <- anandSpWiseLitter %>% mutate (SP.CTON = C../N.) 
  # Average litter quality for each species
  speciesWiseLitterQual <- anandSpWiseLitter %>% 
    group_by(species) %>% summarise(CN.RATIO=mean(SP.CTON))
  
  # Plot-level tree-id and basal area measurements. 
  anandSpDominance <- read.table("F:/Workspace/MSc/Project/scratchPad/Anand_TreeCommunity_SpeciesDominance_plotData.csv",
                                  as.is=T, header=T, sep=",")
  anandSpDominance <- unite(anandSpDominance, SITE.PLOT.ID, site.name, plot.id,
                            remove=F)
  # Fractional contribution of each tree to its plot-level basal area
  anandSpDominance <- anandSpDominance %>% group_by(SITE.PLOT.ID) %>% 
    arrange(desc(basal.area_sq.cm)) %>% 
    mutate(BAS.AREA.FRAC = basal.area_sq.cm/sum(basal.area_sq.cm)) %>%
    mutate(BAS.AREA.CUMFRAC = cumsum(BAS.AREA.FRAC))
  # Biggest trees that contribute to top 70% basal area at plot level
  anand70PercDom <- anandSpDominance %>%
    dplyr::select(-index, -date, -family, -diameter_cm, -height_m, -Remarks) %>%
    filter(BAS.AREA.CUMFRAC <= 0.7)
  
  # LUT for tree id and code used
  anandSpeciesCodes <- read.table("F:/Workspace/MSc/Project/scratchPad/Anand_TreeCommunity_Families.csv",
                                  as.is=T, header=T, sep=",")
  anandSpeciesCodes <- anandSpeciesCodes %>% dplyr::select(code, combine) %>%
    rename(species=combine)
  
  # Retain 
  anand70PercDom <- inner_join(anandSpeciesCodes, anand70PercDom, by="species")
  anand70PercDom <- right_join(speciesWiseLitterQual, anand70PercDom,
                               by= c("species"="code"))
  
  anand70PercDomMySites <- filter (anand70PercDom, 
                                   site.name=="Kokka" |
                                     site.name=="Ruduraguppae" |
                                     site.name=="Arji" |
                                     site.name=="Arpattu")
  commAvgdLittQualMySites <- anand70PercDomMySites %>%
    group_by(SITE.PLOT.ID) %>%
    summarise(SITE.AVG.CN.RATIO = sum(CN.RATIO*BAS.AREA.FRAC, na.rm=T)/0.7,
              SITE.TYPE = first(type))
  noLargeTreesMySitesAnandData <- anand70PercDomMySites %>% group_by(type) %>%
    summarise(NO.LRG.TREES=n())
  commAvgLittQual <- commAvgdLittQualMySites %>% group_by(SITE.TYPE) %>% 
    summarise(LITT.CTON.MEAN= mean(SITE.AVG.CN.RATIO, na.rm=T),
              LITT.CTON.SE=sd(SITE.AVG.CN.RATIO, na.rm=T)/sqrt(n()))
  commAvgLittQual <- bind_cols(data.frame(NAME=c("Comm. Avg.", "Comm. Avg.")),
                               commAvgLittQual)
}

# Processing ...
{
  # merge aboveground & belowground data point-wise
  cCycle <- belowGround %>% 
    dplyr::select(-SITE.ID, -POINT.TYPE, -POINT.IN.SITE, -POINT.NO, -SOIL.WT) %>%
    left_join(aboveGround, ., by="POINT")
  cCycle <- cCycle %>% 
    separate(col="SITE.ID", into=c("SITE.TYPE", "SITE.NO"), sep=2, remove=F) %>%
    dplyr::select(-SITE.NO)
  cCycle$SITE.TYPE <- as.factor(cCycle$SITE.TYPE)
  
  cCycle <- left_join(cCycle, largeTreeDens, by="POINT")
  
  cCycle <- cCycle %>% mutate(LITTER.NP.RATIO = LITTER.N/LITTER.P) %>%
    mutate(LITTER.WT = LITTER.WT/0.09)
  cCycle$SOIL.SIR[cCycle$SOIL.SIR < 0.18] = NA
  
  # Some data exploration visualization for fitting linear models
  sirPred <- dplyr::select(cCycle, TIP, LITTER.CN.RATIO, LITTER.N, LITTER.WT, LITTER.P)
  socPred <- dplyr::select(cCycle, TIP, LITTER.CN.RATIO, LITTER.N, LITTER.WT, LITTER.P,
                    SOIL.SIR, SOIL.N)  
  
  # Variance inflation factors of predictors
  library(usdm)
  vif(sirPred)
  vif(socPred)
  
  # Linear models of the responses
  socFullLM <- lm(SOIL.C ~ LITTER.CN.RATIO + LITTER.P + LITTER.WT + BAS.AREA + SITE.TYPE, data=cCycle)
  
  sirFullLM <- lm(SOIL.SIR ~ LITTER.CN.RATIO + LITTER.P + LITTER.WT + BAS.AREA + SITE.TYPE, data=cCycle)
  
  # Linear mixed effects models of the responses
  library(lme4)
  soilcLittPLMM    <- lmer(SOIL.C ~ LITTER.P + SITE.TYPE +
                             (1|SITE.ID),
                           data=cCycle, REML=FALSE)
  soilcLittQualLMM <- lmer(SOIL.C ~ LITTER.CN.RATIO + LITTER.P + SITE.TYPE + 
                             (1|SITE.ID),
                           data=cCycle, REML=FALSE)
  soilcStdBMassLMM <- lmer(SOIL.C ~ LITTER.WT + BAS.AREA + SITE.TYPE + 
                             (1|SITE.ID),
                           data=cCycle, REML=FALSE)
  soilcFullLMM     <- lmer(SOIL.C ~ LITTER.CN.RATIO + LITTER.P + LITTER.WT + 
                             BAS.AREA + SITE.TYPE +
                             (1|SITE.ID),
                           data=cCycle, REML=FALSE)
  
#   anova(soilcLittPLMM, soilcLittQualLMM, soilcStandBiomassLMM,
#         soilcLittPForTypeLMM, soilcLittQualForTypeLMM, soilcStandBiomassForTypeLMM,
#         soilcFullLMM, soilcFullForTypeLMM)
#   
#   plot(predict(soilcFullLMM), log(cCycle$SOIL.C[!is.na(cCycle$SOIL.C)]),
#        ylab='y = log(Soil C)', xlab='predicted', main='soilcFullLMM y v/s y_hat')
#   abline(0,1)
#   plot(predict(soilcFullLMM), residuals(soilcFullLMM), ylab='residuals', 
#        xlab='predicted', main='soilcFullLMM residuals v/s predicted')
#   abline(0,0)
#   qqnorm(residuals(soilcFullLMM), main='soilcFullLMM residuals')
#   qqline(residuals(soilcFullLMM))
  
  sirLittPLMM    <- lmer(log(SOIL.SIR) ~ LITTER.P + SITE.TYPE +
                           (1|SITE.ID),
                         data=cCycle, REML=FALSE)
  sirLittQualLMM <- lmer(log(SOIL.SIR) ~ LITTER.CN.RATIO + LITTER.P + SITE.TYPE +
                           (1|SITE.ID),
                           data=cCycle, REML=FALSE)
  sirStdBMassLMM <- lmer(log(SOIL.SIR) ~ LITTER.WT + log(BAS.AREA) + SITE.TYPE +
                           (1|SITE.ID),
                               data=cCycle, REML=FALSE)
  sirFullLMM     <- lmer(log(SOIL.SIR) ~ LITTER.CN.RATIO + LITTER.P + LITTER.WT + 
                         log(BAS.AREA) + SITE.TYPE +
                           (1|SITE.ID),
                       data=cCycle, REML=FALSE)
  
#   anova(sirLittPLMM, sirLittQualLMM, sirStandBiomassLMM,
#         sirLittPForTypeLMM, sirLittQualForTypeLMM, sirStandBiomassForTypeLMM,
#         sirFullLMM, sirFullForTypeLMM)
#   
#   plot(predict(sirFullLMM), log(cCycle$SOIL.SIR[!is.na(cCycle$SOIL.SIR)]),
#        ylab='y = log(SIR)', xlab='predicted', main='sirFullLMM y v/s y_hat')
#   abline(0,1)
#   plot(predict(sirFullLMM), residuals(sirFullLMM), ylab='residuals', 
#        xlab='predicted', main='sirFullLMM residuals v/s predicted')
#   abline(0,0)
#   qqnorm(residuals(sirFullLMM), main='sirFullLMM residuals')
#   qqline(residuals(sirFullLMM))
  
  # Modeling with scaling of variables
  cCycle.S <- data.frame(POINT = cCycle$POINT,
                         SITE.TYPE=cCycle$SITE.TYPE,
                         SITE.ID=cCycle$SITE.ID,
                         SOIL.SIR.S=scale(cCycle$SOIL.SIR), 
                         SOIL.C.S=scale(cCycle$SOIL.C),
                         LITTER.CN.RATIO.S=scale(cCycle$LITTER.CN.RATIO),
                         LITTER.P.S=scale(cCycle$LITTER.P),
                         LITTER.WT.S=scale(cCycle$LITTER.WT),
                         BAS.AREA.S=scale(cCycle$BAS.AREA))

  S.soilcLittPOnlyLMM    <- lmer(SOIL.C.S ~ LITTER.P.S +
                               (1|SITE.ID),
                             data=cCycle.S, REML=FALSE)
  S.soilcLittCNOnlyLMM <- lmer(SOIL.C.S ~ LITTER.CN.RATIO.S +
                             (1|SITE.ID),
                           data=cCycle.S, REML=FALSE)
  S.soilcLittWtOnlyLMM <- lmer(SOIL.C.S ~ LITTER.WT.S +
                             (1|SITE.ID),
                           data=cCycle.S, REML=FALSE)
  S.soilcBasAreaOnlyLMM <- lmer(SOIL.C.S ~ BAS.AREA.S +
                              (1|SITE.ID),
                            data=cCycle.S, REML=FALSE)
  S.soilcLittPLMM    <- lmer(SOIL.C.S ~ LITTER.P.S + SITE.TYPE +
                               (1|SITE.ID),
                             data=cCycle.S, REML=FALSE)
  S.soilcLittCNLMM <- lmer(SOIL.C.S ~ LITTER.CN.RATIO.S + SITE.TYPE + 
                               (1|SITE.ID),
                             data=cCycle.S, REML=FALSE)
  S.soilcLittWtLMM <- lmer(SOIL.C.S ~ LITTER.WT.S + SITE.TYPE + 
                             (1|SITE.ID),
                           data=cCycle.S, REML=FALSE)
  S.soilcBasAreaLMM <- lmer(SOIL.C.S ~ BAS.AREA.S + SITE.TYPE + 
                              (1|SITE.ID),
                            data=cCycle.S, REML=FALSE)
  S.soilcLittPnFragLMM <- lmer(SOIL.C.S ~ (LITTER.P.S*SITE.TYPE) +
                                 (1|SITE.ID),
                               data=cCycle.S, REML=FALSE)
  S.soilcLittCNnFragLMM <- lmer(SOIL.C.S ~ (LITTER.CN.RATIO.S*SITE.TYPE) +
                                 (1|SITE.ID),
                               data=cCycle.S, REML=FALSE)
  S.soilcLittWtnFragLMM <- lmer(SOIL.C.S ~ (LITTER.WT.S*SITE.TYPE) +
                                 (1|SITE.ID),
                               data=cCycle.S, REML=FALSE)
  S.soilcBasAreanFragLMM <- lmer(SOIL.C.S ~ (BAS.AREA.S*SITE.TYPE) +
                                 (1|SITE.ID),
                               data=cCycle.S, REML=FALSE)
  S.soilcLittQualLMM <- lmer(SOIL.C.S ~ LITTER.CN.RATIO.S + LITTER.P.S + SITE.TYPE + 
                                   (1|SITE.ID),
                                 data=cCycle.S, REML=FALSE)
  S.soilcStdBMassLMM <- lmer(SOIL.C.S ~ LITTER.WT.S + BAS.AREA.S + SITE.TYPE + 
                               (1|SITE.ID),
                             data=cCycle.S, REML=FALSE)
  S.soilcFullLMM     <- lmer(SOIL.C.S ~ LITTER.CN.RATIO.S + LITTER.P.S + LITTER.WT.S + 
                               BAS.AREA.S + SITE.TYPE + (1|SITE.ID),
                             data=cCycle.S, REML=FALSE)
  summary(model.avg(S.soilcLittCNLMM, S.soilcLittCNnFragLMM,
                    S.soilcLittPLMM, S.soilcLittPnFragLMM,
                    S.soilcLittWtLMM, S.soilcLittWtnFragLMM,
                    S.soilcBasAreaLMM, S.soilcBasAreanFragLMM))
  summary(model.avg(S.soilcFullLMM, 
                    S.soilcStdBMassLMM, S.soilcLittWtLMM, S.soilcBasAreaLMM,
                    S.soilcLittQualLMM, S.soilcLittCNLMM, S.soilcLittPLMM, rank="AIC"))

  par(mfrow = c(2,2))
  plot(predict(S.soilcLittWtLMM), (cCycle.S$SOIL.C.S[!is.na(cCycle.S$SOIL.C.S)]),
       ylab='observed', xlab='predicted', main='Soil C ~ Litter wt * Frag')
  abline(0,1)
  plot(predict(S.soilcLittWtLMM), residuals(S.soilcLittWtLMM), ylab='residuals', 
       xlab='predicted', main='Soil C ~ Litter wt * Frag')
  abline(0,0)
  qqnorm(residuals(S.soilcLittWtLMM), main='Soil C ~ Litter wt * Frag: residuals Q-Q plot')
  qqline(residuals(S.soilcLittWtLMM))
  
  S.sirLittPOnlyLMM    <- lmer(SOIL.SIR.S ~ LITTER.P.S +
                             (1|SITE.ID),
                           data=cCycle.S, REML=FALSE)
  S.sirLittCNOnlyLMM <- lmer(SOIL.SIR.S ~ LITTER.CN.RATIO.S + 
                           (1|SITE.ID),
                         data=cCycle.S, REML=FALSE)
  S.sirLittWtOnlyLMM <- lmer(SOIL.SIR.S ~ LITTER.WT.S + 
                           (1|SITE.ID),
                         data=cCycle.S, REML=FALSE)
  S.sirBasAreaOnlyLMM <- lmer(SOIL.SIR.S ~ BAS.AREA.S + 
                            (1|SITE.ID),
                          data=cCycle.S, REML=FALSE)
  S.sirLittPLMM    <- lmer(SOIL.SIR.S ~ LITTER.P.S + SITE.TYPE +
                             (1|SITE.ID),
                               data=cCycle.S, REML=FALSE)
  S.sirLittCNLMM <- lmer(SOIL.SIR.S ~ LITTER.CN.RATIO.S + SITE.TYPE + 
                               (1|SITE.ID),
                             data=cCycle.S, REML=FALSE)
  S.sirLittWtLMM <- lmer(SOIL.SIR.S ~ LITTER.WT.S + SITE.TYPE + 
                           (1|SITE.ID),
                         data=cCycle.S, REML=FALSE)
  S.sirBasAreaLMM <- lmer(SOIL.SIR.S ~ BAS.AREA.S + SITE.TYPE + 
                            (1|SITE.ID),
                          data=cCycle.S, REML=FALSE)
  S.sirLittPnFragLMM <- lmer(SOIL.SIR.S ~ (LITTER.P.S*SITE.TYPE) +
                                 (1|SITE.ID),
                               data=cCycle.S, REML=FALSE)
  S.sirLittCNnFragLMM <- lmer(SOIL.SIR.S ~ (LITTER.CN.RATIO.S*SITE.TYPE) +
                                  (1|SITE.ID),
                                data=cCycle.S, REML=FALSE)
  S.sirLittWtnFragLMM <- lmer(SOIL.SIR.S ~ (LITTER.WT.S*SITE.TYPE) +
                                  (1|SITE.ID),
                                data=cCycle.S, REML=FALSE)
  S.sirBasAreanFragLMM <- lmer(SOIL.SIR.S ~ (BAS.AREA.S*SITE.TYPE) +
                                   (1|SITE.ID),
                                 data=cCycle.S, REML=FALSE)
  S.sirLittQualLMM <- lmer(SOIL.SIR.S ~ LITTER.CN.RATIO.S + LITTER.P.S + SITE.TYPE + 
                                 (1|SITE.ID),
                               data=cCycle.S, REML=FALSE)
  S.sirStdBMassLMM <- lmer(SOIL.SIR.S ~ LITTER.WT.S + BAS.AREA.S + SITE.TYPE + 
                               (1|SITE.ID),
                             data=cCycle.S, REML=FALSE)
  S.sirFullLMM     <- lmer(SOIL.SIR.S ~ LITTER.CN.RATIO.S + LITTER.P.S + LITTER.WT.S + 
                               BAS.AREA.S + SITE.TYPE +
                               (1|SITE.ID),
                             data=cCycle.S, REML=FALSE)
  summary(model.avg(S.sirLittCNLMM, S.sirLittPLMM,
                    S.sirLittWtLMM, S.sirBasAreaLMM))
  summary(model.avg(S.sirLittCNLMM, S.sirLittCNnFragLMM,
                    S.sirLittPLMM, S.sirLittPnFragLMM,
                    S.sirLittWtLMM, S.sirLittWtnFragLMM,
                    S.sirBasAreaLMM, S.sirBasAreanFragLMM,
                    S.sirLittPOnlyLMM, S.sirLittCNOnlyLMM,
                    S.sirLittWtOnlyLMM, S.sirBasAreaOnlyLMM))
  summary(model.avg(S.sirFullLMM, 
                    S.sirStdBMassLMM, S.sirLittWtLMM, S.sirBasAreaLMM,
                    S.sirLittQualLMM, S.sirLittCNLMM, S.sirLittPLMM, rank="AIC"))
  
  par(mfrow = c(2,2))
  plot(predict(S.sirLittWtLMM), (cCycle.S$LITTER.WT.S[!is.na(cCycle.S$SOIL.SIR.S)]),
         ylab='observed', xlab='predicted', main='SIR ~ Litter wt * Frag')
  abline(0,1)
  plot(predict(S.sirLittWtLMM), residuals(S.sirLittWtLMM), ylab='residuals', 
       xlab='predicted', main='SIR ~ Litter wt * Frag')
  abline(0,0)
  qqnorm(residuals(S.sirLittWtLMM), main='SIR ~ Litter wt * Frag: residuals Q-Q')
  qqline(residuals(S.sirLittWtLMM))

  
  # Checking if all variables are normally distributed
  par(mfrow = c(4,2))
  qqnorm(cCycle$SOIL.C, main='Soil C')
  qqline(cCycle$SOIL.C)
  qqnorm(cCycle$SOIL.SIR, main='Soil SIR')
  qqline(cCycle$SOIL.SIR)
  qqnorm(log(cCycle$SOIL.C), main='log(Soil C)')
  qqline(log(cCycle$SOIL.C))
  qqnorm(log(cCycle$SOIL.SIR), main='log(Soil SIR)')
  qqline(log(cCycle$SOIL.SIR))
  qqnorm(cCycle$LITTER.P, main='Litter P')
  qqline(cCycle$LITTER.P)
  qqnorm(cCycle$LITTER.CN.RATIO, main='Litter C/N')
  qqline(cCycle$LITTER.CN.RATIO)
  qqnorm(cCycle$LITTER.WT, main='Litter Weight')
  qqline(cCycle$LITTER.WT)
  qqnorm(cCycle$BAS.AREA, main='Basal area')
  qqline(cCycle$BAS.AREA)
  
    
  # To check for interaction
  coplot(SOIL.C ~ SOIL.SIR | LITTER.P, ylab = "SOC (%)", data=cCycle,
         xlab = "SIR (microgram C-CO2)",
         panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
          abline(tmp)
          points(x, y) })
  
  # Scatter plots of responses vs each predictor
  
  p1 <- ggplot(cCycle, aes(x=LITTER.CN.RATIO, y=SOIL.SIR, color=SITE.TYPE)) +
    geom_point(size=2) + xlab("Litter C/N") + ylab("Soil SIR")
  
  p2 <- ggplot(cCycle, aes(x=LITTER.P, y=SOIL.SIR, color=SITE.TYPE)) +
    geom_point(size=2) + xlab("Litter P (%)") + ylab("Soil SIR")
  
  p3 <- ggplot(cCycle, aes(x=LITTER.N, y=SOIL.SIR, color=SITE.TYPE)) +
    geom_point(size=2) + xlab("Litter N (%)") + ylab("Soil SIR")
  
  p4 <- ggplot(cCycle, aes(x=LITTER.WT, y=SOIL.SIR, color=SITE.TYPE)) +
    geom_point(size=2) + xlab("Litter weight (gm)") + ylab("Soil SIR")
  
  p5 <- ggplot(cCycle, aes(x=TIP, y=SOIL.SIR, color=SITE.TYPE)) +
    geom_point(size=2) + xlab("Tree influence potential") + ylab("Soil SIR")
  
  p6 <- ggplot(cCycle, aes(x=LITTER.CN.RATIO, y=SOIL.C, color=SITE.TYPE)) +
    geom_point(size=2) + xlab("Litter C/N") + ylab("Soil organic C (%)")
  
  p7 <- ggplot(cCycle, aes(x=LITTER.P, y=SOIL.C, color=SITE.TYPE)) +
    geom_point(size=2) + xlab("Litter P (%)") + ylab("Soil organic C (%)")
  
  p8 <- ggplot(cCycle, aes(x=LITTER.N, y=SOIL.C, color=SITE.TYPE)) +
    geom_point(size=2) + xlab("Litter N (%)") + ylab("Soil organic C (%)")
  
  p9 <- ggplot(cCycle, aes(x=LITTER.WT, y=SOIL.C, color=SITE.TYPE)) +
    geom_point(size=2) + xlab("Litter weight (gm)") + ylab("Soil organic C (%)")
  
  p10 <- ggplot(cCycle, aes(x=TIP, y=SOIL.C, color=SITE.TYPE)) +
    geom_point(size=2) + xlab("Tree influence potential") + ylab("Soil organic C (%)")
  
  multiplot(p1, p2, p3, p4, p5, cols=3)
  multiplot(p6, p7, p8, p9, p10, cols=3)
  
  # To visualize linear trends, site treatment-wise
  ggplot(cCycle, aes(x=SOIL.SIR, y=SOIL.C, color=SITE.TYPE)) +
    geom_point(size=2) + stat_smooth(method=lm)
  
  ggplot(cCycle, aes(x=LITTER.CN.RATIO, y=SOIL.C, color=SITE.TYPE)) +
    geom_point(size=2)
  
  ggplot(cCycle, aes(x=LITTER.WT, y=TIP, color=SITE.TYPE)) +
    geom_point(size=2)
  
  ggplot(cCycle, aes(x=SOIL.SIR, y=ADJ.SOIL.C, color=SITE.TYPE)) +
    geom_point(size=2) + stat_smooth(method=lm)
  
  # Interaction checking
  coplot(SOIL.SIR ~ LITTER.CN.RATIO | TIP, ylab = "SIR", data=cCycle,
                 xlab = "TIP",
                 panel = function(x, y, ...) {
                       tmp <- lm(y ~ x, na.action = na.omit)
                       abline(tmp)
                       points(x, y) })
  
  # Boxplots for broad-level visualization.
  p20 <- cCycle %>% ggplot(aes(y=LITTER.CN.RATIO, x=SITE.ID, color=SITE.TYPE)) +
    geom_boxplot()+ ylab("Litter C/N")
  p21 <- cCycle %>% ggplot(aes(y=LITTER.N, x=SITE.ID, color=SITE.TYPE)) +
    geom_boxplot()+ ylab("Litter N (% DM)")
  p22 <- cCycle %>% ggplot(aes(y=LITTER.P, x=SITE.ID, color=SITE.TYPE)) +
    geom_boxplot()+ ylab("Litter P (% DM)")
  p23 <- cCycle %>% ggplot(aes(y=LITTER.WT, x=SITE.ID, color=SITE.TYPE)) +
    geom_boxplot()+ ylab("Litter weight (gm)")
  p24 <- cCycle %>% ggplot(aes(y=TIP, x=SITE.ID, color=SITE.TYPE)) +
    geom_boxplot()+ ylab("Tree Influence Potential")
  p25 <- cCycle %>% ggplot(aes(y=SOIL.C, x=SITE.ID, color=SITE.TYPE)) +
    geom_boxplot()+ ylab("Soil C (% DM)")
  p26 <- cCycle %>% ggplot(aes(y=SOIL.SIR, x=SITE.ID, color=SITE.TYPE)) +
    geom_boxplot()+ ylab("Soil SIR (mu gm C-CO2 / hr gm)")
    
  multiplot(p25, p20, p21, p24, p26, p22, p23, cols=2)

  # Comparing average values across FR and CT
  soilcMeansPerSite <- cCycle %>% group_by(SITE.ID) %>% 
    summarise(SITE.MEAN = mean(SOIL.C), SITE.TYPE=first(SITE.TYPE)) 
  t.test(soilcMeansPerSite$SITE.MEAN~soilcMeansPerSite$SITE.TYPE)
  soilcMeans <- soilcMeansPerSite %>% group_by(SITE.TYPE) %>%
    summarise(MEAN=mean(SITE.MEAN), SE=sd(SITE.MEAN)/sqrt(n()))
  ggplot(soilcMeans, aes(y=MEAN, x=SITE.TYPE, fill=SITE.TYPE)) + theme_classic() +
    geom_bar(position="dodge", stat="identity", width=0.4) + guides(fill=F) +
    geom_errorbar(aes(ymin=MEAN-(SE),
                      ymax=MEAN+(SE)), width=.07,
                  position=position_dodge(0.9)) +
    scale_x_discrete(labels=
                       c("Contiguous forests \n (N = 6 sites)", "Fragmented forests \n (N = 6 sites)")) +    
    ylab("Soil %C \n") +
    theme(legend.title=element_blank(),
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 18)) +
    scale_fill_grey(start=0.5)
  ggsave("figs/SoilC_CTvsFRCompare.png")
  soilcMeans <- bind_cols(data.frame(QUANT=c("Soil C", "Soil C")), soilcMeans)

  sirMeansPerSite <- cCycle %>% group_by(SITE.ID) %>% 
    summarise(SITE.MEAN = mean(SOIL.SIR, na.rm=T),
              SITE.TYPE = first(SITE.TYPE))
  t.test(sirMeansPerSite$SITE.MEAN~sirMeansPerSite$SITE.TYPE)
  sirMeans <- sirMeansPerSite %>% group_by(SITE.TYPE) %>%
    summarise(MEAN = mean(SITE.MEAN), SE = sd(SITE.MEAN)/sqrt(n()))
  ggplot(sirMeans, aes(y=MEAN, x=SITE.TYPE, fill=SITE.TYPE)) + theme_classic() +
    geom_bar(position="dodge", stat="identity", width=0.4) + guides(fill=F) +
    geom_errorbar(aes(ymin=MEAN-(SE),
                      ymax=MEAN+(SE)), width=.07,
                  position=position_dodge(0.9)) +
    scale_x_discrete(labels=
                       c("Contiguous forests \n (N = 6 sites)", "Fragmented forests \n (N = 6 sites)")) +    
#     ylab(expression('Soil SIR ('*mu gm ~C-CO[2]~ g^{-1} dry wt. soil hr^{-1}*')')) +
    ylab("Soil SIR \n (micro gm C-CO2 \n per gm dry wt. soil per hr) \n") +
    theme(legend.title=element_blank(),
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 18)) +
    scale_fill_grey(start=0.5)
  ggsave("figs/SIR_CTvsFRCompare.png")
  sirMeans <- bind_cols(data.frame(QUANT=c("Soil SIR", "Soil SIR")), sirMeans)

  soilCompare <- bind_rows(soilcMeans, sirMeans)

  littWtMeansPerSite <- cCycle %>% group_by(SITE.ID) %>% 
    summarise(SITE.MEAN = mean(LITTER.WT), SITE.TYPE=first(SITE.TYPE)) 
  t.test(littWtMeansPerSite$SITE.MEAN~littWtMeansPerSite$SITE.TYPE)
  littWtMeans <- littWtMeansPerSite %>% group_by(SITE.TYPE) %>%
    summarise(MEAN=mean(SITE.MEAN), SE=sd(SITE.MEAN)/sqrt(n()))
  ggplot(littWtMeans, aes(y=MEAN, x=SITE.TYPE, fill=SITE.TYPE)) + theme_classic() +
    geom_bar(position="dodge", stat="identity", width=0.4) + guides(fill=F) +
    geom_errorbar(aes(ymin=MEAN-(SE),
                      ymax=MEAN+(SE)), width=.07,
                  position=position_dodge(0.9)) +
    scale_x_discrete(labels=
                       c("Contiguous \n (N = 6 sites)", "Fragment \n (N = 6 sites)")) +    
    ylab("Litter weight (gm per m.sq)\n") +
    theme(legend.title=element_blank(),
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 24),
          axis.text = element_text(size = 24)) +
    scale_fill_grey(start=0.5)
  ggsave("figs/LittWt_CTvsFRCompare.png")
  
  littCNMeansPerSite <- cCycle %>% group_by(SITE.ID) %>% 
    summarise(SITE.MEAN = mean(LITTER.CN.RATIO), SITE.TYPE=first(SITE.TYPE)) 
  t.test(littCNMeansPerSite$SITE.MEAN~littCNMeansPerSite$SITE.TYPE)
  littCNMeans <- littCNMeansPerSite %>% group_by(SITE.TYPE) %>%
    summarise(MEAN=mean(SITE.MEAN), SE=sd(SITE.MEAN)/sqrt(n()))
  ggplot(littCNMeans, aes(y=MEAN, x=SITE.TYPE, fill=SITE.TYPE)) + theme_classic() +
    geom_bar(position="dodge", stat="identity", width=0.4) + guides(fill=F) +
    geom_errorbar(aes(ymin=MEAN-(SE),
                      ymax=MEAN+(SE)), width=.07,
                  position=position_dodge(0.9)) +
    scale_x_discrete(labels=
                       c("Contiguous \n (N = 6 sites)", "Fragment \n (N = 6 sites)")) +    
    ylab("Litter C/N \n") +
    theme(legend.title=element_blank(),
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 24),
          axis.text = element_text(size = 24)) +
    scale_fill_grey(start=0.5)
  ggsave("figs/LittCN_CTvsFRCompare.png")

  basAreaMeansPerSite <- cCycle %>% group_by(SITE.ID) %>% 
    summarise(SITE.MEAN = mean(LITTER.CN.RATIO), SITE.TYPE=first(SITE.TYPE)) 
  t.test(basAreaMeansPerSite$SITE.MEAN~basAreaMeansPerSite$SITE.TYPE)
  basAreaMeans <- basAreaMeansPerSite %>% group_by(SITE.TYPE) %>%
    summarise(MEAN=mean(SITE.MEAN), SE=sd(SITE.MEAN)/sqrt(n()))
  ggplot(basAreaMeans, aes(y=MEAN, x=SITE.TYPE, fill=SITE.TYPE)) + theme_classic() +
    geom_bar(position="dodge", stat="identity", width=0.4) + guides(fill=F) +
    geom_errorbar(aes(ymin=MEAN-(SE),
                      ymax=MEAN+(SE)), width=.07,
                  position=position_dodge(0.9)) +
    scale_x_discrete(labels=
                       c("Contiguous \n (N = 6 sites)", "Fragment \n (N = 6 sites)")) +    
    ylab("Basal area (m. sq. per ha) \n") +
    theme(legend.title=element_blank(),
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 24),
          axis.text = element_text(size = 24)) +
    scale_fill_grey(start=0.5)
  ggsave("figs/BasArea_CTvsFRCompare.png")

  littPMeansPerSite <- cCycle %>% group_by(SITE.ID) %>% 
    summarise(SITE.MEAN = mean(LITTER.P), SITE.TYPE=first(SITE.TYPE)) 
  t.test(littPMeansPerSite$SITE.MEAN~littPMeansPerSite$SITE.TYPE)
  littPMeans <- littPMeansPerSite %>% group_by(SITE.TYPE) %>%
    summarise(MEAN=mean(SITE.MEAN), SE=sd(SITE.MEAN)/sqrt(n()))
  ggplot(littPMeans, aes(y=MEAN, x=SITE.TYPE, fill=SITE.TYPE)) + theme_classic() +
    geom_bar(position="dodge", stat="identity", width=0.4) + guides(fill=F) +
    geom_errorbar(aes(ymin=MEAN-(SE),
                      ymax=MEAN+(SE)), width=.07,
                  position=position_dodge(0.9)) +
    scale_x_discrete(labels=
                       c("Contiguous \n (N = 6 sites)", "Fragment \n (N = 6 sites)")) +    
    ylab("Litter %P \n") +
    theme(legend.title=element_blank(),
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 24),
          axis.text = element_text(size = 24)) +
    scale_fill_grey(start=0.5)
  ggsave("figs/LittP_CTvsFRCompare.png")

  # Responses, site-wise boxplots
  ggplot(cCycle, aes(y=SOIL.C, x=SITE.ID, fill=SITE.TYPE)) +
    theme_classic() +
    geom_boxplot()+ ylab("Soil %C \n") + xlab("\n Site ID") +
    theme(legend.title=element_blank(),
          legend.text = element_text(size = 15),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 15),
          axis.text.x = element_text(angle=45, hjust=1)) +
    scale_fill_grey(labels=c("Contiguous", "Fragment"), start=0.5)
  ggsave("figs/SoilCvsSites.png")
  ggplot(cCycle, aes(y=SOIL.SIR, x=SITE.ID, fill=SITE.TYPE)) +
    theme_classic() +
    geom_boxplot() +
    ylab("Soil SIR \n (micro gm C-CO2 \n per gm dry wt. soil per hr) \n") +
    xlab("\n Site ID") +
    theme(legend.title=element_blank(),
          legend.text = element_text(size = 15),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 15),
          axis.text.x = element_text(angle=45, hjust=1)) +
    scale_fill_grey(labels=c("Contiguous", "Fragment"), start=0.5)
  ggsave("figs/SIRvsSites.png")

  # To plot SIR and soil C v/s litter C/N and P, but site-wise with SEs
  soilNLittMeans <- cCycle %>% group_by(SITE.ID) %>% 
    summarise(SITE.TYPE=first(SITE.TYPE),
              SOIL.C.MEAN = mean(SOIL.C),
              SOIL.C.SE = sd(SOIL.C)/sqrt(n()),
              SIR.MEAN = mean(SOIL.SIR, na.rm=T),
              SIR.SE = sd(SOIL.SIR, na.rm=T)/sqrt(n()),
              LITTER.CN.MEAN = mean(LITTER.CN.RATIO),
              LITTER.CN.SE = sd(LITTER.CN.RATIO)/sqrt(n()),
              LITTER.P.MEAN = mean(LITTER.P),
              LITTER.P.SE = sd(LITTER.P)/sqrt(n()),
              LITTER.WT.MEAN = mean(LITTER.WT),
              LITTER.WT.SE = sd(LITTER.WT)/sqrt(n()),
              BAS.AREA.MEAN = mean(BAS.AREA),
              BAS.AREA.SE = sd(BAS.AREA)/sqrt(n()))
  ggplot(soilNLittMeans, aes(y=SOIL.C.MEAN, x=LITTER.CN.MEAN,
                                 color=SITE.TYPE)) + 
    theme_classic() + geom_point(size=4) +
    geom_errorbarh(aes(xmin=LITTER.CN.MEAN-(LITTER.CN.SE),
                       xmax=LITTER.CN.MEAN+(LITTER.CN.SE))) + 
    geom_errorbar(aes(ymin=SOIL.C.MEAN-(SOIL.C.SE),
                      ymax=SOIL.C.MEAN+(SOIL.C.SE)), width=0) +
    xlab("\n Litter C/N") + ylab("Soil %C \n") +
    theme(legend.title=element_blank(),
          legend.text = element_text(size = 15),
          legend.position = c(0.85, 0.95),
          #legend.background = element_rect(color="grey"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18)) +
    scale_color_grey(labels=c("Contiguous", "Fragment"))
  ggsave("figs/SoilCvsLitterCN_SiteAvgs.png")
  ggplot(soilNLittMeans, aes(y=SIR.MEAN, x=LITTER.CN.MEAN,
                                 color=SITE.TYPE)) + 
    theme_classic() + geom_point(size=4) +
    geom_errorbarh(aes(xmin=LITTER.CN.MEAN-(LITTER.CN.SE),
                       xmax=LITTER.CN.MEAN+(LITTER.CN.SE))) + 
    geom_errorbar(aes(ymin=SIR.MEAN-(SIR.SE),
                      ymax=SIR.MEAN+(SIR.SE)), width=0) +
    xlab("\n Litter C/N") +
    ylab("Soil SIR \n (micro gm C-CO2 \n per gm dry wt. soil per hr) \n") +
    theme(legend.title=element_blank(),
          legend.text = element_text(size = 15),
          legend.position = c(0.85, 0.55),
          #legend.background = element_rect(color="grey"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18)) +
    scale_color_grey(labels=c("Contiguous", "Fragment"))
  ggsave("figs/SIRvsLitterCN_SiteAvgs.png")
  ggplot(soilNLittMeans, aes(y=SOIL.C.MEAN, x=LITTER.P.MEAN,
                                 color=SITE.TYPE)) + 
    theme_classic() + geom_point(size=4) +
    geom_errorbarh(aes(xmin=LITTER.P.MEAN-(LITTER.P.SE),
                       xmax=LITTER.P.MEAN+(LITTER.P.SE))) + 
    geom_errorbar(aes(ymin=SOIL.C.MEAN-(SOIL.C.SE),
                      ymax=SOIL.C.MEAN+(SOIL.C.SE)), width=0) +
    xlab("\n Litter %P") + ylab("Soil %C \n") +
    theme(legend.title=element_blank(),
          legend.text = element_text(size = 15),
          legend.position = c(0.85, 0.15),
          #legend.background = element_rect(color="grey"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18)) +
    scale_color_grey(labels=c("Contiguous", "Fragment"))
  ggsave("figs/SoilCvsLitterP_SiteAvgs.png")
  ggplot(soilNLittMeans, aes(y=SIR.MEAN, x=LITTER.P.MEAN,
                                 color=SITE.TYPE)) + 
    theme_classic() + geom_point(size=4) +
    geom_errorbarh(aes(xmin=LITTER.P.MEAN-(LITTER.P.SE),
                       xmax=LITTER.P.MEAN+(LITTER.P.SE))) + 
    geom_errorbar(aes(ymin=SIR.MEAN-(SIR.SE),
                      ymax=SIR.MEAN+(SIR.SE)), width=0) +
    xlab("\n Litter %P") +
    ylab("Soil SIR \n (micro gm C-CO2 \n per gm dry wt. soil per hr) \n") +
    theme(legend.title=element_blank(),
          legend.text = element_text(size = 15),
          legend.position = c(0.85, 0.25),
          #legend.background = element_rect(color="grey"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18)) +
    scale_color_grey(labels=c("Contiguous", "Fragment"))
  ggsave("figs/SIRvsLitterP_SiteAvgs.png")
  multiplot(p91, p93, p92, p94, cols=2)
  ggplot(soilNLittMeans, aes(y=SOIL.C.MEAN, x=LITTER.WT.MEAN,
                             color=SITE.TYPE)) + 
    theme_classic() + geom_point(size=4) +
    geom_errorbarh(aes(xmin=LITTER.WT.MEAN-(LITTER.WT.SE),
                       xmax=LITTER.WT.MEAN+(LITTER.WT.SE))) + 
    geom_errorbar(aes(ymin=SOIL.C.MEAN-(SOIL.C.SE),
                      ymax=SOIL.C.MEAN+(SOIL.C.SE)), width=0) +
    xlab("\n Litter weight (gm per m.sq)") + ylab("Soil %C \n") +
    theme(legend.title=element_blank(),
          legend.text = element_text(size = 15),
          legend.position = c(0.85, 0.9),
          #legend.background = element_rect(color="grey"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18)) +
    scale_color_grey(labels=c("Contiguous", "Fragment"))
  ggsave("figs/SoilCvsLitterWt_SiteAvgs.png")
  ggplot(soilNLittMeans, aes(y=SIR.MEAN, x=LITTER.WT.MEAN,
                             color=SITE.TYPE)) + 
    theme_classic() + geom_point(size=4) +
    geom_errorbarh(aes(xmin=LITTER.WT.MEAN-(LITTER.WT.SE),
                       xmax=LITTER.WT.MEAN+(LITTER.WT.SE))) + 
    geom_errorbar(aes(ymin=SIR.MEAN-(SIR.SE),
                      ymax=SIR.MEAN+(SIR.SE)), width=0) +
    xlab("\n Litter weight (gm per m.sq)") +
    ylab("Soil SIR \n (micro gm C-CO2 \n per gm dry wt. soil per hr) \n") +
    theme(legend.title=element_blank(),
          legend.text = element_text(size = 15),
          legend.position = c(0.85, 0.9),
          #legend.background = element_rect(color="grey"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18)) +
    scale_color_grey(labels=c("Contiguous", "Fragment"))
  ggsave("figs/SIRvsLitterWt_SiteAvgs.png")
  ggplot(soilNLittMeans, aes(y=SOIL.C.MEAN, x=BAS.AREA.MEAN,
                             color=SITE.TYPE)) + 
    theme_classic() + geom_point(size=4) +
    geom_errorbarh(aes(xmin=BAS.AREA.MEAN-(BAS.AREA.SE),
                       xmax=BAS.AREA.MEAN+(BAS.AREA.SE))) + 
    geom_errorbar(aes(ymin=SOIL.C.MEAN-(SOIL.C.SE),
                      ymax=SOIL.C.MEAN+(SOIL.C.SE)), width=0) +
    xlab("\n Basal area (m sq. per ha)") + ylab("Soil %C \n") +
    theme(legend.title=element_blank(),
          legend.text = element_text(size = 15),
          legend.position = c(0.9, 0.9),
          #legend.background = element_rect(color="grey"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18)) +
    scale_color_grey(labels=c("Contiguous", "Fragment"))
  ggsave("figs/SoilCvsBasArea_SiteAvgs.png")
  ggplot(soilNLittMeans, aes(y=SIR.MEAN, x=BAS.AREA.MEAN,
                             color=SITE.TYPE)) + 
    theme_classic() + geom_point(size=4) +
    geom_errorbarh(aes(xmin=BAS.AREA.MEAN-(BAS.AREA.SE),
                       xmax=BAS.AREA.MEAN+(BAS.AREA.SE))) + 
    geom_errorbar(aes(ymin=SIR.MEAN-(SIR.SE),
                      ymax=SIR.MEAN+(SIR.SE)), width=0) +
    xlab("\n Basal area (m sq. per ha)") +
    ylab("Soil SIR \n (micro gm C-CO2 \n per gm dry wt. soil per hr) \n") +
    theme(legend.title=element_blank(),
          legend.text = element_text(size = 15),
          legend.position = c(0.85, 0.5),
          #legend.background = element_rect(color="grey"),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 18)) +
    scale_color_grey(labels=c("Contiguous", "Fragment"))
  ggsave("figs/SIRvsBasArea_SiteAvgs.png")
  # To see if log-transform helps with normality assumption of variables
  p31 <- ggplot(cCycle, aes(x=log(LITTER.CN.RATIO), y=log(SOIL.SIR), color=SITE.TYPE)) +
    geom_point(size=2) + xlab("log(Litter C/N)") + ylab("log(Soil SIR)")
  
  p32 <- ggplot(cCycle, aes(x=log(LITTER.P), y=log(SOIL.SIR), color=SITE.TYPE)) +
    geom_point(size=2) + xlab("log(Litter P (%))") + ylab("log(Soil SIR)")
  
  p33 <- ggplot(cCycle, aes(x=log(LITTER.N), y=log(SOIL.SIR), color=SITE.TYPE)) +
    geom_point(size=2) + xlab("log(Litter N (%))") + ylab("log(Soil SIR)")
  
  p34 <- ggplot(cCycle, aes(x=log(LITTER.WT), y=log(SOIL.SIR), color=SITE.TYPE)) +
    geom_point(size=2) + xlab("log(Litter weight (gm))") + ylab("log(Soil SIR)")
  
  p35 <- ggplot(cCycle, aes(x=log(TIP), y=log(SOIL.SIR), color=SITE.TYPE)) +
    geom_point(size=2) + xlab("log(TIP)") + ylab("log(Soil SIR)")
  
  p36 <- ggplot(cCycle, aes(x=log(LITTER.CN.RATIO), y=log(SOIL.C), color=SITE.TYPE)) +
    geom_point(size=2) + xlab("log(Litter C/N)") + ylab("log(Soil organic C (%))")
  
  p37 <- ggplot(cCycle, aes(x=log(LITTER.P), y=log(SOIL.C), color=SITE.TYPE)) +
    geom_point(size=2) + xlab("log(Litter P (%))") + ylab(log("Soil organic C (%))")
  
  p38 <- ggplot(cCycle, aes(x=log(LITTER.N), y=log(SOIL.C), color=SITE.TYPE)) +
    geom_point(size=2) + xlab("log(Litter N (%))") + ylab("log(Soil organic C (%))")
  
  p39 <- ggplot(cCycle, aes(x=log(LITTER.WT), y=log(SOIL.C), color=SITE.TYPE)) +
    geom_point(size=2) + xlab("log(Litter weight (gm))") + ylab("log(Soil organic C (%))")
  
  p30 <- ggplot(cCycle, aes(x=log(TIP), y=log(SOIL.C), color=SITE.TYPE)) +
    geom_point(size=2) + xlab("log(TIP)") + ylab("log(Soil organic C (%))")
  
  multiplot(p31, p32, p33, p34, p35, cols=3)
  multiplot(p36, p37, p38, p39, p30, cols=3)
  
  # all sites separately
#   p1 <- ggplot(cCycle, aes(x=LITTER.CN.RATIO, y=SOIL.SIR, color=SITE.TYPE)) +
#     geom_point(size=2)
#   p1 + facet_grid(SITE.TYPE ~ SITE.ID)
  
  # plot fit-lines
  ggplot(cCycle, aes(x=LITTER.CN.RATIO, y=SOIL.SIR, color=SITE.TYPE)) +
    geom_point(size=2) + stat_smooth(method=lm)
  ggplot(cCycle, aes(x=LITTER.CN.RATIO, y=SOIL.C, color=SITE.TYPE)) +
    geom_point(size=2) + stat_smooth(method=lm)
  
  # To see if basal area and tip are correlated
  ggplot(cCycle, aes(x=TIP, y=BAS.AREA, color=SITE.TYPE)) +
    geom_point(size=2) + stat_smooth(method=lm)

  # Comparing community avg v/s my sample avg litter qual
  sampleAvgLittQual <- cCycle %>% filter(SITE.ID == "CT1" |       # Kokka 1
                                           SITE.ID == "CT2" |     # Kokka 2
                                           SITE.ID == "FR1" |     # Rudraguppe
                                           SITE.ID == "FR3" |     # Arji
                                           SITE.ID == "FR6") %>%  # Arapattu
    group_by(SITE.ID) %>% 
    summarise(MEAN= mean(LITTER.CN.RATIO, na.rm=T),
              SE=sd(LITTER.CN.RATIO, na.rm=T)/sqrt(n()),
              SITE.TYPE=first(SITE.TYPE)) %>% group_by(SITE.TYPE) %>%
    summarise(LITT.CTON.MEAN=mean(MEAN), 
              LITT.CTON.SE=sd(MEAN, na.rm=T)/sqrt(n()))
  levels(sampleAvgLittQual$SITE.TYPE)[levels(sampleAvgLittQual$SITE.TYPE)=="FR"] <-
    "FRAGMENT"
  levels(sampleAvgLittQual$SITE.TYPE)[levels(sampleAvgLittQual$SITE.TYPE)=="CT"] <-
    "CONTROL"
  sampleAvgLittQual <- bind_cols(data.frame(NAME=c("Sample Avg.", "Sample Avg.")),
                                 sampleAvgLittQual)
  littQualCompare <- bind_rows(commAvgLittQual, sampleAvgLittQual)
  ggplot(littQualCompare, aes(y=LITT.CTON.MEAN, x=NAME, fill=SITE.TYPE)) +
    theme_classic() +
    geom_bar(position="dodge", stat="identity") +
    geom_errorbar(aes(ymin=LITT.CTON.MEAN-(LITT.CTON.SE),
                      ymax=LITT.CTON.MEAN+(LITT.CTON.SE)), width=.1,
                  position=position_dodge(0.9)) +
    xlab("Estimate type") + ylab("Average litter C/N \n") +
    theme(legend.title=element_blank(),
          legend.text = element_text(size = 15),
          legend.position = c(0.15, 0.9),
          axis.title.x=element_blank(), 
          axis.ticks=element_blank(),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 18)) +
    scale_x_discrete(labels=
                       c("Community \n weighted-average", "Sample average")) +
    scale_fill_grey(labels=c("Contiguous", "Fragment"), start=0.5)
  ggsave("figs/LittQualCompare_CommVsSampleAvg.png")
  
#   cCycle.l1g <- cCycle %>% group_by(SITE.TYPE) %>%
#     do(glance(fitl1 <- lm(SOIL.SIR ~ LITTER.CN.RATIO+LITTER.N, data=., na.action=na.omit)))
# 
#   ggplot(cCycle, aes(x=LITTER.CN.RATIO, y=SOIL.SIR, color=SITE.TYPE)) +
#     geom_point() + ggplot(cCycle.l1, aes(x=LITTER.CN.RATIO, y=SOIL.SIR, color=SITE.TYPE)) +
#     geom_line()

  # Simple visualization of Anand's data
  anandLitter %>% ggplot(aes(y=Litter.CtoN, x=Forest.name, color=Forest.type)) +
    geom_boxplot()+ ylab("Litter C/N")
  anandLitter %>% ggplot(aes(y=Soil.C, x=Forest.name, color=Forest.type)) +
    geom_boxplot()+ ylab("Soil C")
  ggplot(anandLitter, aes(x=Litter.CtoN, y=Soil.C, color=Forest.type)) +
    geom_point(size=2)

  anandLitter %>% filter(Forest.name == "Arapattu" |
                           Forest.name == "Karnad" |
                           Forest.name == "Betolli" |
                           Forest.name == "Arji" |
                           Forest.name == "Perumbadi") %>% 
    ggplot(aes(x=Litter.CtoN, y=Soil.C, color=Forest.type)) + geom_point(size=2)

  # Size class-wise plotting of tree community
  treeDensS1 <- filter(treeDensBySizeClass, SIZE.CLASS=="S1")
  t.test(treeDensS1$TREE.DENS ~ treeDensS1$SITE.TYPE)
  t.test(treeDensS1$PROP.STEMS ~ treeDensS1$SITE.TYPE)
  treeDensS2 <- filter(treeDensBySizeClass, SIZE.CLASS=="S2")
  t.test(treeDensS2$TREE.DENS ~ treeDensS2$SITE.TYPE)
  t.test(treeDensS2$PROP.STEMS ~ treeDensS2$SITE.TYPE)
  treeDensS3 <- filter(treeDensBySizeClass, SIZE.CLASS=="S3")
  t.test(treeDensS3$TREE.DENS ~ treeDensS3$SITE.TYPE)
  t.test(treeDensS3$PROP.STEMS ~ treeDensS3$SITE.TYPE)
  treeDensS4 <- filter(treeDensBySizeClass, SIZE.CLASS=="S4")
  t.test(treeDensS4$TREE.DENS ~ treeDensS4$SITE.TYPE)
  t.test(treeDensS4$PROP.STEMS ~ treeDensS4$SITE.TYPE)
  treeDensS5 <- filter(treeDensBySizeClass, SIZE.CLASS=="S5")
  t.test(treeDensS5$TREE.DENS ~ treeDensS5$SITE.TYPE)
  t.test(treeDensS5$PROP.STEMS ~ treeDensS5$SITE.TYPE)

  p41 <- ggplot(treeDensMeans, aes(y=MEAN, x=SIZE.CLASS, fill=SITE.TYPE)) +
    theme_classic() +
    geom_bar(position="dodge", stat="identity") +
    geom_errorbar(aes(ymin=MEAN-SE, ymax=MEAN+SE), width=.2, position=position_dodge(0.9)) +
    scale_x_discrete(labels=c(paste(binBounds[1]*100,binBounds[2]*100, sep=" - "),
                              paste(binBounds[2]*100,binBounds[3]*100, sep=" - "),
                              paste(binBounds[3]*100,binBounds[4]*100, sep=" - "),
                              paste(binBounds[4]*100,binBounds[5]*100, sep=" - "),
                              paste(binBounds[5]*100,binBounds[6]*100, sep=" - "))) +
    xlab("\n DBH (cm) class") + ylab("Stems/ha \n") + ylim(0,200) +
    theme(legend.position=c(0.85, 0.85),
          legend.title=element_blank(),
          legend.text = element_text(size = 15),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 15)) +
    annotate("text", x = treeDensMeans$SIZE.CLASS[1], 
             y = treeDensMeans$MEAN[1] + (1.5*treeDensMeans$SE[1]),
             label = "*", size=8) +
    annotate("text", x = treeDensMeans$SIZE.CLASS[2], 
             y = treeDensMeans$MEAN[2] + (1.5*treeDensMeans$SE[2]),
             label = "*", size=8) +
    annotate("text", x = treeDensMeans$SIZE.CLASS[3], 
             y = treeDensMeans$MEAN[3] + (1.5*treeDensMeans$SE[3]),
             label = "*", size=8) +
    scale_fill_grey(labels=c("Contiguous", "Fragment"), start=0.5)
  ggsave("figs/TreeSizeClassDistr_StemsPerHa.png")
  p42 <- ggplot(treeStemPropsMeans, aes(y=MEAN, x=SIZE.CLASS, fill=SITE.TYPE)) +
    theme_classic() +
    geom_bar(position="dodge", stat="identity") + 
    geom_errorbar(aes(ymin=MEAN-SE, ymax=MEAN+SE), width=.2, position=position_dodge(0.9)) +
    scale_x_discrete(labels=c(paste(binBounds[1]*100,binBounds[2]*100, sep=" - "),
                              paste(binBounds[2]*100,binBounds[3]*100, sep=" - "),
                              paste(binBounds[3]*100,binBounds[4]*100, sep=" - "),
                              paste(binBounds[4]*100,binBounds[5]*100, sep=" - "),
                              paste(binBounds[5]*100,binBounds[6]*100, sep=" - "))) +
    xlab("\n DBH (cm) class") + ylab("Proportion of stems \n") + ylim(0,0.4) +
    theme(legend.position=c(0.85, 0.85),
          legend.title=element_blank(),
          legend.text = element_text(size = 15),
          axis.line.x = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.line.y = element_line(colour = 'black', size=0.5, linetype='solid'),
          axis.title = element_text(size = 18),
          axis.text = element_text(size = 15)) +
    annotate("text", x = treeStemPropsMeans$SIZE.CLASS[9], 
             y = treeStemPropsMeans$MEAN[9] + (1.5*treeStemPropsMeans$SE[9]),
             label = "**", size=8) +
    scale_fill_grey(labels=c("Contiguous", "Fragment"), start=0.5)
  ggsave("figs/TreeSizeClassDistr_PropStems.png")
  multiplot(p41, p42, cols=1)
}


multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}