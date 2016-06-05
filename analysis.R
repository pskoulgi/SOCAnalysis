library(dplyr)
library(tidyr)
library(reshape2)

# Read-in soil data
{
  soilCN <- read.table("F:/Workspace/MSc/Project/soc-fragments/Pradeep_CoorgSoilCN.csv",
                       as.is=T, header=T, sep=",")
  
  # Clean-up col names, drop needless cols
  soilCN <- select(soilCN, c(SAMPLE, MASS, C_PERC, N_PERC))
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
  soilCN <- select(soilCN, -DUMMY)
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
  soilSIR <- soilSIR %>% rename(POINT=SAMPLE) %>%
    select(POINT, C.RELEASED)
  
  belowGround <- left_join(soilCN, soilSIR, by="POINT")
  belowGround <- rename(belowGround, SOIL.SIR=C.RELEASED)
  
  belowGround <- separate(belowGround, col="POINT", into=c("SITE.ID", "POINT.TYPE", "POINT.ID"), sep=c(3,4), remove=F)
  belowGround <- separate(belowGround, col="SITE.ID", into=c("SITE.TYPE", "SITE.NO"), sep=2, remove=F)
  belowGround$SITE.TYPE <- as.factor(belowGround$SITE.TYPE)
  belowGround$SITE.ID <- as.factor(belowGround$SITE.ID)
}

# Read-in aboveground data
{
  # Litter C, N data
  litterCN <- read.table("F:/Workspace/MSc/Project/soc-fragments/Pradeep_CoorgLitterCN.csv",
                         as.is=T, header=T, sep=",")
  
  litterCN <- select(litterCN, c(SAMPLE, C_PERC, N_PERC))
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
  
  litterP <- select(litterP, c(SAMPLE, P_PERC))
  litterP <- litterP %>% rename(POINT=SAMPLE) %>% rename(LITTER.P=P_PERC) %>%
    filter(grepl("CT|FR", POINT))
  
  # Litter weights
  litterWts <- read.table("F:/Workspace/MSc/Project/soc-fragments/litterProperties.csv",
                          as.is=T, header=T, sep=",")
  litterWts$AVGWT <- rowMeans(litterWts[ ,c("DRYWT1", "DRYWT2", "DRYWT3")], na.rm=TRUE)
  litterAvgWt <- litterWts %>% select(c(POINT, AVGWT)) %>% rename(LITTER.WT = AVGWT)
  
  litter <- left_join(litterCN, litterP, by="POINT")
  litter <- left_join(litter, litterAvgWt, by="POINT")
  
  # Tree influence potential
  treeDistrib <- read.table("F:/Workspace/MSc/Project/soc-fragments/treeInfluencePotentialData.csv",
                            as.is=T, header=T, sep=",")
  
  gbhs <- select(treeDistrib, c(GBH1, GBH2, GBH3, GBH4))
  
  treeDistrib$DBH <- (gbhs/pi)^2 %>% rowSums(na.rm=T) %>% sqrt()
  
  treeDistrib <- filter(treeDistrib, DBH>=0.1)
  
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
  
  aboveGround <- left_join(litter, basArea, by="POINT")
  aboveGround <- left_join(aboveGround, treeDens, by="POINT")
    
}

# Processing ...
{
  # merge aboveground & belowground data point-wise
  cCycle <- left_join(agParams, soilSIRC, by="POINT")
  
  cCycle <- separate(cCycle, col="POINT", 
                     into=c("SITEID", "POINTTYPE", "POINTID"),
                     sep=c(3,4), remove=F)
  cCycle <- separate(cCycle, col="SITEID", into=c("TREATMENT", "SITENO"), sep=2, remove=F)
  cCycle$TREATMENT <- as.factor(cCycle$TREATMENT)
  cCycle$SITEID <- as.factor(cCycle$SITEID)
}
