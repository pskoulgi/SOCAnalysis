library(dplyr)
library(tidyr)
library(reshape2)

# Read-in soil C, N data
{
  soilCN <- read.table("F:/Workspace/MSc/Project/soc-fragments/Pradeep_CoorgSoilCN.csv",
                       as.is=T, header=T, sep=",")
  
  # Clean-up col names, drop needless cols
  soilCN <- rename(soilCN, POINT=SAMPLE)
  soilCN <- select(soilCN, c(POINT, MASS, C_PERC, N_PERC))
  
  # Remove "Blank", "Soil 502-309" standard sample rows
  soilCN <- filter(soilCN, grepl("CT|FR", POINT))
  
  # Copy CT4_PT0 row to CT5_PT0 row since 
  # they are the same point
  t <- soilCN[grep("CT4_PT0",soilCN$POINT), ]
  t$POINT[1] <- "CT5_PT0"
  soilCN <- bind_rows(soilCN, t)
  
  # Group points according to sites
  soilCN <- arrange(soilCN, POINT)
  
  soilCN <- separate(soilCN, col="POINT", into=c("SITEID", "POINTTYPE", "POINTID"), sep=c(3,4))
  soilCN <- separate(soilCN, col="SITEID", into=c("TREATMENT", "SITENO"), sep=2, remove=F)
  soilCN$TREATMENT <- as.factor(soilCN$TREATMENT)
  soilCN$SITEID <- as.factor(soilCN$SITEID)
}

# Read-in aboveground data
{
  # Litter C, N data
  litterCN <- read.table("F:/Workspace/MSc/Project/soc-fragments/Pradeep_CoorgLitterCN.csv",
                         as.is=T, header=T, sep=",")
  
  litterCN <- rename(litterCN, POINT=SAMPLE)
  litterCN <- select(litterCN, c(POINT, MASS, C_PERC, N_PERC))
  
  # Remove "Blank", "EDTA*" standard sample rows
  litterCN <- filter(litterCN, grepl("CT|FR", POINT))
  
  # Group points according to sites
  litterCN <- arrange(litterCN, POINT)
  
  # Calculate C/N ratio
  litterCN <- mutate(litterCN, CNRATIO=C_PERC/N_PERC)
  
  # Litter P data
  litterP <- read.table("F:/Workspace/MSc/Project/soc-fragments/ICP_calculated_P_pradeep.csv",
                        as.is=T, header=T, sep=",")
  
  litterP <- rename(litterP, POINT=SAMPLE)
  litterP <- select(litterP, c(POINT, P_PERC))
  litterP <- filter(litterP, grepl("CT|FR", POINT))
  
  # Litter weights
  litterWts <- read.table("F:/Workspace/MSc/Project/soc-fragments/litterProperties.csv",
                          as.is=T, header=T, sep=",")
  litterWts$AVGWT <- rowMeans(litterWts[ ,c("DRYWT1", "DRYWT2", "DRYWT3")], na.rm=TRUE)
  litterAvgWt <- select(litterWts, c(POINT, AVGWT))
  
  litter <- left_join(litterCN, litterP, by="POINT")
  litter <- left_join(litter, litterAvgWt, by="POINT")
  litter <- select(litter, -MASS)
  
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
  
  litter <- left_join(litter, tip)
  
  # Compute basal area
  treeDistrib$TreeBasArea <- (pi*(treeDistrib$DBH)^2)/4
  basArea <- treeDistrib %>% group_by(POINT) %>% 
    summarise(BASAREA=sum(TreeBasArea, na.rm=T)*10000/(pi*10^2))
  treeDens <- treeDistrib %>% group_by(POINT) %>% 
    summarise(TREEDENS=n()*10000/(pi*10^2))
  
  agParams <- left_join(litter, basArea)
  agParams <- left_join(agParams, treeDens)
  
  # Read-in soil SIR data
  {
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
    soilSIRC <- soilSIR %>% rename(POINT=SAMPLE) %>%
      select(POINT, C.RELEASED)
    
  }
  
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
