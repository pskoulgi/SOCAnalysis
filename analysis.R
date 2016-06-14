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
    
  
  binBounds = c(0.1, 0.15, 0.25, 0.45, 0.95, 2)
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

# Processing ...
{
  # merge aboveground & belowground data point-wise
  cCycle <- belowGround %>% 
    dplyr::select(-SITE.ID, -POINT.TYPE, -POINT.IN.SITE, -POINT.NO, -SOIL.WT) %>%
    left_join(aboveGround, ., by="POINT")
  
  # Some data exploration visualization for fitting linear models
  sirPred <- dplyr::select(cCycle, TIP, LITTER.CN.RATIO, LITTER.N, LITTER.WT, LITTER.P)
  socPred <- dplyr::select(cCycle, TIP, LITTER.CN.RATIO, LITTER.N, LITTER.WT, LITTER.P,
                    SOIL.SIR, SOIL.N)  
  
  # Variance inflation factors of predictors
  library(usdm)
  vif(sirPred)
  vif(socPred)
  
  # To check for interaction
  coplot(SOIL.C ~ SOIL.SIR | LITTER.P, ylab = "SOC (%)", data=cCycle,
         xlab = "SIR (microgram C-CO2)",
         panel = function(x, y, ...) {
         tmp <- lm(y ~ x, na.action = na.omit)
          abline(tmp)
          points(x, y) })
  
  # Scatter plots of responses vs each predictor
  cCycle <- cCycle %>% 
    separate(col="SITE.ID", into=c("SITE.TYPE", "SITE.NO"), sep=2, remove=F) %>%
    select(-SITE.NO)
  cCycle$SITE.TYPE <- as.factor(cCycle$SITE.TYPE)
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
    geom_point(size=2)
  
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
  ggplot(treeDensMeans, aes(y=MEAN, x=SIZE.CLASS, fill=SITE.TYPE)) +
    geom_bar(position="dodge", stat="identity") +
    geom_errorbar(aes(ymin=MEAN-SE, ymax=MEAN+SE), width=.1, position=position_dodge(0.9)) +
    scale_x_discrete(labels=c(paste(binBounds[1]*100,binBounds[2]*100, sep=" - "),
                              paste(binBounds[2]*100,binBounds[3]*100, sep=" - "),
                              paste(binBounds[3]*100,binBounds[4]*100, sep=" - "),
                              paste(binBounds[4]*100,binBounds[5]*100, sep=" - "),
                              paste(binBounds[5]*100,binBounds[6]*100, sep=" - "))) +
    xlab("DBH (cm) class") + ylab("Stems/ha") + ylim(0,200) +
    theme(legend.position=c(0.9, 0.8), legend.title=element_blank()) +
    scale_fill_hue(labels=c("Contiguous", "Fragment"))
  ggplot(treeStemPropsMeans, aes(y=MEAN, x=SIZE.CLASS, fill=SITE.TYPE)) +
    geom_bar(position="dodge", stat="identity") + 
    geom_errorbar(aes(ymin=MEAN-SE, ymax=MEAN+SE), width=.1, position=position_dodge(0.9)) +
    scale_x_discrete(labels=c(paste(binBounds[1]*100,binBounds[2]*100, sep=" - "),
                              paste(binBounds[2]*100,binBounds[3]*100, sep=" - "),
                              paste(binBounds[3]*100,binBounds[4]*100, sep=" - "),
                              paste(binBounds[4]*100,binBounds[5]*100, sep=" - "),
                              paste(binBounds[5]*100,binBounds[6]*100, sep=" - "))) +
    xlab("DBH (cm) class") + ylab("Proportion of stems") + ylim(0,0.4) +
    theme(legend.position = c(0.9, 0.8), legend.title=element_blank()) +
    scale_fill_hue(labels=c("Contiguous", "Fragment"))
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