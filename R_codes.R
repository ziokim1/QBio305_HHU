### Initialization ----
# load packages
library(ggplot2)    # load the ggplot2 library for plotting 
library(car)        # load the car library to make a QQplot
library(qqman)      # load the qqman library to make manhattan plots
library(cowplot)    # plot multigraphs
library(qtl)        # load the qtl library for linkage analysis
options(scipen=999) # remove scientific notation of very low integers   


# load data
setwd("~/Documents/QBIO305_group1_assignment_data") # set the directory to the measurement data

dta.raw <- read.csv("305_ANOVA_QTL_analysis_data.csv")
map <- read.csv("305_owb_input_map.csv")

### ANOVA assumption test ----
dta.singleM <- dta.raw[, c(1,3,155,500,1000)] # select random genetic markers at each chromosome to check for ANOVA assumptions
colnames(dta.singleM) <- c("group", "trait", "allele_CH1", "allele_CH2", "allele_CH3")

# conduct the raw ANOVA by inserting the correct model
model.CH1 <- aov(trait ~ allele_CH1, dta.singleM)
model.CH2 <- aov(trait ~ allele_CH2, dta.singleM)
model.CH3 <- aov(trait ~ allele_CH3, dta.singleM)

# extract the residuals
resid.CH1 <- model.CH1$residuals
resid.CH2 <- model.CH2$residuals
resid.CH3 <- model.CH3$residuals

# plot the histogram
his.CH1 <- ggplot(data = data.frame(1:length(resid.CH1), resid.CH1), 
       aes(x = resid.CH1)) +
  geom_histogram(aes(y =..density..), color="darkgrey", bins=50) +
  stat_function(fun = dnorm, 
                args = list(mean = mean(resid.CH1), 
                            sd = sd(resid.CH1)),
                color="darkred")

his.CH2 <- ggplot(data = data.frame(1:length(resid.CH2), resid.CH2), 
                  aes(x = resid.CH2)) +
  geom_histogram(aes(y =..density..), color="darkgrey", bins=50) +
  stat_function(fun = dnorm, 
                args = list(mean = mean(resid.CH2), 
                            sd = sd(resid.CH2)),
                color="darkred")

his.CH3 <- ggplot(data = data.frame(1:length(resid.CH3), resid.CH3), 
                  aes(x = resid.CH3)) +
  geom_histogram(aes(y =..density..), color="darkgrey", bins=50) +
  stat_function(fun = dnorm, 
                args = list(mean = mean(resid.CH3), 
                            sd = sd(resid.CH3)),
                color="darkred")

his <- plot_grid(his.CH1, his.CH2, his.CH3, labels = c('Marker_Chromosome 1', 'Marker_Chromosome 2', 'Marker_Chromosome 3'), label_size = 18, nrow = 3)
title <- ggdraw() + draw_label("Residual distributions of genetic markers", fontface='bold')
plot_grid(title, his, ncol=1, rel_heights=c(0.1, 1))

# plot qq-plot
qq <- par(mfrow=c(1,3), pty="s", cex.lab=1.3)
q1 <- qqPlot(resid.CH1, envelope = F)
q2 <- qqPlot(resid.CH2, envelope = F)
q3 <- qqPlot(resid.CH3, envelope = F)
mtext("QQ-Plots", side = 3, line = -15, outer = TRUE, cex=2)
on.exit(par(qq))

# check for homoscedasticity and outliers 
rf <- par(mfrow=c(1,3), pty="s", cex.lab=1.3)
plot(model.CH1, which = 1)
plot(model.CH2, which = 1)
plot(model.CH3, which = 1)
mtext("Residuals vs Fitted values", side = 3, line = -15, outer = TRUE, cex=2)
on.exit(par(rf))


### QTL Analysis (Plant Height) ----

# loop through all markers and conduct ANOVAs for each one of them
lis.out <- list()

for(m in 5:ncol(dta.raw)){
  # extract the first two columns (group + spike length) and the m-th column (=marker)
  loop.dta <- dta.raw[,c(2,3,m)]
  # save the name of the current marker
  loop.marker.name <- colnames(loop.dta)[3]
  # in the model in the aov command, we will have to define the independent 
  # variable name. since the marker name changes in each loop, we need to 
  # to change the column names here to have the same marker name in each loop
  colnames(loop.dta) <- c("group","trait", "allele")
  # conduct one-way ANOVA
  loop.aov <- aov(trait ~ allele + group, loop.dta)
  # extract the allele's p-value
  loop.pval.allele <- summary(loop.aov)[[1]][1,5]
  loop.pval.group <- summary(loop.aov)[[1]][2,5]
  # extract the marker's genetic position from the genetic map
  loop.map <- map[which(map$marker == loop.marker.name),]
  # create the output data frame
  loop.df.out <- data.frame(loop.map,
                            pval.allele=loop.pval.allele,
                            pval.group=loop.pval.group)
  # save the output data frame in the list
  lis.out[[m]] <- loop.df.out
}
# combine the loop's output data frames into a single data frame
res.aov <- do.call("rbind", lis.out)

# create a data frame that follows the requirements of the 
# manhattan command of the qqman library (run as is)
dta.plot <- data.frame(CHR = as.numeric(gsub("H","",res.aov$chr)),
                       BP = res.aov$pos,
                       SNP = 1:nrow(res.aov),
                       pval.marker = res.aov$pval.allele,
                       pval.group = res.aov$pval.group)

# calculate the negative logarithm of the Bonferroni corrected significance threshold
sig.threshold.BonfCorrected <- -log(0.05/nrow(dta.plot))

# plot genome-wide p-values for markers
manhattan(dta.plot, genomewideline = sig.threshold.BonfCorrected,
          suggestiveline = F, logp=T, p="pval.marker", type="l", 
          lwd=3, ylab="-log10(p-values)", main="Multiple marker QTL analysis")
text(80,10,  "a_(Bonferroni)", cex=1, pos=3,col="red")

### Linkage mapping (Seed Colour) ----

#read in linkage map and qualitative phenotypes as denoted in marker annotation
owb <- read.cross("csv", ".", "305_linkage_map_data.csv", genotypes=c("a","b"), 
                  alleles=c("a", "b"), crosstype="dh")
#have a look at the complete linkage map
plotMap(owb, main="", show.marker.names=T)

#To map the locus, calculate all pairwise recombination frequencies and LOD scores
seed_colour <- tryallpositions(owb, "seed_colour", error.prob=0)
#show the best linkage to a marker from each chromosome
summary(seed_colour)
#move the locus to the best marker position
owb <- movemarker(owb, "seed_colour", "1H", 147.12)
#update map
plotMap(owb, main="", show.marker.names=T)

# writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
