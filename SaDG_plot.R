#--- Authorship
## AUTHOR : Perales-Paton, Javier
## E-MAIL : jperales@cnio.es
## DATE : Oct2015
#--- About this
## TITLE : Saturation curve of detected genes by down-sampling BAM files
## DESCRIPTION : 
#   Technical QC. comprehensive saturation curve of detected gene expression by increasing read depth.
## DEPENDENCIES :
#   - ggplot2
#   - plyr
## It is a free piece of code. Please, let me know improvements and ideas.



require("ggplot2")
require("plyr")

## Summarizes data.
# Source: http://www.cookbook-r.com/Graphs/Plotting_means_and_error_bars_(ggplot2)/
## Gives count, mean, standard deviation, standard error of the mean, and confidence interval (default 95%).
##   data: a data frame.
##   measurevar: the name of a column that contains the variable to be summariezed
##   groupvars: a vector containing names of columns that contain grouping variables
##   na.rm: a boolean that indicates whether to ignore NA's
##   conf.interval: the percent range of the confidence interval (default is 95%)
summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(N    = length2(xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm),
                     sd   = sd     (xx[[col]], na.rm=na.rm)
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  datac$ci <- datac$se * ciMult
  
  return(datac)
}


# Simulated data
DF <- data.frame(sname=c(rep("S1",10),rep("S2",10)),
                 type=rep(seq(from=0.5,to=2.5,by=0.5),20),
                 Ngenes=rep(log2(seq(from=1,to=5000,by=1000)),20)-runif(min = 1,max=2.5,n=20))
fl <- "/local/jperales/DimitrovS_singleCell/Results/Saturation_curve/C1_saDG_table_v4.tsv";
fl <- "/local/jperales/MullallyA/AM_results/QC/Saturationcurve.tsv"
DF <- read.table(fl,sep = "\t",col.names=c("sname","type","NgenesByCPM","Ngenes"))

DF$type <- DF$type*1e-6

DFc <- summarySE(DF, measurevar="Ngenes",groupvars=c("sname","type"),na.rm = TRUE)


ggplot(DFc, aes(x=type, y=Ngenes, colour=sname, group=sname)) + 
  geom_errorbar(aes(ymin=Ngenes-se, ymax=Ngenes+se), colour="black", width=.1, position=pd) +
  geom_line(position=pd) +
  geom_point(position=pd, size=2, shape=21, fill="white") + # 21 is filled circle
  xlab("Sampled reads (Million)") +
  ylab("No. Detected Genes") +
  ylim(c(0,15000)) +
  scale_x_continuous(limits=c(0, 20),breaks=seq(from=0,to=20,by=1)) +
  scale_colour_hue(name="Sample Name",    # Legend label, use darker colors
                   breaks=unique(DFc$sname),
                   labels=unique(DFc$sname),
                   l=40) +                    # Use darker colors, lightness=40
  expand_limits(y=0) +                        # Expand y range
  theme_bw() + theme(axis.text.y= element_text(size=15),
                     axis.text.x= element_text(size=12))

