#alpha diversity from qiime raw data output

rm(list=ls())
library(ggplot2)

setwd("YOUR_WORKING_DIRECTORY")


###########Figure 2############
alpha_raw <- read.csv("YOUR_WORKING_DIRECTORY/Alpha_raw.csv")
alpha_raw$Compartment <- factor(x = alpha_raw$Compartment,levels = c("Soil","Rhizosphere","Roots","Shoots"))

my_breaks <- function(x) { if (max(x) <12 ) seq(7, 9,11) else seq(0, 1, 2) }

p <- ggboxplot(alpha_raw,x="Antibiotics", y="shannon",color="Compartment",palette = "jco")+
  facet_wrap(~Compartment,scales="free_y",nrow=1)+
  stat_compare_means(fontface="italic",label.x =1 )+
  labs(x="Antibiotics in Irrigation Water", y=expression("Shannon Diversity")) + 
  ggtitle("Alpha Diversity Statistics")+
  theme(legend.position="top",axis.text=element_text(size=16),axis.title=element_text(size=16),plot.title = element_text(size=22, hjust=0.5),text = element_text(size = 16),axis.text.x = element_text(angle = 0))
#scale_y_continuous(breaks = scales::pretty_breaks(4), limits = c(7, 10.5)) 
#The line above is for axis adjustment
dat_text <- data.frame(
  label = c("ASD=30,000","ASD=30,000", "ASD=10,000", "ASD=500"),
  Compartment = c("Soil","Rhizosphere","Roots","Shoots")
)

ALPHAFIG <- p + geom_text(
  data    = dat_text,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -0.1,
  vjust   = -1
)

#https://community.rstudio.com/t/separate-axis-breaks-for-facets/10352
#https://fishandwhistle.net/post/2018/modifying-facet-scales-in-ggplot2/
#Figure panels mannually adjusted in publication for axis consistency 


#########Figure S2############
rm(list=ls())#Clean work space
setwd("YOUR_WORKING_DIRECTORY")
alpha_raw <- read.csv("YOUR_WORKING_DIRECTORY/Alpha_raw.csv")

alpha_raw$Time <- gsub("Day1","Day25",alpha_raw$Time)
alpha_raw$Time <- gsub("\\<Day2\\>","Day27",alpha_raw$Time)
alpha_raw$Time <- gsub("\\<Day3\\>","Day30",alpha_raw$Time)
alpha_raw$Time <- gsub("\\<Day4\\>","Day35",alpha_raw$Time)

alpha_raw$Compartment <- factor(x = alpha_raw$Compartment,levels = c("Soil","Rhizosphere","Roots","Shoots"))

p <- ggboxplot(alpha_raw,x="Time", y="shannon",color="Compartment",palette = "jco")+facet_wrap(~Compartment,scales="free_y",nrow=1)+stat_compare_means(fontface="italic",label.x =1.2 )+labs(x="Harvest Day", y=expression("Shannon Diversity")) + ggtitle("Alpha Diversity Statistics")+theme(legend.position="top",axis.text=element_text(size=16),axis.title=element_text(size=16),plot.title = element_text(size=22, hjust=0.5),text = element_text(size = 16),axis.text.x = element_text(angle = 30))


dat_text <- data.frame(
  label = c("ASD=30,000", "ASD=30,000","ASD=10,000", "ASD=500"),
  Compartment = c("Soil","Rhizosphere","Roots","Shoots")
)

p + geom_text(
  data    = dat_text,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -0.1,
  vjust   = -1
)

```

#Days divide antibiotics yes or no
```{r}
alpha_raw_yes <- alpha_raw[ which(alpha_raw$Antibiotics=='Yes'), ]

alpha_raw_yes$Compartment <- factor(x = alpha_raw_yes$Compartment,levels = c("Soil","Rhizosphere","Roots","Shoots"))

p2 <- ggboxplot(alpha_raw_yes,x="Time", y="shannon",color="Compartment",palette = "jco")+facet_wrap(~Compartment,scales="free_y",nrow=1)+stat_compare_means(fontface="italic",label.x =1.2 )+labs(x="Harvest Day", y=expression("Shannon Diversity")) + ggtitle("Alpha Diversity Statistics_Yes")+theme(legend.position="top",axis.text=element_text(size=16),axis.title=element_text(size=16),plot.title = element_text(size=22, hjust=0.5),text = element_text(size = 16),axis.text.x = element_text(angle = 30))+geom_text(
  data    = dat_text,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -0.1,
  vjust   = -1)+scale_y_continuous(breaks = scales::pretty_breaks(4), limits = c(0, 2)) 


dat_text <- data.frame(
  label = c("ASD=30,000", "ASD=30,000","ASD=10,000", "ASD=500"),
  Compartment = c("Soil","Rhizosphere","Roots","Shoots"))
dat_text$Compartment <- factor(x = dat_text$Compartment,levels = c("Soil","Rhizosphere","Roots","Shoots"))



alpha_raw_no <- alpha_raw[ which(alpha_raw$Antibiotics=='No'), ]
p3 <- ggboxplot(alpha_raw_no,x="Time", y="shannon",color="Compartment",palette = "jco")+facet_wrap(~Compartment,scales="free_y",nrow=1)+stat_compare_means(fontface="italic",label.x =1.2 )+labs(x="Harvest Day", y=expression("Shannon Diversity")) + ggtitle("Alpha Diversity Statistics_No")+theme(legend.position="top",axis.text=element_text(size=16),axis.title=element_text(size=16),plot.title = element_text(size=22, hjust=0.5),text = element_text(size = 16),axis.text.x = element_text(angle = 30))+geom_text(
  data    = dat_text,
  mapping = aes(x = -Inf, y = -Inf, label = label),
  hjust   = -0.1,
  vjust   = -1)
#+scale_y_continuous(breaks = scales::pretty_breaks(4), limits = c(7, 10.5)) 
#Figure panels mannually adjusted in publication for axis consistency 

