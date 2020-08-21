library(ggplot2)
library(reshape2)
library(jcolors)
library(ggpubr)


setwd("E:/figures/")

#### overall suitability line plot
suit <- read.csv(file="overall suit.csv",header=T,sep=",")


p<- ggplot(suit, aes(x=year, y=suitable, group=RCP, color=RCP)) + 
  geom_line(size=1.5) #+
  #geom_point() 


p+labs(title="Projected Percentage of Suitable Sites")+
  theme_minimal() +
  scale_color_jcolors(palette = "default") +
  scale_y_continuous(labels = scales::percent) +
  theme(legend.title = element_blank()) +
  xlab("year")+ylab("% of reef sites")





#### suitability each variable
suitvar <- read.csv(file="suitvars.csv",header=T, sep=",")

lp <- ggplot(suitvar, aes(x=year, y=suitable, group=variable))+
  geom_line(aes(color=variable), size=2)

lp + facet_wrap(. ~ RCP, scales="free_x",nrow=1)+
  theme_light() +
  scale_color_brewer(palette='Set2')+
  scale_y_continuous(labels = scales::percent) +
  xlab("year")+ylab("% of reef sites")

  
###
# historic = read.csv(file="historic.csv",header=T, sep=",")
#   
# hist <- ggline(historic, x='year', y='suitable',group = 'variable')
# hist



