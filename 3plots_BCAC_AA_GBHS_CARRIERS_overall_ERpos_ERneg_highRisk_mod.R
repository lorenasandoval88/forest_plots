require(gridExtra)
library(ggplot2)
library(gridExtra)
#install.packages("dplyr")
library(dplyr)

getwd()
setwd("C:/Users/sandovall2/Box/GBHS/GBHS Collaborator Folder/GBHS Projects/Genetic Susceptibility/GBHS BRIDGES/Lorena Plotting Results/GBHS_AA_BCAC_CARRIERS_high_risk_and_mod/updateddata06172021")
#colors in hexadecimal code http://www.sthda.com/english/wiki/colors-in-r
#shapes in geom_point http://www.sthda.com/english/wiki/be-awesome-in-ggplot2-a-practical-guide-to-be-highly-effective-r-software-and-data-visualization
# OVERALL--------------------------------------------------------------------
data_a <- data.frame(read.csv("input_BCAC_AA_GBHS_CARRIERS_overall_06172021.csv")) 

data_a$group <- factor(data_a$group, levels=c("Carriers", "BCAC - Pathogenic missense variants","BCAC - PTV","African American","GBHS"))#, labels=c("MBB", "MAA", "MCC"))

data1_a<- data_a[1:2]
data2_a<- log(data_a[3:5])
data3_a<-data_a[6:13]
data_a<- cbind(data1_a,data2_a,data3_a)
#ORDER of genes top to bottom
data_a$gene <- factor(data_a$gene[1:4], levels = data_a$gene[1:4])
data3a <- data2_a %>% mutate(LCL_l_1 = ifelse(lowerLimit < -2, oddsRatio  -.3, NA),
                             UCL_l_1 = ifelse(upperLimit ==data2_a[1,3],data2_a[1,3]- oddsRatio, NA))#position of row outside limit

data3_a<- data3a[4:5]
data3_a$LCL_l_1<-as.numeric(data3_a$LCL_l_1)
data_a<- cbind(data_a,data3_a)
# ER positive--------------------------------------------------------------------
data_b <- data.frame(read.csv("input_BCAC_AA_GBHS_CARRIERS_ERpos_06172021.csv")) 
data_b$group <- factor(data_b$group, levels=c("Carriers", "BCAC - Pathogenic missense variants","BCAC - PTV","African American","GBHS"))#, labels=c("MBB", "MAA", "MCC"))

data1_b<- data_b[1:2]
data2_b<- log(data_b[3:5])
data3_b<-data_b[6:13]
data_b<- cbind(data1_b,data2_b,data3_b)
#ORDER of genes top to bottom
data_b$gene <- factor(data_b$gene[1:4], levels = data_b$gene[1:4])
data3b <- data2_b %>% mutate(LCL_l_1 = ifelse(lowerLimit < -2, oddsRatio  -.3, NA),
                             UCL_l_1 = ifelse(upperLimit ==data2_b[1,3],data2_b[1,3]- oddsRatio, NA))#position of row outside limit

data3_b<- data3b[4:5]
data3_b$LCL_l_1<-as.numeric(data3_b$LCL_l_1)
data_b<- cbind(data_b,data3_b)
# ER negative--------------------------------------------------------------------
data_c <- data.frame(read.csv("input_BCAC_AA_GBHS_CARRIERS_ERneg_06172021.csv")) 
data_c$group <- factor(data_a$group, levels=c("Carriers", "BCAC - Pathogenic missense variants","BCAC - PTV","African American","GBHS"))#, labels=c("MBB", "MAA", "MCC"))

data1_c<- data_c[1:2]
data2_c<- log(data_c[3:5])
data3_c<-data_c[6:13]
data_c<- cbind(data1_c,data2_c,data3_c)
#ORDER of genes top to bottom
data_c$gene <- factor(data_c$gene[1:4], levels = data_c$gene[1:4])

# add arrow segments, for lower bound, we want to cut off at 0.3, replace lowest GBHS susceptibility OR for -1.203973
# log(0.3) = -1.203973
data2_c[8,2] <- -1.203973
data_c[8,4] <- -1.203973
data3c <- data2_c %>% mutate(LCL_l_1 = ifelse(lowerLimit ==data2_c[8,2], data2_c[8,2] - oddsRatio, NA),
                             UCL_l_1 = ifelse(upperLimit ==data2_c[1,3],data2_c[1,3]- oddsRatio, NA))#position of row outside limit
data3c
data3_c<- data3c[4:5]
data3_c$LCL_l_1<-as.numeric(data3_c$LCL_l_1)
data_c<- cbind(data_c,data3_c)

########## Function to limit to one decimal point, nolonger need, see line 86
#scaleFUN1 <- function(x) sprintf("%.1f", x)
#scaleFUN0 <- function(x) sprintf("%.0f", x)

# OVERALL -------------------------------------------------
#ATM, BARD1,CHECK2,RAD51C, RAD51D, TP53 IN BCAC
a = ggplot(data=data_a,
           aes(x = group,y = oddsRatio, ymin = lowerLimit, ymax = upperLimit ))+
  geom_pointrange(aes(col=group),stroke=0)+
  geom_hline(aes(fill=group),yintercept =0, linetype=2)+
  xlab('Gene')+ ylab("Overall Odds Ratio (95% Cl)")+
  geom_errorbar(aes(ymin=lowerLimit, ymax=upperLimit,col=group),width=0,cex=0.6)+
  facet_wrap(~gene,strip.position="left",nrow=15,scales = "free_y") +
  coord_trans(y = scales::exp_trans()) +
  
  ########### x axis ticks and labels
  scale_y_continuous(breaks = (c(-1.203973,0,1.609438,3.218876, 4.60517, 5.521461)), # removed 5.010635 (150 tick)                      
                     minor_breaks = NULL,
                     labels = c(0.3, 1, 5, 25, 100, 250))+
  
  expand_limits(y=c(-1.203973,5.521461))+
  
  theme(plot.title=element_text(size=9,face="bold"),
        plot.margin = unit(c(1.2,0.5,0.1,0,0.08), "cm"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=12,face="bold"),
        legend.position = c(4.5,0.85), #USE WHEN THERE ARE 3 PLOTS SIDE-BY-SIDE
        # legend.position = "none",
        # legend.title = element_blank(),
        
        # y axis gene labels listed horizontally,  strip.text.y.left
        strip.text.y= element_text(hjust=0,vjust = 1,angle=0,size =11))+
  
  geom_segment(aes(x = group, xend = group, y = oddsRatio, yend = oddsRatio + UCL_l_1),
               colour="#CC99FF", arrow= arrow(length = unit(0.3, "cm")),size=0.6)+
  scale_color_manual(values = c("GBHS"="#FF9999", "African American"="#CC99FF","BCAC - PTV"="#0066FF","BCAC - Pathogenic missense variants"="#33CC66", "Carriers"="turquoise2"))+
  #                            GBHS RED, AA PURPLE, NIGERIA DARK BLUE
  # add OR bars and points
  geom_point(aes(shape=as.factor(gear)),colour="#CC99FF",size=0.3)+
  geom_point(aes(shape=as.factor(gear2)),colour="#FF9999",size=0.3)+
  geom_point(aes(shape=as.factor(gear3)),colour="#0066FF",size=0.3)+
  geom_point(aes(shape=as.factor(gear4)),colour="#33CC66",size=0.3)+
  geom_point(aes(shape=as.factor(gear5)),colour="turquoise2",size=0.3)+
  
  coord_flip()
a

# ER positive -------------------------------------------------
#ATM, BARD1,CHECK2,RAD51C, RAD51D, TP53 IN BCAC
b = ggplot(data=data_b,
           aes(x = group,y = oddsRatio, ymin = lowerLimit, ymax = upperLimit ))+
  geom_pointrange(aes(col=group),stroke=0)+
  
  geom_hline(aes(fill=group),yintercept =0, linetype=2)+
  labs(y="ER+ Odds Ratio (95% Cl)",x="")+
  geom_errorbar(aes(ymin=lowerLimit, ymax=upperLimit,col=group),width=0,cex=0.6)+
  facet_wrap(~gene,strip.position="left",nrow=15,scales = "free_y") +
  coord_trans(y = scales::exp_trans()) +
  
  ########### x axis ticks and labels
  scale_y_continuous(breaks = (c(-1.203973,0,1.609438,3.218876, 4.60517, 5.521461)), # removed 5.010635 (150 tick)                      
                     minor_breaks = NULL,
                     labels = c(0.3, 1, 5, 25, 100, 250))+
  
  expand_limits(y=c(-1.203973,5.521461))+
  
  theme(plot.title=element_text(size=9,face="bold"),
        plot.margin = unit(c(1.2,0.5,0.1,0,0.08), "cm"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=12,face="bold"),
        legend.position = "none",
        #legend.text=element_blank(),
        
        # y axis gene labels listed horizontally,  strip.text.y.left
        strip.text.y= element_text(hjust=0,vjust = 1,angle=0,size =11))+
  
  geom_segment(aes(x = group, xend = group, y = oddsRatio, yend = oddsRatio + UCL_l_1),
               colour="#CC99FF", arrow= arrow(length = unit(0.3, "cm")),size=0.6)+
  scale_color_manual(values = c("GBHS"="#FF9999", "African American"="#CC99FF","BCAC - PTV"="#0066FF","BCAC - Pathogenic missense variants"="#33CC66","Carriers"="turquoise2"))+
  #                            GBHS RED, AA PURPLE, NIGERIA DARK BLUE
  
  geom_point(aes(shape=as.factor(gear)),colour="#CC99FF",size=0.3)+
  geom_point(aes(shape=as.factor(gear2)),colour="#FF9999",size=0.3)+
  geom_point(aes(shape=as.factor(gear3)),colour="#0066FF",size=0.3)+
  geom_point(aes(shape=as.factor(gear4)),colour="#33CC66",size=0.3)+
  geom_point(aes(shape=as.factor(gear5)),colour="turquoise2",size=0.3)+
  
  coord_flip()
b
# ER negative -------------------------------------------------
#ATM, BARD1,CHECK2,RAD51C, RAD51D, TP53 IN BCAC 
c = ggplot(data=data_c,
           aes(x = group,y = oddsRatio, ymin = lowerLimit, ymax = upperLimit ))+
  geom_pointrange(aes(col=group),stroke=0)+
  
  geom_hline(aes(fill=group),yintercept =0, linetype=2)+
  labs(y="ER- Odds Ratio (95% Cl)",x="")+
  geom_errorbar(aes(ymin=lowerLimit, ymax=upperLimit,col=group),width=0,cex=0.6)+ 
  facet_wrap(~gene,strip.position="left",nrow=15,scales = "free_y") +
  coord_trans(y = scales::exp_trans()) +
  
  ########### x axis ticks and labels
  scale_y_continuous(breaks = (c(-1.203973,0,1.609438,3.218876, 4.60517, 5.521461)), # removed 5.010635 (150 tick)                      
                     minor_breaks = NULL,
                     labels = c(0.3, 1, 5, 25, 100, 250))+
  
  expand_limits(y=c(-1.203973,5.521461))+
  
  theme(plot.title=element_text(size=9,face="bold"),
        plot.margin = unit(c(1.2,0.5,0.1,0,0.08), "cm"),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_text(size=12),
        axis.title=element_text(size=12,face="bold"),
        legend.position = "none",
        #legend.text=element_blank(),
        
        # y axis gene labels listed horizontally,  strip.text.y.left
        strip.text.y= element_text(hjust=0,vjust = 1,angle=0,size =11))+
  
  geom_segment(aes(x = group, xend = group, y = oddsRatio, yend = oddsRatio + UCL_l_1), 
               colour="#CC99FF", arrow= arrow(length = unit(0.3, "cm")),size=0.6)+
  geom_segment(aes(x = group, xend = group, y = oddsRatio, yend = oddsRatio + LCL_l_1), 
               colour="#FF9999", arrow= arrow(length = unit(0.3, "cm")),size=0.6)+
  scale_color_manual(values = c("GBHS"="#FF9999", "African American"="#CC99FF","BCAC - PTV"="#0066FF","BCAC - Pathogenic missense variants"="#33CC66","Carriers"="turquoise2"))+
  #                            GBHS RED, AA PURPLE, NIGERIA DARK BLUE
  
  geom_point(aes(shape=as.factor(gear)),colour="#CC99FF",size=0.3)+
  geom_point(aes(shape=as.factor(gear2)),colour="#FF9999",size=0.3)+
  geom_point(aes(shape=as.factor(gear3)),colour="#0066FF",size=0.3)+
  geom_point(aes(shape=as.factor(gear4)),colour="#33CC66",size=0.3)+
  geom_point(aes(shape=as.factor(gear5)),colour="turquoise2",size=0.3)+
  
  coord_flip()
#c= c +  geom_point(aes(gene,oddsRatio, shape=as.factor(gear)))
c
grid.arrange(a,b,c,ncol=5)
