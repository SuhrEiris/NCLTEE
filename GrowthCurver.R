## Growth curve of wildtypes ##

#### Packages ####

library("dplyr")
library("tidyr")
library("ggplot2")
library("devtools")
library("BiocManager")
library("readr")
library("wesanderson")
library("ggpubr")
#remotes::install_github("Russel88/growthcurver", force = TRUE)
library('growthcurver')
library("readxl")
library("stringr")
library('gridExtra')

#### Data ####

dir()

#Import data
biorep1 <- as.data.frame(read_excel("GK_f2_BioRep1_10032021.xlsx", sheet = 3))


# Plot growth curves
# Make data longformat and discard blanks
biorep1_longf <- biorep1[,1:70] %>% 
  gather(Strain, OD, 2:70) %>%
  separate(Strain, into=c("Strain", "dilution", "TechRep"), sep= "_")

#Overview
ggplot(biorep1_longf, aes(x=Time, y = OD, color = Strain))+
  geom_point()+
  facet_grid(Strain~dilution)

#Only WT and D in 100x
biorep1_longf_reduced <- biorep1_longf %>%
  filter(!grepl("Wino", Strain)) %>% 
  filter(!grepl("1000x", dilution))

#Overview
ggplot(biorep1_longf_reduced, aes(x=Time, y = OD, color = TechRep))+
  geom_point()+
  facet_grid(.~Strain)

#summarise 
sum_reduced = biorep1_longf_reduced %>%
  group_by(Strain, Time) %>%
  summarise( mean.OD = mean(OD),
             sd.OD = sd(OD),
             count.tech = n())

#Plot
curves = ggplot(sum_reduced, aes(x=Time, y = mean.OD, color = Strain))+
  geom_line(size = 1)+
  geom_ribbon(aes(ymax = mean.OD + sd.OD, ymin = mean.OD - sd.OD, fill = Strain), alpha = 0.2, color = NA) +
  #facet_grid(.~Strain) + 
  scale_fill_manual(values=c("#046C9A", "#D69C4E"))+
  scale_color_manual(values=c("#046C9A", "#D69C4E"))+
  theme_bw()+
  xlab("\n Time (min)")+
  ylab(" OD (590nm)")+
  theme(axis.text.x = element_text(size=18))+
  theme(axis.title = element_text(size = 14))+
  theme(axis.text = element_text(size = 14, color = "Black")) + 
  ggtitle("Growth curves")


#Growth curver
biorep1_gc <- SummarizeGrowthByPlate(biorep1, plot_fit = TRUE, plot_file = "Biorep1_raw.pdf")

#With background correction
biorep1_bg <- biorep1[,1:71]
biorep1_bgcorrect <- SummarizeGrowthByPlate(biorep1_bg, bg_correct = "blank", plot_fit = TRUE, plot_file = "Biorep1.pdf")

#Add biorep column
# biorep1_bgcorrect = biorep1_bgcorrect %>% 
#   mutate(BioRep = 1)


# data wrangling (get correct isolate names and techrep in sep column)
data <- biorep1_bgcorrect %>% separate(sample, into= c("Strain", "dilution", "TechRep"), sep = "_")
data1 <- data[1:69,] #remove NA
data2 <- data1[-58,] # remove questionale fit

# Save data as csv file
write.csv(data2, "GrowthCurve_SummaryData_10032021.csv")



#  sd of all togther
data.gr <- data2 %>%
  group_by(Strain, dilution) %>%
  summarise(count = n(),
            gen.time = mean(t_gen),
            sd.gen = sd(t_gen),
            carry.c = mean(k),
            sd.carry = sd(k),
            growth.rate = mean(r),
            sd.rate = sd(r),
            start.od = mean(n0),
            sd.start = sd(n0))

# data for only wildtype and D 

data.gr1 = data.gr %>%
  filter(!grepl("Wino", Strain)) %>% 
  filter(!grepl("1000x", dilution))


#### Plot ####

# Generation time
gen = ggplot(data.gr1, aes(x= Strain, y= gen.time, color = Strain ))+
  geom_point(size = 5)+
  geom_errorbar(aes(ymin=gen.time + sd.gen, ymax=gen.time-sd.gen), width=0.2, size=1.5) + 
  theme_bw() +
  ylab ("min\n")+
  xlab("\n")+
  ggtitle("Generation time")+
  scale_color_manual(values=c("#046C9A", "#D69C4E"))+
  #facet_wrap(dilution~., scales = "free_x") +
  theme(axis.text.x = element_text(size=18))+
  theme(axis.title = element_text(size = 14))+
  theme(axis.text = element_text(size = 14, color = "Black"))

# CarryC
carry = ggplot(data.gr1, aes(x= Strain, y= carry.c, color = Strain))+
  geom_point(size = 4)+
  geom_errorbar(aes(ymin=carry.c + sd.carry, ymax=carry.c-sd.carry), width=0.2, size=1.5) + 
  theme_bw() +
  ylab ("\n")+
  xlab("")+
  ggtitle("Carrying capacity") +
  scale_color_manual(values=c("#046C9A", "#D69C4E"))+
  theme(axis.text.x = element_text(size=18))+
  theme(axis.title = element_text(size = 14))+
  theme(axis.text = element_text(size = 14, color = "Black"))


# Starting OD estimated with the log func
# start = ggplot(data.gr1, aes(x= Strain, y= start.od, color = Strain))+
#   geom_point(size = 5)+
#   geom_errorbar(aes(ymin= start.od + sd.start, ymax= start.od-sd.start), width=0.2, size=1.5) + 
#   theme_bw() +
#   ylab ("OD590nm\n")+
#   xlab("\n")+
#   ggtitle("Starting OD")+
#   scale_color_manual(values=c("#046C9A", "#D69C4E"))+
#   theme(axis.title = element_text(size = 24), axis.text = element_text(size = 13, color = "Black"))

# Get correct start OD
# correct starting point (need to remoce D outliers)
time0 <- subset(biorep1_longf_reduced, Time =="0")

sum_start = time0 %>%
  group_by(Strain) %>%
  summarise(mean.start = mean(OD),
            sd.start = sd(OD))

start = ggplot(sum_start, aes(x= Strain, y= mean.start, color = Strain))+
  geom_point(size = 5)+
  geom_errorbar(aes(ymin= mean.start + sd.start, ymax= mean.start-sd.start), width=0.2, size=1.5) + 
  theme_bw() +
  ylab ("OD (590nm)\n")+
  xlab("\n")+
  ggtitle("Starting OD")+
  scale_color_manual(values=c("#046C9A", "#D69C4E"))+
  theme(axis.text.x = element_text(size=18))+
  theme(axis.title = element_text(size = 14))+
  theme(axis.text = element_text(size = 14, color = "Black"))


#grid_arrange
grid.arrange(start, carry, gen, curves, ncol = 2)

# T-test 

#devtools::install_github("kassambara/ggpubr")
library("ggpubr") 

# Test day for only 100x

ttest.data = data2 %>%
  filter(!grepl("Wino", Strain)) %>% 
  filter(!grepl("1000x", dilution))

# Assump 1) are the samples independent: yes, they are from different cultures
# Assump 2) are the data from each group normal distributed = do Shapiro-Wilk normality test

with(ttest.data, shapiro.test(t_gen[Strain == "WT"])) # P > 0.0.5 = normal dist, use parametric
with(ttest.data, shapiro.test(t_gen[Strain == "D"])) 

# Assump 3) do the the two pop have the same variances? Use F-test for homogeneity

ftest <- var.test(OD ~ Strain, data = time0)
ftest 
#p > 0.05, no sig difference between variances of the two groups, 
#we can use t.test, assus equal variance

# T.test 
res <- t.test(OD ~ Strain, data = time0, var.equal= FALSE)
res$p.value

#start, p = 0.5156169
#gen, p = 0.0001253141
#carry, k = 8.641081e-08 


