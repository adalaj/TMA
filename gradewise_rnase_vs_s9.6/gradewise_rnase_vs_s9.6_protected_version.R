# ------------------------------------------------------------------------------
# This code is part of paper: RNA:DNA hybrids determine tissue pathological heterogeneity and aggressiveness in clinically diagnosed cancers. 
# Author: Jyoti Devendra Adala under supervision of Dr. Vladimir A Kuznetsov 
# For updates and contributions, visit : https://github.com/adalaj

# Purpose: This script analyzes the distribution of signal intensities (DAB: Nucleus: Mean) for two experimental conditions—S9.6 and 
#RNase III+H—across various grades (Normal, Grade 1, Grade 2, and Grade 3) in ICNST
#The outputs include density plots, statistical summaries, and files for further analysis.
# ------------------------------------------------------------------------------



#JA: load required libraries
library(data.table)
library(tidyverse)


#JA: Step1: Load dataset and preprocessing

rnase <- fread("BR720_IDC_only_RNase_treated_osk_No_negative.csv", sep = ",", header = TRUE)
#JA:169708


unique(rnase$Grade)
#JA:[1] "1" "1--2" 2"    "3" "2--3" "-"

#JA:select desired columns
rnase<- rnase %>% select("Image", "Name", "Class", "TMA core", "Hematoxylin: Nucleus: Mean", 
                         "DAB: Nucleus: Mean" ,"Pathology diagnosis", "Grade")

#JA:removing all negative values from fig8a input
rnase<- rnase %>% filter(`DAB: Nucleus: Mean` > 0) #JA: there is no negative values in this dataset



#JA:Grade-wise Splitting
grade1 <- rnase %>% filter(Grade == "1")
#JA:1142

grade2<- rnase %>% filter(Grade == "2")
#JA:97891

grade3<- rnase %>% filter(Grade == "3")
#JA:68541


#JA:inset normal file for RNase treated:
normal<- fread("BR720_RNase treated_N_No Negative.csv", header = TRUE, sep = ",")
range(normal$`DAB: Nucleus: Mean`) #JA: no negative values
#JA:[1] 0.0342 0.9654

#JA:select desired columns
normal <- normal %>% select("Image", "Name", "Class", "TMA core", "Hematoxylin: Nucleus: Mean", 
                            "DAB: Nucleus: Mean" ,"Pathology diagnosis", "Grade")
#JA:260

#JA:Setting a column as identifier
grade1$Identifier <- "Grade 1"
grade2$Identifier <- "Grade 2"
grade3$Identifier <- "Grade 3"
normal$Identifier <- "Normal"


#JA: Step2: Data Transformation

#JA:Combing all three dataset as one for analyzing them together
fig8a_input<- rbind(grade1, grade2, grade3, normal)
#JA:1142+97891+68541+260 = 167834

#JA:First find the range for binning
range<- range(fig8a_input$`DAB: Nucleus: Mean`)
#JA:[1] 0.0248  1.1394

a<- range[[1]][1]
b<- range[[2]][1]
bin_breaks_75 <- seq(a, b, length.out=76)

#JA:add bin range in a new column

fig8a_input <- fig8a_input %>%
  dplyr::mutate(Bin = cut(`DAB: Nucleus: Mean`, breaks = bin_breaks_75, include.lowest = TRUE))
fig8a_input$Identifier <- factor(fig8a_input$Identifier, levels = c("Normal", "Grade 1", "Grade 2", "Grade 3"))

fwrite(fig8a_input, "BR720_fig8a_grade1_vs_2_vs_3_vs_normal_RNase_density_input.csv")

#JA: Calculate the mean of nucleus signal intensities in that range not necessarily same as mean of midpoint 
fig8a_graph_input <- fig8a_input %>%
  group_by(Bin, Identifier) %>%
  summarise(
    Mean_Nucleus_signal_intensity = round(mean(`DAB: Nucleus: Mean`, na.rm = TRUE),2), #JA:mean of nucleus intensity calculation within that unique group
    Frequency = n()
  )



#JA:Displaying bin midpoints
fig8a_graph_input$Bin_char <- gsub("\\[|\\]|\\(|\\)", "", fig8a_graph_input$Bin)
fig8a_graph_input <- fig8a_graph_input %>%
  separate(Bin_char, c('lower_limit', 'upper_limit'), sep = ",")
fig8a_graph_input$lower_limit<- as.numeric(fig8a_graph_input$lower_limit)
fig8a_graph_input$upper_limit<- as.numeric(fig8a_graph_input$upper_limit)
fig8a_graph_input<- fig8a_graph_input %>% mutate(Bin_Midpoint = (lower_limit + upper_limit) /2)


#JA:rearrange the columns in graoh input
fig8a_graph_input<- fig8a_graph_input %>% select(Bin, upper_limit, lower_limit, Bin_Midpoint, Mean_Nucleus_signal_intensity, Frequency, Identifier)



#JA: plotting average cells count compared between core wise in each bin
#JA:for that find how many cores we have in each identifier set

normal_core<- length(unique(normal$`TMA core`)) #JA:1
grade1_core <- length(unique(grade1$`TMA core`)) #JA:1
grade2_core<- length(unique(grade2$`TMA core`)) #JA:37
grade3_core<- length(unique(grade3$`TMA core`)) #JA:15


#JA: Step3: Normalization:
#JA: calculating average cells count compared between core wise in each bin
fig8a_graph_input <- fig8a_graph_input %>%
  mutate(average_frequency = case_when(
    Identifier == "Normal" ~ Frequency / normal_core,
    Identifier == "Grade 1" ~ Frequency / grade1_core,
    Identifier == "Grade 2" ~ Frequency / grade2_core,
    Identifier == "Grade 3" ~ Frequency / grade3_core,
    TRUE ~ NA_real_  #JA: Handles cases where Identifier do not match to either "Adjacent",  "grade2", or "grade3"  
  ))

#JA: calculate normalized average frequency                                                                                               
fig8a_graph_input<- fig8a_graph_input %>% 
  mutate(norm_average_frequency = case_when(
    Identifier == "Normal" ~ average_frequency / sum(fig8a_graph_input$average_frequency[fig8a_graph_input$Identifier == "Normal"]),
    Identifier == "Grade 1" ~ average_frequency / sum(fig8a_graph_input$average_frequency[fig8a_graph_input$Identifier == "Grade 1"]),
    Identifier == "Grade 2" ~ average_frequency / sum(fig8a_graph_input$average_frequency[fig8a_graph_input$Identifier == "Grade 2"]),
    Identifier == "Grade 3" ~ average_frequency / sum(fig8a_graph_input$average_frequency[fig8a_graph_input$Identifier == "Grade 3"]),
    TRUE ~ NA_real_  #JA: Handles cases where Identifier do not match to either "Adjacent",  "grade2", or "grade3"
  ))

fig8a_graph_input$Identifier <- factor(fig8a_graph_input$Identifier, levels = c("Normal", "Grade 1", "Grade 2", "Grade 3"))
fwrite(fig8a_graph_input, "BR720_fig8a_grade1_vs_2_vs_3_vs_normal_RNase_geom_point.csv")


#JA: Step4: Visualization
max_value <- round(max(fig8a_graph_input$Bin_Midpoint), 1)



fig8a_graph<- ggplot()+
  geom_density(
    #JA: density plot for "S9.6 in normal"
    data = fig8a_input %>% filter(Identifier == "Normal"),
    aes(
      #JA: x is values before binning
      #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8a_graph_input$norm_average_frequency[fig8a_graph_input$Identifier == "Normal"]) * diff(bin_breaks_75)[1],
      color = "Normal"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
  #JA:plot for grade 1

  geom_density(
    data = fig8a_input %>% filter(Identifier == "Grade 1"),
    aes(
      #JA: x is values before binning
      #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8a_graph_input$norm_average_frequency[fig8a_graph_input$Identifier == "Grade 1"]) * diff(bin_breaks_75)[1],
      color = "Grade 1"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+

  #JA:plot for Grade2
  geom_density(
    data = fig8a_input %>% filter(Identifier == "Grade 2"),
    aes(
      #JA: x is values before binning
      #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8a_graph_input$norm_average_frequency[fig8a_graph_input$Identifier == "Grade 2"]) * diff(bin_breaks_75)[1],
      color = "Grade 2"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
  #JA:plot for grade3
  geom_density(
    data = fig8a_input %>% filter(Identifier == "Grade 3"),
    aes(
      #JA: x is values before binning
      #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8a_graph_input$norm_average_frequency[fig8a_graph_input$Identifier == "Grade 3"]) * diff(bin_breaks_75)[1],
      color = "Grade 3"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
  
  geom_point(
    data = fig8a_graph_input,
    aes(x = Mean_Nucleus_signal_intensity, y = norm_average_frequency, color = Identifier),
    size = 3
  )+
  
  scale_color_manual(
    values = c("Normal" = "blue", "Grade 1" = "darkolivegreen", "Grade 2" = "red", "Grade 3" = "black")
  )+ 
  
  scale_x_continuous(breaks= seq(0, max_value, by= 0.1))+
  
  labs(
    title = "S9.6 signal intensity across different grades in ICNST",
    subtitle ="RNase III+ RNase H",
    x = "DAB intensity mean per nucleus",
    y = "Normalized Frequency",
    color = "Legends"
  ) +
  
  #JA:scale_x_log10()+
  #JA:scale_y_log10()+
  
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(), #JA:facet wrap is not useful because of lot of tma_cores
        legend.position = "top")



#JA:change bw to 0.05
ggsave("BR720_fig8a_ICNST_grade1_vs_2_vs_3_vs_normal_RNase_core_wise_cell_avg_per_bin_normalised_0.05bw.tiff", 
       plot = fig8a_graph, width=15,height=10, dpi=600)

#JA:change bw to 0.02
ggsave("BR720_fig8a_ICNST_grade1_vs_2_vs_3_vs_normal_RNase_core_wise_cell_avg_per_bin_normalised_0.02bw.tiff", 
       plot = fig8a_graph, width=15,height=10, dpi=600)

#JA:Step5: Statistical analysis:
anova_result <- aov(`DAB: Nucleus: Mean` ~ Identifier, data = fig8a_input)
summary (anova_result)
#JA:Df Sum Sq Mean Sq F value Pr(>F)    
#JA:Identifier       3    250   83.35    3036 <2e-16 ***
#JA:Residuals   167830   4608    0.03                   
#JA:---
#JA:Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_result <- TukeyHSD(anova_result) #JA:This output is from the Tukey HSD (Honest Significant Difference) test, which shows the results of pairwise comparisons between the means of the groups in your ANOVA model. Here's a breakdown of each part:

print(tukey_result)
#JA:Tukey multiple comparisons of means
#JA:95% family-wise confidence level

#JA:Fit: aov(formula = `DAB: Nucleus: Mean` ~ Identifier, data = fig8a_input)

#JA:$Identifier
#JA:diff         lwr         upr    p adj
#JA:Grade 1-Normal  -0.02196239 -0.05121523 0.007290439 0.215911
#JA:Grade 2-Normal   0.08690702  0.06047057 0.113343455 0.000000
#JA:Grade 3-Normal   0.16151903  0.13506760 0.187970461 0.000000
#JA:Grade 2-Grade 1  0.10886941  0.09619875 0.121540068 0.000000
#JA:Grade 3-Grade 1  0.18348143  0.17077952 0.196183330 0.000000
#JA:Grade 3-Grade 2  0.07461202  0.07249177 0.076732259 0.000000


#JA:repeat the same for S9.6 dataset too.
#JA:part II: Grade 1 vs 2 vs 3 vs normal with S9.6(Nucleus)
S9.6 <- fread("BR720_IDC_only_S9.6_positive_osk_No_negative.csv", sep = ",", header = TRUE)
#JA:269381


unique(S9.6$Grade)
#JA:[1] "1" "1--2" 2"    "3" "2--3" "-"


S9.6<- S9.6 %>% select("Image", "Name", "Class", "TMA core", "Hematoxylin: Nucleus: Mean", 
                       "DAB: Nucleus: Mean" ,"Pathology diagnosis", "Grade")

#JA:removing all negative values from fig8a_S9.6 input
S9.6<- S9.6 %>% filter(`DAB: Nucleus: Mean` > 0) #JA: there is no negative values in this dataset



#JA:filtering for grade wise
grade1 <- S9.6 %>% filter(Grade == "1")
#JA:2115

grade2<- S9.6 %>% filter(Grade == "2")
#JA:152703

grade3<- S9.6 %>% filter(Grade == "3")
#JA:98688

normal<- fread("BR720_S9.6_N_No negative.csv", header = TRUE, sep = ",")
range(normal$`DAB: Nucleus: Mean`) #JA: no negative values
#JA:[1] 0.0409 0.7994

normal <- normal %>% select("Image", "Name", "Class", "TMA core", "Hematoxylin: Nucleus: Mean", 
                            "DAB: Nucleus: Mean" ,"Pathology diagnosis", "Grade")
#JA:789

#JA:Setting a column as identifier
grade1$Identifier <- "Grade 1"
grade2$Identifier <- "Grade 2"
grade3$Identifier <- "Grade 3"
normal$Identifier <- "Normal"


#JA:Combing all three dataset as one for analyzing them together
fig8a_S9.6_input<- rbind(grade1, grade2, grade3, normal)
#JA:2115+152703+98688+789 = 254295

#JA:First find the range for binning
range<- range(fig8a_S9.6_input$`DAB: Nucleus: Mean`)
#JA:[1]  0.0119 1.3881

a<- range[[1]][1]
b<- range[[2]][1]
bin_breaks_75 <- seq(a, b, length.out=76)

#JA:add bin range in a new column

fig8a_S9.6_input <- fig8a_S9.6_input %>%
  dplyr::mutate(Bin = cut(`DAB: Nucleus: Mean`, breaks = bin_breaks_75, include.lowest = TRUE))
fig8a_S9.6_input$Identifier <- factor(fig8a_S9.6_input$Identifier, levels = c("Normal", "Grade 1", "Grade 2", "Grade 3"))

fwrite(fig8a_S9.6_input, "BR720_fig8a_grade1_vs_2_vs_3_vs_normal_S9.6_density_input.csv")

#JA: Calculate the mean of nucleus signal intensities in that range not necessarily same as mean of midpoint 
fig8a_S9.6_graph_input <- fig8a_S9.6_input %>%
  group_by(Bin, Identifier) %>%
  summarise(
    Mean_Nucleus_signal_intensity = round(mean(`DAB: Nucleus: Mean`, na.rm = TRUE),2), #JA:mean of nucleus intensity calculation within that unique group
    Frequency = n()
  )



#JA:Displaying bin midpoints
fig8a_S9.6_graph_input$Bin_char <- gsub("\\[|\\]|\\(|\\)", "", fig8a_S9.6_graph_input$Bin)
fig8a_S9.6_graph_input <- fig8a_S9.6_graph_input %>%
  separate(Bin_char, c('lower_limit', 'upper_limit'), sep = ",")
fig8a_S9.6_graph_input$lower_limit<- as.numeric(fig8a_S9.6_graph_input$lower_limit)
fig8a_S9.6_graph_input$upper_limit<- as.numeric(fig8a_S9.6_graph_input$upper_limit)
fig8a_S9.6_graph_input<- fig8a_S9.6_graph_input %>% mutate(Bin_Midpoint = (lower_limit + upper_limit) /2)


#JA:rearrange the columns in graoh input
fig8a_S9.6_graph_input<- fig8a_S9.6_graph_input %>% select(Bin, upper_limit, lower_limit, Bin_Midpoint, Mean_Nucleus_signal_intensity, Frequency, Identifier)



#JA: plotting average cells count compared between core wise in each bin
#JA:for that find how many cores we have in each identifier set

normal_core<- length(unique(normal$`TMA core`)) #JA:1
grade1_core <- length(unique(grade1$`TMA core`)) #JA:1
grade2_core<- length(unique(grade2$`TMA core`)) #JA:37
grade3_core<- length(unique(grade3$`TMA core`)) #JA:15

#JA: calculating average cells count compared between core wise in each bin
fig8a_S9.6_graph_input <- fig8a_S9.6_graph_input %>%
  mutate(average_frequency = case_when(
    Identifier == "Normal" ~ Frequency / normal_core,
    Identifier == "Grade 1" ~ Frequency / grade1_core,
    Identifier == "Grade 2" ~ Frequency / grade2_core,
    Identifier == "Grade 3" ~ Frequency / grade3_core,
    TRUE ~ NA_real_  #JA: Handles cases where Identifier do not match to either "Adjacent",  "grade2", or "grade3"  
  ))

#JA:calculate normalized average frequency                                                                                               
fig8a_S9.6_graph_input<- fig8a_S9.6_graph_input %>% 
  mutate(norm_average_frequency = case_when(
    Identifier == "Normal" ~ average_frequency / sum(fig8a_S9.6_graph_input$average_frequency[fig8a_S9.6_graph_input$Identifier == "Normal"]),
    Identifier == "Grade 1" ~ average_frequency / sum(fig8a_S9.6_graph_input$average_frequency[fig8a_S9.6_graph_input$Identifier == "Grade 1"]),
    Identifier == "Grade 2" ~ average_frequency / sum(fig8a_S9.6_graph_input$average_frequency[fig8a_S9.6_graph_input$Identifier == "Grade 2"]),
    Identifier == "Grade 3" ~ average_frequency / sum(fig8a_S9.6_graph_input$average_frequency[fig8a_S9.6_graph_input$Identifier == "Grade 3"]),
    TRUE ~ NA_real_  #JA: Handles cases where Identifier do not match to either "Adjacent",  "grade2", or "grade3"
  ))

fig8a_S9.6_graph_input$Identifier <- factor(fig8a_S9.6_graph_input$Identifier, levels = c("Normal", "Grade 1", "Grade 2", "Grade 3"))
fwrite(fig8a_S9.6_graph_input, "BR720_fig8a_grade1_vs_2_vs_3_vs_normal_S9.6_geom_point.csv")

max_value <- round(max(fig8a_S9.6_graph_input$Bin_Midpoint), 1)



fig8a_S9.6_graph<- ggplot()+
  geom_density(
    #JA: density plot for "S9.6 in normal"
    data = fig8a_S9.6_input %>% filter(Identifier == "Normal"),
    aes(
      #JA: x is values before binning
      #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8a_S9.6_graph_input$norm_average_frequency[fig8a_S9.6_graph_input$Identifier == "Normal"]) * diff(bin_breaks_75)[1],
      color = "Normal"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
  #JA:plot for grade 1
  
  geom_density(
    data = fig8a_S9.6_input %>% filter(Identifier == "Grade 1"),
    aes(
      #JA: x is values before binning
      #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8a_S9.6_graph_input$norm_average_frequency[fig8a_S9.6_graph_input$Identifier == "Grade 1"]) * diff(bin_breaks_75)[1],
      color = "Grade 1"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
  #JA:plot for Grade2
  geom_density(
    data = fig8a_S9.6_input %>% filter(Identifier == "Grade 2"),
    aes(
      #JA: x is values before binning
      #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8a_S9.6_graph_input$norm_average_frequency[fig8a_S9.6_graph_input$Identifier == "Grade 2"]) * diff(bin_breaks_75)[1],
      color = "Grade 2"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
  #JA:plot for grade3
  geom_density(
    data = fig8a_S9.6_input %>% filter(Identifier == "Grade 3"),
    aes(
      #JA: x is values before binning
      #JA: y is probability density scaled to total number of records * bin_width
      x = `DAB: Nucleus: Mean`,
      y = ..density..* sum(fig8a_S9.6_graph_input$norm_average_frequency[fig8a_S9.6_graph_input$Identifier == "Grade 3"]) * diff(bin_breaks_75)[1],
      color = "Grade 3"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
  
  geom_point(
    data = fig8a_S9.6_graph_input,
    aes(x = Mean_Nucleus_signal_intensity, y = norm_average_frequency, color = Identifier),
    size = 3
  )+
  
  scale_color_manual(
    values = c("Normal" = "blue", "Grade 1" = "darkolivegreen", "Grade 2" = "red", "Grade 3" = "black")
  )+ 
  
  scale_x_continuous(breaks= seq(0, max_value, by= 0.1))+
  
  labs(
    title = "S9.6 signal intensity across different grades in ICNST",
    subtitle ="S9.6 Only",
    x = "DAB intensity mean per nucleus",
    y = "Normalized Frequency",
    color = "Legends"
  ) +
  
  #JA:scale_x_log10()+
  #JA:scale_y_log10()+
  
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(), #JA:facet wrap is not useful because of lot of tma_cores
        legend.position = "top")


ggsave("BR720_fig8a_ICNST_grade1_vs_2_vs_3_vs_normal_S9.6_core_wise_cell_avg_per_bin_normalised_0.05bw.tiff",
       plot = fig8a_S9.6_graph, width=15,height=10, dpi=600)


#JA:change bw to 0.02
ggsave("BR720_fig8a_ICNST_grade1_vs_2_vs_3_vs_normal_S9.6_core_wise_cell_avg_per_bin_normalised_0.02bw.tiff",
       plot = fig8a_S9.6_graph, width=15,height=10, dpi=600)

#JA:Statistics:
anova_result <- aov(`DAB: Nucleus: Mean` ~ Identifier, data = fig8a_S9.6_input)
summary (anova_result)
#JA:                Df Sum Sq Mean Sq F value Pr(>F)    
#JA:Identifier       3    405  134.93    3831 <2e-16 ***
#JA: Residuals   254291   8957    0.04                   
#JA:---
#JA:Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

tukey_result <- TukeyHSD(anova_result) #JA:This output is from the Tukey HSD (Honest Significant Difference) test, which shows the results of pairwise comparisons between the means of the groups in your ANOVA model. Here's a breakdown of each part:

print(tukey_result)
#JA:Tukey multiple comparisons of means
#JA:95% family-wise confidence level

#JA:Fit: aov(formula = `DAB: Nucleus: Mean` ~ Identifier, data = fig8a_S9.6_input)

#JA:$Identifier
#JA:diff         lwr         upr    p adj
#JA:Grade 1-Normal  -0.05969786 -0.07981120 -0.03958451     0
#JA:Grade 2-Normal   0.04053575  0.02332657  0.05774494     0
#JA:Grade 3-Normal   0.11856581  0.10133244  0.13579919     0
#JA:Grade 2-Grade 1  0.10023361  0.08967732  0.11078990     0
#JA:Grade 3-Grade 1  0.17826367  0.16766798  0.18885936     0
#JA:Grade 3-Grade 2  0.07803006  0.07606082  0.07999930     0



