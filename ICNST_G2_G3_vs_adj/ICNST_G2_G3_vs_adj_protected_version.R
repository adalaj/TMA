# ------------------------------------------------------------------------------
# This code is part of paper: RNA:DNA hybrids determine tissue pathological heterogeneity and aggressiveness in clinically diagnosed cancers. 
# Author: Jyoti Devendra Adala under supervision of Dr. Vladimir A Kuznetsov 
# For updates and contributions, visit : https://github.com/adalaj

# Purpose: This code analyzes the S9.6 signal intensity across different grades in invasive carcinoma of no special type (ICNST). 
#The script processes the dataset, normalizes it, and visualizes the results to compare nucleus signal intensity distributions across three categories: "Adjacent normal tissue," "Grade 2," and "Grade 3."
# ------------------------------------------------------------------------------


#JA: load required libraries
library(data.table)
library(tidyverse)

##JA: Step 1: Data Preprocessing:
#JA:load the dataset
BC08112 <- fread("BC081120f_MERGED_SORTED_No Negative_OSK.csv")
#670253


#JA: First find the range for binning to remove negative values
range(BC08112$`Nucleus: DAB OD mean`)
#[1] -0.4313  1.2655

#JA: removing all negative values 
BC08112<- BC08112 %>% filter(`Nucleus: DAB OD mean` > 0)
#661271


#JA: filtering for two or more criteria based on pathology and grade
adjacent <- BC08112 %>% filter(`Pathology diagnosis` %in% c("Adjacent normal breast tissue (fibrofatty tissue)",  "Adjacent normal breast tissue", "Adjacent normal breast tissue (adenosis)"))
#JA:15900 rows

grade2<- BC08112 %>% filter(Grade == "2")
#JA:382122

grade3<- BC08112 %>% filter(Grade == "3")
#JA: 249094


#JA: Setting a column as identifier
adjacent$Identifier <- "Adjacent"
grade2$Identifier <- "Grade 2"
grade3$Identifier <- "Grade 3"

#JA:Combing all three dataset as one for analyzing them together
fig7_input<- rbind(adjacent, grade2, grade3)
#15900+382122+249094 = 647116



##JA: Step2: Binning and aggregation:
#JA: First find the range for binning
range<- range(fig7_input$`Nucleus: DAB OD mean`)
#[1] 0.0001  1.2655


a<- range[[1]][1]
b<- range[[2]][1]
bin_breaks_75 <- seq(a, b, length.out=76)

#JA: add bin range in a new column
fig7_input <- fig7_input %>%
  dplyr::mutate(Bin = cut(`Nucleus: DAB OD mean`, breaks = bin_breaks_75, include.lowest = TRUE))

#JA: saving file
fwrite(fig7_input, "BC08112_Fig7_S9.6_ICNST_grade2_vs_3_vs_adj_density_input.csv")

#JA: Calculate the mean of nucleus signal intensities in that range not necessarily same as mean of midpoint 
fig7_graph_input <- fig7_input %>%
  group_by(Bin, Identifier) %>%
  summarise(
    Mean_Nucleus_signal_intensity = round(mean(`Nucleus: DAB OD mean`, na.rm = TRUE),2), #JA:mean of nucleus intensity calculation within that unique group
    Frequency = n()
  )



#JA: Displaying bin midpoints
fig7_graph_input$Bin_char <- gsub("\\[|\\]|\\(|\\)", "", fig7_graph_input$Bin)
fig7_graph_input <- fig7_graph_input %>%
  separate(Bin_char, c('lower_limit', 'upper_limit'), sep = ",")
fig7_graph_input$lower_limit<- as.numeric(fig7_graph_input$lower_limit)
fig7_graph_input$upper_limit<- as.numeric(fig7_graph_input$upper_limit)
fig7_graph_input<- fig7_graph_input %>% mutate(Bin_Midpoint = (lower_limit + upper_limit) /2)


#JA: rearrange the columns in graph input
fig7_graph_input<- fig7_graph_input %>% select(Bin, upper_limit, lower_limit, Bin_Midpoint, Mean_Nucleus_signal_intensity, Frequency, Identifier)



#JA: plotting average cells count compared between core wise in each bin
#for that find how many cores we have in each identifier set

adj_core<- length(unique(adjacent$`TMA core`))
grade2_core<- length(unique(grade2$`TMA core`))
grade3_core<- length(unique(grade3$`TMA core`))



##JA: Step3: Normalization
#JA: calculating normalized average frequency of cells count within each bin for different categories. 
fig7_graph_input <- fig7_graph_input %>%
  mutate(average_frequency = case_when(
    Identifier == "Adjacent" ~ Frequency / adj_core,
    Identifier == "Grade 2" ~ Frequency / grade2_core,
    Identifier == "Grade 3" ~ Frequency / grade3_core,
    TRUE ~ NA_real_  #JA:Handles cases where Identifier do not match to either "Adjacent",  "grade2", or "grade3"  
    ))

##JA:calculate normalized average frequency                                                                                               
fig7_graph_input<- fig7_graph_input %>% 
  mutate(norm_average_frequency = case_when(
    Identifier == "Adjacent" ~ average_frequency / sum(fig7_graph_input$average_frequency[fig7_graph_input$Identifier == "Adjacent"]),
    Identifier == "Grade 2" ~ average_frequency / sum(fig7_graph_input$average_frequency[fig7_graph_input$Identifier == "Grade 2"]),
    Identifier == "Grade 3" ~ average_frequency / sum(fig7_graph_input$average_frequency[fig7_graph_input$Identifier == "Grade 3"]),
    TRUE ~ NA_real_  # JA: Handles cases where Identifier do not match to either "Adjacent",  "grade2", or "grade3"
  ))

#JA: saving the file
fwrite(fig7_graph_input, "BC08112_Fig7_S9.6_ICNST_grade2_vs_3_vs_adj_geom_point.csv")




##JA: Step4: Visualization
#JA: max value
max_value <- round(max(fig7_graph_input$Bin_Midpoint), 1)

#JA: plotting of graph begins
fig7_graph<- ggplot()+
  geom_density(
    #JA: density plot for "S9.6 in adjacent"
    data = fig7_input %>% filter(Identifier == "Adjacent"),
    aes(
      #JA: x is values before binning
      #JA: y is probability density scaled to total number of records * bin_width
      x = `Nucleus: DAB OD mean`,
      y = ..density..* sum(fig7_graph_input$norm_average_frequency[fig7_graph_input$Identifier == "Adjacent"]) * diff(bin_breaks_75)[1],
      color = "Adjacent"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
  
  #JA:plot for Grade2
  geom_density(
    data = fig7_input %>% filter(Identifier == "Grade 2"),
    aes(
      # JA: x is values before binning
      # JA: y is probability density scaled to total number of records * bin_width
      x = `Nucleus: DAB OD mean`,
      y = ..density..* sum(fig7_graph_input$norm_average_frequency[fig7_graph_input$Identifier == "Grade 2"]) * diff(bin_breaks_75)[1],
      color = "Grade 2"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
  #JA:plot for grade3
  geom_density(
    data = fig7_input %>% filter(Identifier == "Grade 3"),
    aes(
      # JA: x is values before binning
      # JA: y is probability density scaled to total number of records * bin_width
      x = `Nucleus: DAB OD mean`,
      y = ..density..* sum(fig7_graph_input$norm_average_frequency[fig7_graph_input$Identifier == "Grade 3"]) * diff(bin_breaks_75)[1],
      color = "Grade 3"
    ),
    size = 1,
    kernel = "gaussian",
    bw = 0.05)+
  
  
  geom_point(
    data = fig7_graph_input,
    aes(x = Mean_Nucleus_signal_intensity, y = norm_average_frequency, color = Identifier),
    size = 3
  )+
  
  scale_color_manual(
    values = c("Adjacent" = "blue", "Grade 2" = "red", "Grade 3" = "black")
  )+ 
  
  scale_x_continuous(breaks= seq(0, max_value, by= 0.1))+
  
  labs(
    title = "S9.6 signal intensity across different grades in ICNST",
    x = "DAB intensity mean per nucleus",
    y = "Normalized Frequency",
    color = "Legends"
  ) +
  
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        panel.grid = element_blank(),
        panel.background = element_blank(), 
        legend.position = "top")


ggsave("BC08112_Fig7_S9.6_ICNST_grade2_vs_3_vs_adj_core_wise_cell_avg_per_bin_normalised.tiff", 
       plot = fig7_graph, width=11,height=10, dpi=600)

##JA: Step6: Statiscal Analysis.
#JA: Calculate Statistics:
#anova_result <- aov(`Nucleus: DAB OD mean` ~ Identifier, data = fig7_input)
#summary (anova_result)
#Df Sum Sq Mean Sq F value Pr(>F)    
#Identifier       2   1313   656.7   16082 <2e-16 ***
#Residuals   647113  26425     0.0                   
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

