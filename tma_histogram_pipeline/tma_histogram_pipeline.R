# ------------------------------------------------------------------------------
#This code is part of paper: RNA:DNA hybrids determine tissue pathological heterogeneity and aggressiveness in clinically diagnosed cancers. 
#Author: Jyoti Devendra Adala under supervision of Dr. Vladimir A Kuznetsov 
#For updates and contributions, visit : https://github.com/adalaj

#Purpose: To create two or more histogram plots using the parameters provided by users. 
#JA: The parameters are Nucleus: DAB OD mean and Cytoplasm: DAB OD mean are greater than zero, plot first graph
#JA: and next is The parameters are Nucleus: DAB OD mean and Cytoplasm: DAB OD mean are greater than zero, and Nucleus: Max caliper is <25 along with Cell: DAB OD mean >0.18, plot next graph

#JA: input is provided by users in characters "14895_TMA1.2401b_Breast_core_G-1_S9.6.csv"
#JA: user will paste this R script in the same working directory where files are there.
#The outputs include density plots, statistical summaries, and files for further analysis.
# ------------------------------------------------------------------------------




#JA:user needs to install below package

#JA:install.packages("tidyverse")
library(tidyverse)

#JA:install.packages("data.table")
library(data.table)

#JA:install.packages("fitdistrplus")
library(fitdistrplus)

args<- commandArgs(trailingOnly = TRUE)
input <- args[1]                #JA: First argument: CSV file name
#JA:input<- "14895_TMA1.2401b_Breast_core_G-1_S9.6.csv"

bin_size <- as.numeric(args[2])       #JA: Second argument: Bin size
#JA:user will provide bin size
#JA:bin_size <- c(75)


col1<- args[3]
#JA:col1<- c("Nucleus: DAB OD mean")

col1_cutoff <- as.numeric(args[4])  #JA: Third argument: cutoff for col1
#JA:col1_cutoff <- 0

col2<- args[5]
#JA:col2<- c("Cytoplasm: DAB OD mean")

col2_cutoff <- as.numeric(args[6])  #JA: Fourth argument: cutoff for col2
#JA:col2_cutoff <- 0

graph_title<- as.character(args[7])

col3<- args[8]
#JA:col3<- c("Nucleus: Max caliper")

col3_cutoff <- as.numeric(args[9])  #JA: Fifth argument: cutoff for col3
#JA:col3_cutoff<- 25

col4<- args[10]
#JA:col4<- c("Cell: DAB OD mean")

col4_cutoff <- as.numeric(args[11])  #JA: Sixth argument: cutoff for col4
#JA:col4_cutoff<- 0.18

graph_title_after_cutoff <- as.character(args[12])

print(paste("Below parameters are defined by user for file", input))

print(bin_size)
print(col1)
print(col1_cutoff)

print(col2)
print(col2_cutoff)

print(col3)
print(col3_cutoff)

print(col4)
print(col4_cutoff)

#JA: Cutoff values for filtering parameters

#JA: Define the working directory as current directory (or customize this)
setwd(getwd())



filename<- fread(input, sep = ",", header = TRUE)
#JA:select columns: provided by users:

column_names<- c("V1","TMA core", "Nucleus: DAB OD mean", "Cytoplasm: DAB OD mean", "Nucleus: Max caliper", "Cell: DAB OD mean")

tma<- filename %>% dplyr::select(all_of(column_names))

bin_size_new<- bin_size+1


#JA:when nucleus and cytoplasm DAB OD mean is >0 and no additional filters are applied


positive_tma<- tma[
  tma[[col1]]>0 &
    tma[[col2]]>0, 
]

fwrite(positive_tma, paste("TMA with positive cells", input, ".csv"))
print (paste("Numbers of positive cells:",  nrow(positive_tma), "of",input, sep = " "))

negative_tma<- tma[
  !(tma[[col1]]>0 &
      tma[[col2]]>0),
]

fwrite(negative_tma, paste("TMA with negative cells", input, ".csv"))
print (paste("Numbers of negative cells:",  nrow(negative_tma), "of", input, sep = " "))

nucleus_dab <- positive_tma$`Nucleus: DAB OD mean`
cytoplasm_dab <- positive_tma$`Cytoplasm: DAB OD mean`

print("Calculating statistics")

mean_col1 <- round(mean(positive_tma[[col1]]), 2)
sd_col1 <- round(sd(positive_tma[[col1]]), 2)
median_col1 <- round(median(positive_tma[[col1]]), 2)

print(paste("Mean of", col1, "is", mean_col1))
print(paste("Standard deviation of", col1, "is", sd_col1))
print(paste("Median of", col1, "is",median_col1))


mean_col2 <- round(mean(positive_tma[[col2]]), 2)
sd_col2 <- round(sd(positive_tma[[col2]]), 2)
median_col2 <- round(median(positive_tma[[col2]]), 2)

print(paste("Mean of", col2, "is", mean_col2))
print(paste("Standard deviation of", col2, "is", sd_col2))
print(paste("Median of", col2, "is", median_col2))

#JA: Calculating KS test
ks_test_result <- ks.test(nucleus_dab, cytoplasm_dab)
ks_statistic <- round(ks_test_result$statistic, 5)
ks_pvalue <- format.pval(ks_test_result$p.value, digits = 2, eps = .Machine$double.eps)
print(paste("Asymptotic two-sample Kolmogorov-Smirnov test calculated between", col1, "and", col2, "is: D =", ks_statistic, "with the p-value of", ks_pvalue, "of", input))


#JA:Calculating t test

t_test_result<- t.test(nucleus_dab, cytoplasm_dab)
#JA: Extract t-statistic (t), t-degree of freedom(df) and p-value
t_statistic <- round(t_test_result$statistic, 5)
t_df<- round(t_test_result$parameter,5)
t_pvalue <- format.pval(t_test_result$p.value, digits = 2, eps = .Machine$double.eps)

print(paste("Welch Two Sample t-test between", col1, "and", col2, "of", input,"is as follows:"), sep="")
print(t_test_result)



#JA:Calculating ANOVA:
tma_long <- positive_tma %>%
  pivot_longer(cols = c(`Nucleus: DAB OD mean`, `Cytoplasm: DAB OD mean`), 
               names_to = "Variable", values_to = "Value")

anova_result <- aov(Value ~Variable, tma_long)
print("Summary of ANOVA results:")
print(summary(anova_result))


#JA: Calculating Euclidean distance
euclidean_distance <- round(sqrt(sum((nucleus_dab - cytoplasm_dab)^2)),2)
print(paste("Euclidean_distance calculated between Nucleus DAB OD mean and cytoplasm DAB OD mean", "of", input, "is", euclidean_distance, sep = " "))

a1 <- min(min(positive_tma[[col1]]), min(positive_tma[[col2]]))
print(paste("minimum values between", col1, "and", col2, "is", a1, sep=" "))


b1 <- max(max(positive_tma[[col1]]), max(positive_tma[[col2]]))
print(paste("maximum value between", col1, "and", col2, "is", b1, sep=" "))

print(paste(bin_size,"bins are created between range", a1, "and", b1))

bin_break<- seq(a1, b1, length.out=bin_size_new)

print("Saving all statistics in a table format")

result_table <- data.frame(
  Metric = c("Mean of Col1", "Standard Deviation of Col1", "Median of Col1",
             "Mean of Col2", "Standard Deviation of Col2", "Median of Col2",
             "KS Statistic(D)", "KS P-value", "T-test(t)","T-test(df)","T-test(Pvalue)", "Euclidean Distance", 
             "Min Value", "Max Value", "Number of Bins", "Bin Range", "Number of positive cells", "Number of negative cells", "Total number of cells"),
  Value = c(mean_col1, sd_col1, median_col1,
            mean_col2, sd_col2, median_col2,
            ks_statistic, ks_pvalue, t_statistic, t_df, t_pvalue, euclidean_distance,
           a1, b1, bin_size_new, paste(round(bin_break[1],2), "-", round(bin_break[length(bin_break)],2)), nrow(positive_tma), nrow(negative_tma), nrow(tma))
)


print("Printing statistics table and ANOVA results are printed")

fwrite(result_table, paste("Statistics_table_for_", graph_title, ".csv"))


hist_col1<- hist(positive_tma[[col1]], breaks = bin_break, plot = FALSE)
hist_col2 <- hist(positive_tma[[col2]], breaks = bin_break, plot = FALSE)

hist_col1_df <- data.frame(
  midpoints = (head(hist_col1$breaks, -1) + tail(hist_col1$breaks, -1)) / 2,
  counts = hist_col1$counts,
  Variable = col1
)

hist_col2_df <- data.frame(
  midpoints = (head(hist_col2$breaks, -1) + tail(hist_col2$breaks, -1)) / 2,
  counts = hist_col2$counts,
  Variable = col2
)

#JA:#JA:define Normalized values 
hist_col1_df<- hist_col1_df %>% mutate(norm_counts = counts/sum(counts))
hist_col2_df<- hist_col2_df %>% mutate(norm_counts = counts/sum(counts))

core_final<- cbind(hist_col1_df, hist_col2_df)

#JA:#JA:saving graph input
fwrite(core_final, paste("graphinput no cutoff", input, sep = " "))
fwrite(tma_long, paste("geom density graphinput no cutoff", input, sep = " "))


#JA:#JA:defining graph axis
max_param_value <- max(max(positive_tma[[col1]]), max(positive_tma[[col2]]))
graph_max_value <- round(max_param_value, 1)+0.1

#JA:#JA:#JA:#JA:#JA:Plotting graphs with absolute frequency: 
print(paste("Plotting the graph between", col1, "and", col2, "when only Nucleus and Cytoplasm DAB OD mean has positive values"))

print("Printing probability density function with absolute frequency")
two_histo_approximation <- ggplot() +
 geom_bar(data = hist_col1_df, 
           aes(x = midpoints, y = counts, fill = "col1"), 
           stat = "identity", color = "purple", alpha = 0.7) + 
  geom_bar(data = hist_col2_df, 
           aes(x= midpoints, y= counts, fill= "col2"), 
           stat= "identity", color = "darkgreen", alpha= 0.7) +
  geom_density(data= tma_long %>% filter (Variable == col1), 
               aes(x = Value, y = ..density.. * sum(hist_col1_df$counts) * diff(bin_break)[1], color = "col1"), size = 1, kernel = "gaussian", bw = 0.09) +
  geom_density(data= tma_long %>% filter (Variable == col2), 
               aes(x = Value, y = ..density.. * sum(hist_col2_df$counts) * diff(bin_break)[1], color = "col2"), size = 1, kernel = "gaussian", bw = 0.09) +
  coord_cartesian(xlim = c(0, graph_max_value)) +
  scale_x_continuous(breaks = seq(0, graph_max_value, by = 0.2)) +
  labs(title = graph_title, subtitle = "Value distribution fitted by probability density function", x = "DAB OD mean", y = "Frequency") +
  scale_fill_manual(name = "Legend", values = c("col1" = "purple", "col2" = "darkgreen"), labels= c(col1, col2)) +
  scale_color_manual(name = "Legend", values = c("col1" = "purple", "col2" = "darkgreen"),labels= c(col1, col2)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        legend.position = "top")

ggsave(paste("probability density graph no cutoff", input, ".tiff", sep = " "), 
       plot = two_histo_approximation, width=11,height=10, dpi=600)



print("Printing probability density function with Normalized frequency")
two_histo_approximation_norm <- ggplot() +
  geom_bar(data = hist_col1_df, 
           aes(x = midpoints, y = norm_counts, fill = "col1"), 
           stat = "identity", color = "purple", alpha = 0.7) + 
  geom_bar(data = hist_col2_df, 
           aes(x= midpoints, y= norm_counts, fill= "col2"), 
           stat= "identity", color = "darkgreen", alpha= 0.7) +
  geom_density(data= tma_long %>% filter (Variable == col1), 
               aes(x = Value, y = ..density.. * sum(hist_col1_df$norm_counts) * diff(bin_break)[1], color = "col1"), size = 1, kernel = "gaussian", bw = 0.09) +
  geom_density(data= tma_long %>% filter (Variable == col2), 
               aes(x = Value, y = ..density.. * sum(hist_col2_df$norm_counts) * diff(bin_break)[1], color = "col2"), size = 1, kernel = "gaussian", bw = 0.09) +
  coord_cartesian(xlim = c(0, graph_max_value)) +
  scale_x_continuous(breaks = seq(0, graph_max_value, by = 0.2)) +
  labs(title = graph_title, subtitle = "Value distribution fitted by Weibull probability function", x = "DAB OD mean", y = "Normalized frequency") +
  scale_fill_manual(name = "Legend", values = c("col1" = "purple", "col2" = "darkgreen"), labels= c(col1, col2)) +
  scale_color_manual(name = "Legend", values = c("col1" = "purple", "col2" = "darkgreen"),labels= c(col1, col2)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        legend.position = "top")

ggsave(paste("probability density graph Normalized no cutoff", input, ".tiff", sep = " "), 
       plot = two_histo_approximation_norm, width=11,height=10, dpi=600)



print("printing weibull probability function with absolute frequency")
#JA:#JA:defining other parameters needed for plotting weibull plotting
fit_weibull_col1 <- fitdist(positive_tma[[col1]], "weibull")
fit_weibull_col2 <- fitdist(positive_tma[[col2]], "weibull")

bin_size_new_col1 <- hist_col1_df$midpoints
bin_size_new_col2 <- hist_col2_df$midpoints

#JA:#JA:all are list values hence not saved 
weibull_density_col1 <- dweibull(bin_size_new_col1, shape = fit_weibull_col1$estimate["shape"], scale = fit_weibull_col1$estimate["scale"])
weibull_density_col2 <- dweibull(bin_size_new_col2, shape = fit_weibull_col2$estimate["shape"], scale = fit_weibull_col2$estimate["scale"])

col1_area <- sum(hist_col1_df$counts * diff(bin_break)[1])
col2_area <- sum(hist_col2_df$counts * diff(bin_break)[1])



paste("Printing and saving the statistics for Weibull")

print(fit_weibull_col1)
print(fit_weibull_col2)


fit_weibull_dt <- data.table(
  Parameter = c("shape",  "shape_SD", "scale","scale_SD","method", "log-likelihood", "AIC (Akaike Information Criterion)", "BIC (Bayesian Information Criterion)", "sample_size", "distname"),
  col1 = c(
    fit_weibull_col1$estimate["shape"],
    fit_weibull_col1$sd["shape"],
    fit_weibull_col1$estimate["scale"],
    fit_weibull_col1$sd["scale"],
    fit_weibull_col1$method,
    fit_weibull_col1$loglik,
    fit_weibull_col1$aic,
    fit_weibull_col1$bic,
    fit_weibull_col1$n,
    fit_weibull_col1$distname
  )
)

colnames(fit_weibull_dt)[2] <- paste(graph_title, col1, sep = "_")

  





#JA:#JA:
fit_weibull_dt2 <- data.table(
  col2 = c(
    fit_weibull_col2$estimate["shape"],
    fit_weibull_col2$sd["shape"],
    fit_weibull_col2$estimate["scale"],
    fit_weibull_col2$sd["scale"],
    fit_weibull_col2$method,
    fit_weibull_col2$loglik,
    fit_weibull_col2$aic,
    fit_weibull_col2$bic,
    fit_weibull_col2$n,
    fit_weibull_col2$distname
  )
)


colnames(fit_weibull_dt2)[1] <- paste(graph_title, col2, sep = "_")


weibull_parameters<- cbind(fit_weibull_dt, fit_weibull_dt2)
fwrite(weibull_parameters, paste("Weibull_parameters_", graph_title, ".csv"))

#JA:#JA:calculating and save the correlation matrix
cor_matrix_dt <- as.data.table(fit_weibull_col1$cor)
cor_matrix_dt$Identifier <- paste(graph_title, col1, sep = "_")

cor_matrix_dt2 <- as.data.table(fit_weibull_col2$cor)
cor_matrix_dt2$Identifier <- paste(graph_title, col2, sep = "_")

cor_matrix_parameters <- rbind(cor_matrix_dt, cor_matrix_dt2)
fwrite(cor_matrix_parameters, paste("Correlation_matrix_for_Weibull_parameters_", graph_title, ".csv"))

#JA:Save the variance-covariance matrix
vcov_matrix_dt <- as.data.table(fit_weibull_col1$vcov)
vcov_matrix_dt$Identifier <- paste(graph_title, col1, sep = "_")

vcov_matrix_dt2 <- as.data.table(fit_weibull_col2$vcov)
vcov_matrix_dt2$Identifier <- paste(graph_title, col2, sep = "_")

vcov_matrix_parameters <- rbind(vcov_matrix_dt, vcov_matrix_dt2)
fwrite(vcov_matrix_parameters, paste("Variance_covariance_matrix_for_Weibull_parameters_", graph_title, ".csv"))



#JA:#JA:plotting the graphs
two_histo_weibull<- ggplot() +
  geom_bar(data = hist_col1_df, 
           aes(x = midpoints, y = counts, fill = "col1"), 
           stat = "identity", color = "purple", alpha = 0.7) + 
  geom_bar(data = hist_col2_df, 
           aes(x= midpoints, y= counts, fill= "col2"), 
           stat= "identity", color = "darkgreen", alpha= 0.7) +
  geom_line(data= hist_col1_df, 
            aes(x = bin_size_new_col1, y = weibull_density_col1*col1_area, color = "col1"), size = 1) + 
  geom_line(data= hist_col2_df, 
            aes(x = bin_size_new_col2, y = weibull_density_col2*col2_area, color = "col2"), size = 1) +
  coord_cartesian(xlim = c(0, graph_max_value)) +
  scale_x_continuous(breaks = seq(0, graph_max_value, by = 0.2)) +
  labs(title = graph_title, subtitle = "Value distribution fitted by Weibull probability function", x = "DAB OD mean", y = "Frequency") +
  scale_fill_manual(name = "Legend", values = c("col1" = "purple", "col2" = "darkgreen"), labels= c(col1, col2)) +
  scale_color_manual(name = "Legend", values = c("col1" = "purple", "col2" = "darkgreen"),labels= c(col1, col2)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        legend.position = "top")

ggsave(paste("Weibull graph no cutoff", input, ".tiff", sep = " "), 
       plot = two_histo_weibull, width=11,height=10, dpi=600)


print("printing weibull probability function with Normalized frequency")
norm_col1_area <- sum(hist_col1_df$norm_counts * diff(bin_break)[1])
norm_col2_area <- sum(hist_col2_df$norm_counts * diff(bin_break)[1])

two_histo_weibull_norm<- ggplot() +
  geom_bar(data = hist_col1_df, 
           aes(x = midpoints, y = norm_counts, fill = "col1"), 
           stat = "identity", color = "purple", alpha = 0.7) + 
  geom_bar(data = hist_col2_df, 
           aes(x= midpoints, y= norm_counts, fill= "col2"), 
           stat= "identity", color = "darkgreen", alpha= 0.7) +
  geom_line(data= hist_col1_df, 
            aes(x = bin_size_new_col1, y = weibull_density_col1*norm_col1_area, color = "col1"), size = 1) + 
  geom_line(data= hist_col2_df, 
            aes(x = bin_size_new_col2, y = weibull_density_col2*norm_col2_area, color = "col2"), size = 1) +
  coord_cartesian(xlim = c(0, graph_max_value)) +
  scale_x_continuous(breaks = seq(0, graph_max_value, by = 0.2)) +
  labs(title = graph_title, subtitle = "Value distribution fitted by Weibull probability function", x = "DAB OD mean", y = "Normalized Frequency") +
  scale_fill_manual(name = "Legend", values = c("col1" = "purple", "col2" = "darkgreen"), labels= c(col1, col2)) +
  scale_color_manual(name = "Legend", values = c("col1" = "purple", "col2" = "darkgreen"),labels= c(col1, col2)) +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        text = element_text(size = 30),
        axis.line = element_line(color = "black"),
        legend.position = "top")

ggsave(paste("Weibull graph Normalized no cutoff", input, ".tiff", sep = " "), 
       plot = two_histo_weibull_norm, width=11,height=10, dpi=600)

print("Printing complete check your output directory folder")



#JA:#JA:#JA:if user provide cutoff value then proceed with this part of analysis

if((!is.null(col3) && !is.na(col3)) || (!is.null(col4) && !is.na(col4))) {
  col3_condition <- TRUE
  col4_condition <- TRUE #JA:Check if user have defined cell_dab_od_mean, otherwise default that is no filtering at all
  print("Executing the next part of the code because Nucleus: Max caliper or Cell: DAB OD mean cutoffs are provided")
  
  if (!is.null(col3) && !is.na(col3)) {
    col3_condition <- tma[[col3]] < col3_cutoff
    print("col3 cutoff is defined and accepted")
  } else{
    print("col3 cutoff is not defined and analysis complete")
  }
  
  
  
  if (!is.null(col4) && !is.na(col4)) {
    col4_condition <- tma[[col4]] > col4_cutoff
    print("col4 cutoff is defined and accepted")
  } else {
    print("col4 cutoff is not defined and analysis complete")
  }
  
  filtered_tma<- tma[
    tma[[col1]] > col1_cutoff &
      tma[[col2]] > col2_cutoff &
      col3_condition &
      col4_condition, 
  ]
  
  
  fwrite(filtered_tma, "TMA matched with user specified parameter cutoffs.csv")
  print (paste("Numbers of cells after user specified cutoff:",  nrow(filtered_tma), sep = " "))
  
  removed_tma<- tma[
    !tma[[col1]]>col1_cutoff &
      tma[[col2]]>col2_cutoff &
      col3_condition &
      col4_condition, 
  ]
  
  fwrite(removed_tma, "TMA that do not matched with user specified parameter cutoffs.csv")
  filtered_nucleus_dab <- filtered_tma$`Nucleus: DAB OD mean`
  filtered_cytoplasm_dab <- filtered_tma$`Cytoplasm: DAB OD mean`
  
  print("Calculating statistics")
  
  print(paste("Mean of", col1, "is", round(mean(filtered_tma[[col1]]),2)))
  print(paste("Standard deviation of", col1, "is", round(sd(filtered_tma[[col1]]),2)))
  print(paste("Median of", col1, "is", round(median(filtered_tma[[col1]]),2)))
  
  print(paste("Mean of", col2, "is", round(mean(filtered_tma[[col2]]),2)))
  print(paste("Standard deviation of", col2, "is", round(sd(filtered_tma[[col2]]),2)))
  print(paste("Median of", col2, "is", round(median(filtered_tma[[col2]]),2)))
  
  #JA: Perform KS test
  ks_test_result <- ks.test(filtered_nucleus_dab, filtered_cytoplasm_dab)
  ks_statistic <- round(ks_test_result$statistic, 5)
  ks_pvalue <- format.pval(ks_test_result$p.value, digits = 2, eps = .Machine$double.eps)
  print(paste("After cutoff Kolmogorov Smirnov test calculated between Nucleus DAB OD mean and cytoplasm DAB OD mean is: D =", ks_statistic, "with the p-value of", ks_pvalue))
  
  #JA:Calculating ANOVA:
  filtered_tma_long <- removed_tma %>%
    pivot_longer(cols = c(`Nucleus: DAB OD mean`, `Cytoplasm: DAB OD mean`), 
                 names_to = "Variable", values_to = "Value")
  
  anova_result <- aov(Value ~Variable, filtered_tma_long)
  
  print("Summary of ANOVA results:")
  print(summary(anova_result))
  
  #JA: Compute the Euclidean distance
  euclidean_distance <- round(sqrt(sum((filtered_nucleus_dab - filtered_cytoplasm_dab)^2)),2)
  
  print(paste("After cutoff Euclidean_distance calculated between Nucleus DAB OD mean and cytoplasm DAB OD mean is:", euclidean_distance, sep = " "))
  
  
  a2 <- min(min(filtered_tma[[col1]]), min(filtered_tma[[col2]]))
  print(paste("minimum values between", col1, "and", col2, "when cutoff are defined", "is", a2, sep=" "))
  
  
  b2 <- max(max(filtered_tma[[col1]]), max(filtered_tma[[col2]]))
  print(paste("maximum value between", col1, "and", col2, "when cutoff are defined", "is", b2, sep=" "))
  
  print(paste(bin_size,"bins are created between range", a2, "and", b2))
  
  bin_break<- seq(a2, b2, length.out=bin_size_new)
  
  hist_col1<- hist(filtered_tma[[col1]], breaks = bin_break, plot = FALSE)
  hist_col2 <- hist(filtered_tma[[col2]], breaks = bin_break, plot = FALSE)
  
  hist_col1_df <- data.frame(
    midpoints = (head(hist_col1$breaks, -1) + tail(hist_col1$breaks, -1)) / 2,
    counts = hist_col1$counts,
    Variable = col1
  )
  
  hist_col2_df <- data.frame(
    midpoints = (head(hist_col2$breaks, -1) + tail(hist_col2$breaks, -1)) / 2,
    counts = hist_col2$counts,
    Variable = col2
  )
  
  #JA:#JA:define Normalized values 
  hist_col1_df<- hist_col1_df %>% mutate(norm_counts = counts/sum(counts))
  hist_col2_df<- hist_col2_df %>% mutate(norm_counts = counts/sum(counts))
  
  core_final<- cbind(hist_col1_df, hist_col2_df)
  
  #JA:#JA:saving graph input
  fwrite(core_final, paste("graphinput after cutoff", input, sep = " "))
  fwrite(filtered_tma_long, paste("geom density graphinput after cutoff", input, sep = " "))
  
  
  #JA:#JA:defining graph axis
  max_param_value <- max(max(filtered_tma[[col1]]), max(filtered_tma[[col2]]))
  graph_max_value <- round(max_param_value, 1)+0.1
  
  
  #JA:#JA:#JA:#JA:#JA:Plotting graphs with absolute frequency: 
  print(paste("Plotting the graph between", col1, "and", col2, "when only Nucleus and Cytoplasm DAB OD mean has positive values after using user provided cut-off values"))
  
  print("Printing probability density function with absolute frequency")
  two_histo_approximation <- ggplot() +
    geom_bar(data = hist_col1_df, 
             aes(x = midpoints, y = counts, fill = "col1"), 
             stat = "identity", color = "purple", alpha = 0.7) + 
    geom_bar(data = hist_col2_df, 
             aes(x= midpoints, y= counts, fill= "col2"), 
             stat= "identity", color = "darkgreen", alpha= 0.7) +
    geom_density(data= tma_long %>% filter (Variable == col1), 
                 aes(x = Value, y = ..density.. * sum(hist_col1_df$counts) * diff(bin_break)[1], color = "col1"), size = 1, kernel = "gaussian", bw = 0.09) +
    geom_density(data= tma_long %>% filter (Variable == col2), 
                 aes(x = Value, y = ..density.. * sum(hist_col2_df$counts) * diff(bin_break)[1], color = "col2"), size = 1, kernel = "gaussian", bw = 0.09) +
    coord_cartesian(xlim = c(0, graph_max_value)) +
    scale_x_continuous(breaks = seq(0, graph_max_value, by = 0.2)) +
    labs(title = graph_title_after_cutoff, subtitle = "Value distribution fitted by probability density function", x = "DAB OD mean", y = "Frequency") +
    scale_fill_manual(name = "Legend", values = c("col1" = "purple", "col2" = "darkgreen"), labels= c(col1, col2)) +
    scale_color_manual(name = "Legend", values = c("col1" = "purple", "col2" = "darkgreen"),labels= c(col1, col2)) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          legend.position = "top")
  
  ggsave(paste("probability density graph after cutoff", input, ".tiff", sep = " "), 
         plot = two_histo_approximation, width=11,height=10, dpi=600)
  
  
  
  print("Printing probability density function with Normalized frequency")
  two_histo_approximation_norm <- ggplot() +
    geom_bar(data = hist_col1_df, 
             aes(x = midpoints, y = norm_counts, fill = "col1"), 
             stat = "identity", color = "purple", alpha = 0.7) + 
    geom_bar(data = hist_col2_df, 
             aes(x= midpoints, y= norm_counts, fill= "col2"), 
             stat= "identity", color = "darkgreen", alpha= 0.7) +
    geom_density(data= tma_long %>% filter (Variable == col1), 
                 aes(x = Value, y = ..density.. * sum(hist_col1_df$norm_counts) * diff(bin_break)[1], color = "col1"), size = 1, kernel = "gaussian", bw = 0.09) +
    geom_density(data= tma_long %>% filter (Variable == col2), 
                 aes(x = Value, y = ..density.. * sum(hist_col2_df$norm_counts) * diff(bin_break)[1], color = "col2"), size = 1, kernel = "gaussian", bw = 0.09) +
    coord_cartesian(xlim = c(0, graph_max_value)) +
    scale_x_continuous(breaks = seq(0, graph_max_value, by = 0.2)) +
    labs(title = graph_title_after_cutoff, subtitle = "Value distribution fitted by Weibull probability function", x = "DAB OD mean", y = "Normalized frequency") +
    scale_fill_manual(name = "Legend", values = c("col1" = "purple", "col2" = "darkgreen"), labels= c(col1, col2)) +
    scale_color_manual(name = "Legend", values = c("col1" = "purple", "col2" = "darkgreen"),labels= c(col1, col2)) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          legend.position = "top")
  
  ggsave(paste("probability density graph Normalized after cutoff", input, ".tiff", sep = " "), 
         plot = two_histo_approximation_norm, width=11,height=10, dpi=600)
  
  
  
  print("printing weibull probability function with absolute frequency after user specified cut-off values")
  #JA:#JA:defining other parameters needed for plotting weibull plotting
  
  fit_weibull_col1 <- fitdist(filtered_tma[[col1]], "weibull")
  fit_weibull_col2 <- fitdist(filtered_tma[[col2]], "weibull")
  
  bin_size_new_col1 <- hist_col1_df$midpoints
  bin_size_new_col2 <- hist_col2_df$midpoints
  
  weibull_density_col1 <- dweibull(bin_size_new_col1, shape = fit_weibull_col1$estimate["shape"], scale = fit_weibull_col1$estimate["scale"])
  weibull_density_col2 <- dweibull(bin_size_new_col2, shape = fit_weibull_col2$estimate["shape"], scale = fit_weibull_col2$estimate["scale"])
  
  col1_area <- sum(hist_col1_df$counts * diff(bin_break)[1])
  col2_area <- sum(hist_col2_df$counts * diff(bin_break)[1])
  
  
  two_histo_weibull<- ggplot() +
    geom_bar(data = hist_col1_df, 
             aes(x = midpoints, y = counts, fill = "col1"), 
             stat = "identity", color = "purple", alpha = 0.7) + 
    geom_bar(data = hist_col2_df, 
             aes(x= midpoints, y= counts, fill= "col2"), 
             stat= "identity", color = "darkgreen", alpha= 0.7) +
    geom_line(data= hist_col1_df, 
              aes(x = bin_size_new_col1, y = weibull_density_col1*col1_area, color = "col1"), size = 1) + 
    geom_line(data= hist_col2_df, 
              aes(x = bin_size_new_col2, y = weibull_density_col2*col2_area, color = "col2"), size = 1) +
    coord_cartesian(xlim = c(0, graph_max_value)) +
    scale_x_continuous(breaks = seq(0, graph_max_value, by = 0.2)) +
    labs(title = graph_title_after_cutoff, subtitle = "Value distribution fitted by Weibull probability function", x = "DAB OD mean", y = "Frequency") +
    scale_fill_manual(name = "Legend", values = c("col1" = "purple", "col2" = "darkgreen"), labels= c(col1, col2)) +
    scale_color_manual(name = "Legend", values = c("col1" = "purple", "col2" = "darkgreen"),labels= c(col1, col2)) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          legend.position = "top")
  
  ggsave(paste("Weibull graph after cutoff", input, ".tiff", sep = " "), 
         plot = two_histo_weibull, width=11,height=10, dpi=600)
  
  
  print("printing weibull probability function with Normalized frequency")
  norm_col1_area <- sum(hist_col1_df$norm_counts * diff(bin_break)[1])
  norm_col2_area <- sum(hist_col2_df$norm_counts * diff(bin_break)[1])
  
  two_histo_weibull_norm<- ggplot() +
    geom_bar(data = hist_col1_df, 
             aes(x = midpoints, y = norm_counts, fill = "col1"), 
             stat = "identity", color = "purple", alpha = 0.7) + 
    geom_bar(data = hist_col2_df, 
             aes(x= midpoints, y= norm_counts, fill= "col2"), 
             stat= "identity", color = "darkgreen", alpha= 0.7) +
    geom_line(data= hist_col1_df, 
              aes(x = bin_size_new_col1, y = weibull_density_col1*norm_col1_area, color = "col1"), size = 1) + 
    geom_line(data= hist_col2_df, 
              aes(x = bin_size_new_col2, y = weibull_density_col2*norm_col2_area, color = "col2"), size = 1) +
    coord_cartesian(xlim = c(0, graph_max_value)) +
    scale_x_continuous(breaks = seq(0, graph_max_value, by = 0.2)) +
    labs(title = graph_title_after_cutoff, subtitle = "Value distribution fitted by Weibull probability function", x = "DAB OD mean", y = "Frequency") +
    scale_fill_manual(name = "Legend", values = c("col1" = "purple", "col2" = "darkgreen"), labels= c(col1, col2)) +
    scale_color_manual(name = "Legend", values = c("col1" = "purple", "col2" = "darkgreen"),labels= c(col1, col2)) +
    theme(plot.title = element_text(hjust = 0.5),
          plot.subtitle = element_text(hjust = 0.5),
          text = element_text(size = 30),
          axis.line = element_line(color = "black"),
          legend.position = "top")
  
  ggsave(paste("Weibull graph Normalized after cutoff", input, ".tiff", sep = " "), 
         plot = two_histo_weibull_norm, width=11,height=10, dpi=600)
  
  print("Printing complete check your output directory folder")
  
  print(paste("Analysis complete for", input))
} else{
  print("No execution of the next part of the code because Nucleus: Max caliper or Cell: DAB OD mean cutoffs are not provided by the user")
}

#JA:{ if curve needs to start from zero
#JA:max_col1_value<- round(max(filtered_tma[[col1]]), 1)
#JA:max_col2_value<- round(max(filtered_tma[[col2]]), 1)

#JA:min_parameter_value<- min(max_col1_value, max_col2_value)

#JA:bin_size_new_col1 <- seq_len(nrow(hist_col1_df))

#JA:if (max_col1_value > min_parameter_value){
#JA:bin_size_new_col1<- seq(from = 0, to = min_parameter_value, length.out = nrow(hist_col1_df))
#JA:} else{
#JA:bin_size_new_col1<- seq(from = 0, to = max_col1_value, length.out = nrow(hist_col1_df))
#JA:}

#JA:bin_size_new_col2 <- seq_len(nrow(hist_col2_df))

#JA:bin_size_new_col2<- seq(from = 0, to = max_col2_value, length.out = nrow(hist_col2_df))
#JA:if (max_col2_value > min_parameter_value){
#JA:bin_size_new_col2<- seq(from = 0, to = min_parameter_value, length.out = nrow(hist_col2_df))
#JA:} else{
#JA:bin_size_new_col2<- seq(from = 0, to = max_col2_value, length.out = nrow(hist_col2_df))
#JA:}







