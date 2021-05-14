

library(readr)
library(DescTools)
library(ggplot2)
library(igraph)
library(tidyr)
library(tidyverse)

#entropy_function, after singular val transformed

entropy <- function (sv, bins, normalizeEntropy = TRUE){
  N = length(sv)
  svn = sv/max(sv)
  
  counts = .bincode(svn,breaks = seq(0,1,1/bins),include.lowest = TRUE)
  p = counts/sum(counts)
  
  entropy = sum(-p*log(p))
  
  if(normalizeEntropy){
    entropy = entropy/log(length(p))
  }
  
  return(entropy)
}

Gini <- function(sv){
  coef = Gini(temp)
  return(coef)
  
}


df = data.frame("power" = NULL, "entropy" = NULL, "bins" = NULL)
powers = c(1.0, 2.0, 3.0)
bins = 50


for (power in powers){
  print(paste0("power = ",power))
  for (ii in 1:50){
    g <- sample_pa(1000, m=5, power = power, directed = FALSE)
    mat = as_adjacency_matrix(g, sparse = FALSE)
    #mat = mat + t(mat)
    stopifnot(isSymmetric(mat))
    stopifnot(max(mat) == 1)
    evals = eigen(mat, only.values = TRUE)$values
    
    vect = evals
    vect[vect<1] = 1
    vect = log(vect)
    
    S = entropy(vect, bins, normalizeEntropy = TRUE)
    
    current = data.frame("power" = power, "entropy" =  S, "bins" = bins)
    
    df = rbind(df,current)
  }
}




#various value of gamma in BA model for n=50

df = data.frame("power" = NULL, "entropy" = NULL, "bins" = NULL)
powers = seq(0,2,0.1)
bins = 50
N = 50


for (power in powers){
  print(paste0("power = ",power))
  for (ii in 1:N){
    g <- sample_pa(1000, m=5, power = power, directed = FALSE)
    mat = as_adjacency_matrix(g, sparse = FALSE)
    #mat = mat + t(mat)
    stopifnot(isSymmetric(mat))
    stopifnot(max(mat) == 1)
    evals = eigen(mat, only.values = TRUE)$values
    
    vect = evals
    vect[vect<1] = 1
    vect = log(vect)
    
    S = entropy(vect, bins, normalizeEntropy = TRUE)
    
    current = data.frame("power" = power, "entropy" =  S, "bins" = bins)
    
    df = rbind(df,current)
  }
}

fileName = paste0("simulations/sims_BA_N_", N, ".csv")
write.csv(df, file=fileName, row.names = FALSE)


#same for other N


bin_vect = c(10,20, 100, 200)

for (bins in bin_vect){
  df = data.frame("power" = NULL, "entropy" = NULL, "bins" = NULL)
  powers = seq(0,2,0.1)
  #bins = 50
  N = 50
  
  
  for (power in powers){
    print(paste0("power = ",power))
    for (ii in 1:N){
      g <- sample_pa(1000, m=5, power = power, directed = FALSE)
      mat = as_adjacency_matrix(g, sparse = FALSE)
      #mat = mat + t(mat)
      stopifnot(isSymmetric(mat))
      stopifnot(max(mat) == 1)
      evals = eigen(mat, only.values = TRUE)$values
      
      vect = evals
      vect[vect<1] = 1
      vect = log(vect)
      
      S = entropy(vect, bins, normalizeEntropy = TRUE)
      
      current = data.frame("power" = power, "entropy" =  S, "bins" = bins)
      
      df = rbind(df,current)
    }
  }
  
  fileName = paste0("simulations/sims_BA_N_", N,"_bins_",bins, ".csv")
  write.csv(df, file=fileName, row.names = FALSE)
}


#visualize case
df_fromFile = as_tibble(read.csv("simulations/sims_BA_N_50.csv")) 

df_tidy_mean <- df_fromFile %>%
  filter(!is.na(power)) %>%
  group_by(power, bins) %>%
  summarise(n = n(),
            mean = mean(entropy),
            median = median(entropy),
            sd = sd(entropy)) %>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)

ggplot(df_tidy_mean, aes(x=power, y=mean, color = bins)) +
  geom_line() +
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=bins),color="grey70",alpha=0.4) +
  labs(y = "Spectral Entropy", x = 'Power', main = "Spectral Entropy of Albert-Barabási Model") +
  ggtitle("Spectral Entropy of Albert-Barabási Model") + theme(legend.position = "none")

ggsave("ABmodel_spectral_entropy.pdf")

df = data.frame("power" = NULL, "entropy" = NULL, "bins" = NULL)
powers = seq(0,2,0.1)
bins = 10
N = 50


for (power in powers){
  print(paste0("power = ",power))
  for (ii in 1:N){
    g <- sample_pa(200, m=1, power = power, directed = FALSE)
    mat = as_adjacency_matrix(g, sparse = FALSE)
    #mat = mat + t(mat)
    stopifnot(isSymmetric(mat))
    stopifnot(max(mat) == 1)
    evals = eigen(mat, only.values = TRUE)$values
    
    vect = evals
    vect[vect<1] = 1
    vect = log(vect)
    
    S = entropy(vect, bins, normalizeEntropy = TRUE)
    
    current = data.frame("power" = power, "entropy" =  S, "bins" = bins)
    
    df = rbind(df,current)
  }
}

df_tidy_mean <- df %>%
  filter(!is.na(power)) %>%
  group_by(power, bins) %>%
  summarise(n = n(),
            mean = mean(entropy),
            median = median(entropy),
            sd = sd(entropy)) %>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)

ggplot(df_tidy_mean, aes(x=power, y=mean, color = bins)) +
  geom_line() +
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=bins),color="grey70",alpha=0.4) +
  labs(y = "Spectral Entropy", x = 'Power', main = "Spectral Entropy of Albert-Barabási Model") +
  ggtitle("Spectral Entropy of Albert-Barabási Model") + theme(legend.position = "none")


df_fromFile = as_tibble(read.csv("simulations/sims_BA_N_50.csv")) 

df_tidy_mean <- df_fromFile %>%
  filter(!is.na(power)) %>%
  group_by(power, bins) %>%
  summarise(n = n(),
            mean = mean(entropy),
            median = median(entropy),
            sd = sd(entropy)) %>%
  mutate(sem = sd / sqrt(n - 1),
         CI_lower = mean + qt((1-0.95)/2, n - 1) * sem,
         CI_upper = mean - qt((1-0.95)/2, n - 1) * sem)

ggplot(df_tidy_mean, aes(x=power, y=mean, color = bins)) +
  geom_line() +
  geom_ribbon(aes(ymin=CI_lower,ymax=CI_upper,fill=bins),color="grey70",alpha=0.4) +
  labs(y = "Spectral Entropy", x = expression('PA Power'~alpha)) +theme_bw()
#+ theme(legend.position = "none")+  ggtitle("Spectral Entropy of Nonlinear PA Model") 
ggsave("PAmodel_spectral_entropy.pdf")

df_ASE_read = read.csv("simulations/dataframe_for_PAplots_df_ASE.csv")
df_ASE_noAbs_read = read.csv("simulations/dataframe_for_PAplots_df_ASE_noAbsolute.csv")
df_overall_read = read.csv("simulations/dataframe_for_PAplots_df_overall.csv")

ggplot(data= df_overall_read, aes( x = power, y = entropy,
                                   color = factor(power),
                                   shape = factor(entropy_label),
                                   size = 20)) +
  geom_point() +
  scale_color_discrete() + guides(size = FALSE) +
  labs(x ="",y = "Spectral Entropy", color = expression("Power"~alpha), shape = "Spectral\n Entropy\n Rank")+theme_bw()

#ggsave("plot_PA_ASE_top.pdf")

ggplot(data = df_ASE_read, aes(x = X, y = Y,
                               color = factor(power),
                               shape = factor(entropy_label))) +
  geom_point(size = 3, alpha = .5) + facet_grid(factor(entropy_label,levels = c(50,25,1))~power) +
  scale_color_discrete() + guides(size=FALSE, alpha = FALSE) +
  labs(x =expression(abs(X[1])),y = expression(abs(X[2])), color = expression("Power"~alpha), shape = "Spectral\n Entropy\n Rank")+theme_bw()

#ggsave("plot_PA_ASE_bottom.pdf")

ggplot(data = df_ASE_noAbs_read, aes(x = X, y = Y,
                                     color = factor(power),
                                     shape = factor(entropy_label))) +
  geom_point(size = 3, alpha = .5) + facet_grid(factor(entropy_label,levels = c(50,25,1))~power) +
  scale_color_discrete() + guides(size=FALSE, alpha = FALSE) +
  labs(x =expression(X[1]),y = expression(X[2]), color = expression("Power"~alpha), shape = "Spectral\n Entropy\n Rank")+theme_bw()

#ggsave("plot_PA_ASE_noAbsolute.pdf")