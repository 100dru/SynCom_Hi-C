setwd("C:/Users/shatadru.2/OneDrive - The Ohio State University/LAB/Hi-C/informatic_pipeline/6th_test/metacc")

library(tidyverse)
library(dplyr)
library(ggplot2)
library(readxl)
library(readr)
library(ggpubr)
library(RColorBrewer)
library(viridisLite)
library(gplots)

all_phage_host <- read_csv("all_mock_conc_zscore_sens_spec.csv")


# Plots


### subsetting all the mocks

mock_1B1 <- all_phage_host %>% filter(Mock == '1-B1')
mock_1B3 <- all_phage_host %>% filter(Mock == '1-B3')
mock_1D1 <- all_phage_host %>% filter(Mock == '1-D1')
mock_1D2 <- all_phage_host %>% filter(Mock == '1-D2')
mock_1D3 <- all_phage_host %>% filter(Mock == '1-D3')
mock_1N1 <- all_phage_host %>% filter(Mock == '1-N1')
mock_1N2 <- all_phage_host %>% filter(Mock == '1-N2')
mock_1N3 <- all_phage_host %>% filter(Mock == '1-N3')
mock_2N1 <- all_phage_host %>% filter(Mock == '2_N1_S1')
mock_2N2 <- all_phage_host %>% filter(Mock == '2_N2_S2')
mock_2N3 <- all_phage_host %>% filter(Mock == '2_N3_S3')
mock_3N1 <- all_phage_host %>% filter(Mock == '3_N1_S4')
mock_3N3 <- all_phage_host %>% filter(Mock == '3_N3_S6')


#heatmaps

mock_control <- read.csv("mock_control.csv")

##heatmaps without threshold
mock_control$Presence <- as.factor(mock_control$Presence)

mock_1B1$Presence <- ifelse(!is.na(mock_1B1$Score), 1, 0)
mock_1B1$Presence <- as.factor(mock_1B1$Presence)
mock_1B3$Presence <- ifelse(!is.na(mock_1B3$Score), 1, 0)
mock_1B3$Presence <- as.factor(mock_1B3$Presence)
mock_1D1$Presence <- ifelse(!is.na(mock_1D1$Score), 1, 0)
mock_1D1$Presence <- as.factor(mock_1D1$Presence)
mock_1D2$Presence <- ifelse(!is.na(mock_1D2$Score), 1, 0)
mock_1D2$Presence <- as.factor(mock_1D2$Presence)
mock_1D3$Presence <- ifelse(!is.na(mock_1D3$Score), 1, 0)
mock_1D3$Presence <- as.factor(mock_1D3$Presence)
mock_1N1$Presence <- ifelse(!is.na(mock_1N1$Score), 1, 0)
mock_1N1$Presence <- as.factor(mock_1N1$Presence)
mock_1N2$Presence <- ifelse(!is.na(mock_1N2$Score), 1, 0)
mock_1N2$Presence <- as.factor(mock_1N2$Presence)
mock_1N3$Presence <- ifelse(!is.na(mock_1N3$Score), 1, 0)
mock_1N3$Presence <- as.factor(mock_1N3$Presence)
mock_2N1$Presence <- ifelse(!is.na(mock_2N1$Score), 1, 0)
mock_2N1$Presence <- as.factor(mock_2N1$Presence)
mock_2N2$Presence <- ifelse(!is.na(mock_2N2$Score), 1, 0)
mock_2N2$Presence <- as.factor(mock_2N2$Presence)
mock_2N3$Presence <- ifelse(!is.na(mock_2N3$Score), 1, 0)
mock_2N3$Presence <- as.factor(mock_2N3$Presence)
mock_3N1$Presence <- ifelse(!is.na(mock_3N1$Score), 1, 0)
mock_3N1$Presence <- as.factor(mock_3N1$Presence)
mock_3N3$Presence <- ifelse(!is.na(mock_3N3$Score), 1, 0)
mock_3N3$Presence <- as.factor(mock_3N3$Presence)


# List of datasets
datasets <- list(mock_1B1, mock_1B3, mock_1D1, mock_1D2, mock_1D3, mock_1N1, mock_1N2, mock_1N3, mock_2N1, mock_2N2, mock_2N3, mock_3N1, mock_3N3)

# Function to create the TileColor column and plot
create_plot <- function(data, name) {
  data$TileColor <- with(data, 
                         ifelse(Presence == 0, "white",
                                ifelse(Presence == 1 & Prediction == "Correct", "black", 
                                       ifelse(Presence == 1 & Prediction == "Incorrect", "gray", NA))))
  
  plot <- ggplot(data, aes(y = contig1_name, x = contig2_name)) +
    geom_tile(aes(fill = TileColor), colour = "white") +
    scale_fill_identity() +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  
  return(plot)
}

# Apply function to all datasets
plots <- lapply(seq_along(datasets), function(i) {
  dataset_name <- paste0("mock_", c("1B1", "1B3", "1D1", "1D2", "1D3", "1N1", "1N2", "1N3", "2N1", "2N2", "2N3", "3N1", "3N3")[i])
  create_plot(datasets[[i]], dataset_name)
})

# Extract individual plots if needed
tile_1B1 <- plots[[1]]
tile_1B3 <- plots[[2]]
tile_1D1 <- plots[[3]]
tile_1D2 <- plots[[4]]
tile_1D3 <- plots[[5]]
tile_1N1 <- plots[[6]]
tile_1N2 <- plots[[7]]
tile_1N3 <- plots[[8]]
tile_2N1 <- plots[[9]]
tile_2N2 <- plots[[10]]
tile_2N3 <- plots[[11]]
tile_3N1 <- plots[[12]]
tile_3N3 <- plots[[13]]
# Display a plot as an example
tile_3N3

tile_control <- ggplot(mock_control, aes(y = contig1_name, x = contig2_name)) +
  geom_tile(aes(fill = Presence), colour = "white") +
  scale_fill_manual(values = c("white", "black"))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

tile_control
ggsave("tile_control.pdf", tile_control)

# Convert Prediction to numeric if it is categorical
mock_1N1$Prediction_numeric <- as.numeric(as.factor(mock_1N1$Prediction))
mock_1N2$Prediction_numeric <- as.numeric(as.factor(mock_1N2$Prediction))
mock_1N3$Prediction_numeric <- as.numeric(as.factor(mock_1N3$Prediction))
mock_1B1$Prediction_numeric <- as.numeric(as.factor(mock_1B1$Prediction))
mock_1B3$Prediction_numeric <- as.numeric(as.factor(mock_1B3$Prediction))
mock_1D1$Prediction_numeric <- as.numeric(as.factor(mock_1D1$Prediction))
mock_1D2$Prediction_numeric <- as.numeric(as.factor(mock_1D2$Prediction))
mock_1D3$Prediction_numeric <- as.numeric(as.factor(mock_1D3$Prediction))
mock_2N1$Prediction_numeric <- as.numeric(as.factor(mock_2N1$Prediction))
mock_2N2$Prediction_numeric <- as.numeric(as.factor(mock_2N2$Prediction))
mock_2N3$Prediction_numeric <- as.numeric(as.factor(mock_2N3$Prediction))
mock_3N1$Prediction_numeric <- as.numeric(as.factor(mock_3N1$Prediction))
mock_3N3$Prediction_numeric <- as.numeric(as.factor(mock_3N3$Prediction))

##heatmaps with threshold at Z_Score > 1

mock_1N1$aboveZThresh <- ifelse(mock_1N1$z_score >= 0.5, 1, 0)
mock_1N1$aboveZThresh <- as.factor(mock_1N1$aboveZThresh)

mock_1N2$aboveZThresh <- ifelse(mock_1N2$z_score >= 0.5, 1, 0)
mock_1N2$aboveZThresh <- as.factor(mock_1N2$aboveZThresh)

mock_1N3$aboveZThresh <- ifelse(mock_1N3$z_score >= 0.5, 1, 0)
mock_1N3$aboveZThresh <- as.factor(mock_1N3$aboveZThresh)

mock_1B1$aboveZThresh <- ifelse(mock_1B1$z_score >= 0.5, 1, 0)
mock_1B1$aboveZThresh <- as.factor(mock_1B1$aboveZThresh)

mock_1B3$aboveZThresh <- ifelse(mock_1B3$z_score >= 0.5, 1, 0)
mock_1B3$aboveZThresh <- as.factor(mock_1B3$aboveZThresh)

mock_1D1$aboveZThresh <- ifelse(mock_1D1$z_score >= 0.5, 1, 0)
mock_1D1$aboveZThresh <- as.factor(mock_1D1$aboveZThresh)

mock_1D2$aboveZThresh <- ifelse(mock_1D2$z_score >= 0.5, 1, 0)
mock_1D2$aboveZThresh <- as.factor(mock_1D2$aboveZThresh)

mock_1D3$aboveZThresh <- ifelse(mock_1D3$z_score >= 0.5, 1, 0)
mock_1D3$aboveZThresh <- as.factor(mock_1D3$aboveZThresh)

mock_2N1$aboveZThresh <- ifelse(mock_2N1$z_score >= 0.5, 1, 0)
mock_2N1$aboveZThresh <- as.factor(mock_2N1$aboveZThresh)

mock_2N2$aboveZThresh <- ifelse(mock_2N2$z_score >= 0.5, 1, 0)
mock_2N2$aboveZThresh <- as.factor(mock_2N2$aboveZThresh)

mock_2N3$aboveZThresh <- ifelse(mock_2N3$z_score >= 0.5, 1, 0)
mock_2N3$aboveZThresh <- as.factor(mock_2N3$aboveZThresh)

mock_3N1$aboveZThresh <- ifelse(mock_3N1$z_score >= 0.5, 1, 0)
mock_3N1$aboveZThresh <- as.factor(mock_3N1$aboveZThresh)

mock_3N3$aboveZThresh <- ifelse(mock_3N3$z_score >= 0.5, 1, 0)
mock_3N3$aboveZThresh <- as.factor(mock_3N3$aboveZThresh)


# List of datasets
datasets <- list(mock_1B1, mock_1B3, mock_1D1, mock_1D2, mock_1D3, mock_1N1, mock_1N2, mock_1N3, mock_2N1, mock_2N2, mock_2N3, mock_3N1, mock_3N3)

# Function to create the TileColor column and plot
create_plot <- function(data, name) {
  data$TileColor <- with(data, 
                         ifelse(Presence == 0, "white",
                                ifelse(Presence == 1 & Prediction == "Correct" & aboveZThresh == 1, "black", 
                                       ifelse(Presence == 1 & Prediction == "Incorrect" & aboveZThresh == 1, "gray", NA))))
  
  plot <- ggplot(data, aes(y = contig1_name, x = contig2_name)) +
    geom_tile(aes(fill = TileColor), colour = "white") +
    scale_fill_identity() +
    theme_classic() +
    theme(axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
  return(plot)
}

# Apply function to all datasets
plots <- lapply(seq_along(datasets), function(i) {
  dataset_name <- paste0("mock_", c("1B1", "1B3", "1D1", "1D2", "1D3", "1N1", "1N2", "1N3", "2N1", "2N2", "2N3", "3N1", "3N3")[i])
  create_plot(datasets[[i]], dataset_name)
})

# Extract individual plots if needed

z_tile_1B1 <- plots[[1]]
z_tile_1B3 <- plots[[2]]
z_tile_1D1 <- plots[[3]]
z_tile_1D2 <- plots[[4]]
z_tile_1D3 <- plots[[5]]
z_tile_1N1 <- plots[[6]]
z_tile_1N2 <- plots[[7]]
z_tile_1N3 <- plots[[8]]
z_tile_2N1 <- plots[[9]]
z_tile_2N2 <- plots[[10]]
z_tile_2N3 <- plots[[11]]
z_tile_3N1 <- plots[[12]]
z_tile_3N3 <- plots[[13]]

# Display a plot as an example
z_tile_1N1
z_tile_1N2
z_tile_1N3
#ggsave("threshz_tile_all.png", plot = threshz_tile_all, width = 16, height = 9, units = "in", dpi = 300)

#ggsave("threshz_tile_all.pdf", plot = threshz_tile_all, width = 16, height = 9, units = "in", dpi = 300)
#ggsave("thresh_tile_all.png", plot = thresh_tile_all, width = 16, height = 9, units = "in", dpi = 300)


#Dot Plots for the mocks
create_dot_plot <- function(data) {
  ggplot(data=data, aes(x=Score, y=contig1_name, color=Prediction)) +
    scale_colour_manual(values=setNames(c('black', 'gray'), c("Correct", "Incorrect"))) +
    geom_point(position=position_jitter(width=0, height=0.2), size=4) +
    labs(x="Contact Score", 
         y="Hosts", 
         color="Prediction", 
         size="Mock Community") +
    theme_classic() +
    theme(
          axis.text.y = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
}


dot_1B1 <- create_dot_plot(mock_1B1)
dot_1B3 <- create_dot_plot(mock_1B3)
dot_1D1 <- create_dot_plot(mock_1D1)
dot_1D2 <- create_dot_plot(mock_1D2)
dot_1D3 <- create_dot_plot(mock_1D3)
dot_1N1 <- create_dot_plot(mock_1N1)
dot_1N2 <- create_dot_plot(mock_1N2)
dot_1N3 <- create_dot_plot(mock_1N3)
dot_2N1 <- create_dot_plot(mock_2N1)
dot_2N2 <- create_dot_plot(mock_2N2)
dot_2N3 <- create_dot_plot(mock_2N3)
dot_3N1 <- create_dot_plot(mock_3N1)
dot_3N3 <- create_dot_plot(mock_3N3)


### plot score dot plot against zscore

create_z_dot_plot <- function(data) {
  ggplot(data = data, aes(x = z_score, y = contig1_name, color = Prediction)) +
    scale_colour_manual(values = setNames(c('black', 'gray'), c("Correct", "Incorrect"))) +
geom_point(position = position_jitter(width = 0, height = 0.2), size = 4) +
  #geom_vline(xintercept = quantile(data$Score, 0.8), color = "azure4") +
  xlim(-1., 4.8) +
  labs(x = "Z Score", 
       y = "Hosts", 
       color = "Prediction", 
       size = "Mock Community") +
  theme_classic() +
    theme(
          axis.text.y = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank())
}

# Example usage:
z_dot_1B1 <- create_z_dot_plot(mock_1B1)
z_dot_1B3 <- create_z_dot_plot(mock_1B3)
z_dot_1D1 <- create_z_dot_plot(mock_1D1)
z_dot_1D2 <- create_z_dot_plot(mock_1D2)
z_dot_1D3 <- create_z_dot_plot(mock_1D3)
z_dot_1N1 <- create_z_dot_plot(mock_1N1)
z_dot_1N2 <- create_z_dot_plot(mock_1N2)
z_dot_1N3 <- create_z_dot_plot(mock_1N3)
z_dot_2N1 <- create_z_dot_plot(mock_2N1)
z_dot_2N2 <- create_z_dot_plot(mock_2N2)
z_dot_2N3 <- create_z_dot_plot(mock_2N3)
z_dot_3N1 <- create_z_dot_plot(mock_3N1)
z_dot_3N3 <- create_z_dot_plot(mock_3N3)

z_dot_1B1
z_dot_1N3

### figure paper arrangement
Fig_2_before_binary <- ggarrange(tile_1N1, tile_1N2, tile_1N3, tile_2N1, tile_2N2, tile_2N3, tile_3N1, tile_3N3,
                              labels = c("Mock 1.1", "Mock 1.2", "Mock 1.3", "Mock 2.1", "Mock 2.2", "Mock 2.3", "Mock 3.1", "Mock 3.2"),  
                              common.legend = TRUE, legend = "bottom",
                              ncol = 4, nrow = 2)


Fig_2_before_binary

Fig_2_before_dot <- ggarrange(dot_1N1, dot_1N2, dot_1N3, dot_2N1, dot_2N2, dot_2N3, dot_3N1, dot_3N3,
                                 labels = c("Mock 1.1", "Mock 1.2", "Mock 1.3", "Mock 2.1", "Mock 2.2", "Mock 2.3", "Mock 3.1", "Mock 3.2"),  
                                 common.legend = TRUE, legend = "bottom",
                                 ncol = 4, nrow = 2)


Fig_2_before_dot

Figure_2_before <- ggarrange(Fig_2_before_binary, Fig_2_before_dot, ncol = 1, nrow = 2)
Figure_2_before


## version 2

Fig_2_mock1_before <- ggarrange(tile_1N1,dot_1N1, tile_1N2,  dot_1N2, tile_1N3, dot_1N3,
                         labels = c("Mock 1.1", "", "Mock 1.2", "", "Mock 1.3", ""),  
                         common.legend = TRUE, legend = "bottom",
                         ncol = 6, nrow = 1)

Fig_2_mock1_after <- ggarrange(z_tile_1N1,z_dot_1N1, z_tile_1N2,  z_dot_1N2, z_tile_1N3, z_dot_1N3, 
                               labels = c("Mock 1.1", "", "Mock 1.2", "", "Mock 1.3", ""),  
                               common.legend = TRUE, legend = "bottom",
                               ncol = 6, nrow = 1)

Fig_2_mock1 <- ggarrange(Fig_2_mock1_before, Fig_2_mock1_after,common.legend = TRUE, legend = "bottom", ncol = 1, nrow = 2)
Fig_2_mock1

ggsave("C:/Users/shatadru.2/OneDrive - The Ohio State University/LAB/Hi-C/paper_unit/Figures/Figure_2/Figure_2_Mock_1.pdf", 
       plot = Fig_2_mock1, width = 16, height = 9, units = "in", dpi = 300)


#ggsave("C:/Users/shatadru.2/OneDrive - The Ohio State University/LAB/Hi-C/paper_unit/Figures/Figure_2/Figure_2_before.pdf", 
       plot = Figure_2_before, width = 16, height = 9, units = "in", dpi = 300)

Fig_2_after_binary <-  ggarrange(z_tile_1N1, z_tile_1N2, z_tile_1N3, z_tile_2N1, z_tile_2N2, z_tile_2N3, z_tile_3N1, z_tile_3N3,
                                 labels = c("Mock 1.1", "Mock 1.2", "Mock 1.3", "Mock 2.1", "Mock 2.2", "Mock 2.3", "Mock 3.1", "Mock 3.2"),  
                                 common.legend = TRUE, legend = "bottom",
                                 ncol = 4, nrow = 2)
Fig_2_after_binary




Fig_2_after_dot <-  ggarrange(z_dot_1N1, z_dot_1N2, z_dot_1N3, z_dot_2N1, z_dot_2N2, z_dot_2N3, z_dot_3N1, z_dot_3N3,
                                 labels = c("Mock 1.1", "Mock 1.2", "Mock 1.3", "Mock 2.1", "Mock 2.2", "Mock 2.3", "Mock 3.1", "Mock 3.2"),  
                                 common.legend = TRUE, legend = "bottom",
                                 ncol = 4, nrow = 2)
Fig_2_after_dot

Figure_2_after <- ggarrange(Fig_2_after_binary, Fig_2_after_dot, ncol = 1, nrow = 2)
Figure_2_after

#ggsave("C:/Users/shatadru.2/OneDrive - The Ohio State University/LAB/Hi-C/paper_unit/Figures/Figure_2/Figure_2_after.pdf", 
      # plot = Figure_2_after, width = 16, height = 9, units = "in", dpi = 300)

## supplementary for the preservative tests
sup_before_binary <- ggarrange(tile_1B1, tile_1B3, tile_1D1, tile_1D2, tile_1D3,
                                 labels = c("Mock B.1", "Mock B.2", "Mock D.1", "Mock D.2", "Mock D.3"),  
                                 common.legend = TRUE, legend = "bottom",
                                 ncol = 3, nrow = 2)


sup_before_binary

sup_before_dot <- ggarrange(dot_1B1, dot_1B3, dot_1D1, dot_1D2, dot_1D3,
                              labels = c("Mock B.1", "Mock B.2", "Mock D.1", "Mock D.2", "Mock D.3"),  
                              common.legend = TRUE, legend = "bottom",
                              ncol = 3, nrow = 2)


sup_before_dot

sup_before <- ggarrange(sup_before_binary, sup_before_dot, ncol = 1, nrow = 2)
sup_before




#ggsave("C:/Users/shatadru.2/OneDrive - The Ohio State University/LAB/Hi-C/paper_unit/Figures/Figure_2/Sup_4_before.pdf", 
   #    plot = sup_before, width = 16, height = 9, units = "in", dpi = 300)

sup_after_binary <-  ggarrange(z_tile_1B1, z_tile_1B3, z_tile_1D1, z_tile_1D2, z_tile_1D3,
                                 labels = c("Mock B.1", "Mock B.2", "Mock D.1", "Mock D.2", "Mock D.3"),  
                                 common.legend = TRUE, legend = "bottom",
                                 ncol = 3, nrow = 2)
sup_after_binary




sup_after_dot <-  ggarrange(z_dot_1B1, z_dot_1B3, z_dot_1D1, z_dot_1D2, z_dot_1D3,
                              labels = c("Mock B.1", "Mock B.2", "Mock D.1", "Mock D.2", "Mock D.3"),  
                              common.legend = TRUE, legend = "bottom",
                              ncol = 3, nrow = 2)
sup_after_dot

sup_after <- ggarrange(sup_after_binary, sup_after_dot, ncol = 1, nrow = 2)
sup_after

#ggsave("C:/Users/shatadru.2/OneDrive - The Ohio State University/LAB/Hi-C/paper_unit/Figures/Figure_2/Sup_4_after.pdf", 
       plot = sup_after, width = 16, height = 9, units = "in", dpi = 300)

##vrsion 2

Fig_3_preservatives <- ggarrange( z_tile_1D1, z_dot_1D1, z_tile_1D2, z_dot_1D2, z_tile_1D3, z_dot_1D3, 
                                  z_tile_1B1, z_dot_1B1, z_tile_1B3,  z_dot_1B3, 
                                 labels = c("Mock D.1", "", "Mock D.2","", "Mock D.3","", "Mock B.1","", "Mock B.2","" ),  
                                 common.legend = TRUE, legend = "bottom",
                                 ncol = 6, nrow = 2)
  
Fig_3_preservatives

#ggsave("C:/Users/shatadru.2/OneDrive - The Ohio State University/LAB/Hi-C/paper_unit/Figures/Figure_3/Fig_3_preservatives.pdf", 
     #  plot = Fig_3_preservatives, width = 16, height = 9, units = "in", dpi = 300)



n1_phage_score <- ggplot(mock_1N1,aes(x=contig2_name, y=Score, color=Prediction))+
  geom_point()+
  scale_colour_manual(values = setNames(c('black', 'gray'), c("Correct", "Incorrect"))) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

n2_phage_score <- ggplot(mock_1N2,aes(x=contig2_name, y=Score, color=Prediction))+
  geom_point()+
  scale_colour_manual(values = setNames(c('black', 'gray'), c("Correct", "Incorrect"))) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

n3_phage_score <- ggplot(mock_1N3,aes(x=contig2_name, y=Score, color=Prediction))+
  geom_point()+
  scale_colour_manual(values = setNames(c('black', 'gray'), c("Correct", "Incorrect"))) +
  theme_classic()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

n_phage_score <- ggarrange(n1_phage_score, n2_phage_score, n3_phage_score, 
                           labels = c("Mock 1.1", "Mock 1.2", "Mock 1.3"), common.legend = TRUE)

n_phage_score


###### figure 4

Fig_4_limit <- ggarrange( z_tile_2N1, z_dot_2N1, z_tile_2N2, z_dot_2N2, z_tile_2N3, z_dot_2N3, 
                                  z_tile_3N1, z_dot_3N1, z_tile_3N3,  z_dot_3N3, 
                                  labels = c("Mock 2.1", "", "Mock 2.2","", "Mock 2.3","", "Mock 3.1","", "Mock 3.2","" ),  
                                  common.legend = TRUE, legend = "bottom",
                                  ncol = 6, nrow = 2)

ggsave("C:/Users/shatadru.2/OneDrive - The Ohio State University/LAB/Hi-C/paper_unit/Figures/Figure_4/Fig_4_limit.pdf", 
  plot = Fig_4_limit, width = 16, height = 9, units = "in", dpi = 300)
