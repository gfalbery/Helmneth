
# Running Frequentist GAMS

library(mgcv); library(ggregplot); library(tidyverse)

BAMList <- DataList <- list()

for(r in 1:length(WormGroups)){
  
  print(r)
  
  DataList[[WormGroups[r]]] <- FinalHostMatrix %>% filter(!is.na(WormGroups[r]))
  
  Formula = as.formula(paste0(WormGroups[r],"Binary", "~ t2(Space, scale(Phylo)) + s(DietSim) + Eaten"))
  
  BAMList[[WormGroups[r]]] <- bam(Formula,# + Eaten, # + 
                             #(1 + mmc(dom, dom_2) + mmc(d_cites_s1, d_cites_s2) | mm(Sp, Sp2)),
                             data = DataList[[WormGroups[r]]], 
                             family = binomial())
  
}

# Frequentist GAM Output ####

SpaceRange <- seq(from = min(FinalHostMatrix$Space),
                  to = max(FinalHostMatrix$Space),
                  length = 100) %>%
  c(mean(FinalHostMatrix$Space))

PhyloRange <- seq(from = min(FinalHostMatrix$Phylo),
                  to = max(FinalHostMatrix$Phylo),
                  length = 100) %>%
  c(mean(FinalHostMatrix$Phylo))

DietRange <- seq(from = min(FinalHostMatrix$DietSim),
                 to = max(FinalHostMatrix$DietSim),
                 length = 100) %>%
  c(mean(FinalHostMatrix$DietSim))

GAMPredDF <- expand.grid(Space = SpaceRange,
                         Phylo = PhyloRange,
                         DietSim = DietRange,
                         Eaten = c(0,1)
)

GAMPredDF <- GAMPredDF %>% mutate(SpaceQ = cut(Space, quantile(Space, 0:10/10),include.lowest = T, labels = 1:10),
                                  PhyloQ = cut(Phylo, quantile(Phylo, 0:10/10),include.lowest = T, labels = 1:10))

GAMPredDF[,paste0(WormGroups,"Fit")] <- lapply(BAMList, function(a) logistic(predict(a, newdata = GAMPredDF))) %>% 
  bind_cols

FitSlopeTime <- gather(GAMPredDF, key = "Model", value = "Estimate", paste0(WormGroups,"Fit"))

ggplot(FitSlopeTime %>% filter(DietSim == 0), aes(Space, Phylo)) + 
  geom_tile(aes(fill = Estimate)) +
  geom_point(data = FinalHostMatrix) + 
  labs(y = "Predicted Value") +
  facet_wrap(Eaten ~ Model, labeller = labeller(c("Not Eaten" = "0", "Eaten" = "1")))

ggplot(FitSlopeTime %>% 
         filter(Space == dplyr::last(SpaceRange),
                Phylo == dplyr::last(PhyloRange)), aes(DietSim, Estimate)) + 
  geom_line(aes(colour = Eaten, group = as.factor(Eaten)), alpha = 0.5) +
  labs(y = "Predicted Value") +
  facet_wrap(~ Model, labeller = labeller(Eaten = c("1" = "Eaten", "0" = "Not Eaten")), ncol = 2) +
  ggsave("Figures/Helminth_Diet.jpeg", units = "mm", width = 150, height = 150)

ggplot(FitSlopeTime %>% 
         filter(DietSim == 0, Eaten == 1), aes(Space, Estimate)) + 
  geom_line(aes(colour = Phylo, group = as.factor(Phylo)), alpha = 0.5) +
  labs(y = "Predicted Value") +
  facet_wrap(~ Model, labeller = labeller(Eaten = c("1" = "Eaten", "0" = "Not Eaten")), ncol = 2) +
  ggsave("Figures/Helminth_Space.jpeg", units = "mm", width = 200, height = 200)

ggplot(FitSlopeTime %>% 
         filter(DietSim == 0, Eaten == 1), aes(Phylo, Estimate)) + 
  geom_line(aes(colour = Space, group = as.factor(Space)), alpha = 0.5) +
  labs(y = "Predicted Value") +
  facet_wrap(~ Model, labeller = labeller(Eaten = c("1" = "Eaten", "0" = "Not Eaten")), ncol = 2) +
  ggsave("Figures/Helminth_Phylo.jpeg", units = "mm", width = 200, height = 200)

FitSlopeTime %>%
  filter(DietSim = last(unique(DietSim)))

BarGraph(FitSlopeTime, "Model", "Estimate", "Eaten") +
  ggsave()
