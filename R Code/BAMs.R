
# Running Frequentist GAMS

library(mgcv); library(ggregplot); library(tidyverse)

BAMList <- DataList <- PPList <- list()

Covar <- c("s(Phylo, by = ordered(Gz))",
           "t2(Phylo, Space, by = ordered(!Gz))",
           "s(DietSim)",
           "Eaten")#,
#"MinCites",
#"Domestic",
#"Spp")

for(r in 1:length(WormGroups)){
  
  print(WormGroups[r])
  
  DataList[[WormGroups[r]]] <- FinalHostMatrix %>% filter(!is.na(WormGroups[r])) %>% droplevels
  
  DataList[[WormGroups[r]]]$Sp <- factor(DataList[[WormGroups[r]]]$Sp, levels = sort(union(DataList[[WormGroups[r]]]$Sp,DataList[[WormGroups[r]]]$Sp2)))
  DataList[[WormGroups[r]]]$Sp2 <- factor(DataList[[WormGroups[r]]]$Sp2, levels = sort(union(DataList[[WormGroups[r]]]$Sp,DataList[[WormGroups[r]]]$Sp2)))
  
  DataList[[WormGroups[r]]] <- DataList[[WormGroups[r]]] %>% slice(order(Sp, Sp2))
  
  #MZ1 <- model.matrix( ~ Sp - 1, data = DataList[[WormGroups[r]]]) %>% as.matrix
  #MZ2 <- model.matrix( ~ Sp2 - 1, data = DataList[[WormGroups[r]]]) %>% as.matrix
  
  #SppMatrix = MZ1 + MZ2
  
  #DataList[[Resps[[r]]]]$Spp <- SppMatrix
  
  #PPList[[WormGroups[r]]] <- list(Spp = list(rank = nlevels(DataList[[WormGroups[r]]]$Sp), 
  #                                           diag(nlevels(DataList[[WormGroups[r]]]$Sp))))
  
  Formula = as.formula(paste0(WormGroups[r],"Binary",
                              " ~ ",
                              paste(Covar, collapse = " + ")
  ))
  
  BAMList[[WormGroups[r]]] <- bam(Formula,
                                  data = DataList[[WormGroups[r]]], 
                                  family = binomial(),
                                  #paraPen = PPList[[WormGroups[r]]], 
                                  select = T)
  
}

# Frequentist GAM Output ####

SpaceRange <- seq(from = 0,
                  to = 1,
                  length = 101) %>%
  c(mean(FinalHostMatrix$Space))

PhyloRange <- seq(from = 0,
                  to = 1,
                  length = 101) %>%
  c(mean(FinalHostMatrix$Phylo))

DietRange <- seq(from = 0,
                 to = 1,
                 length = 101) %>%
  c(mean(FinalHostMatrix$DietSim))

GAMPredDF <- expand.grid(Space = SpaceRange,
                         Phylo = PhyloRange,
                         DietSim = DietRange,
                         Eaten = c(0,1)
) %>% mutate(Gz = as.numeric(Space==0))

GAMPredDF <- GAMPredDF %>% mutate(SpaceQ = cut(Space, quantile(Space, 0:10/10),include.lowest = T, labels = 1:10),
                                  PhyloQ = cut(Phylo, quantile(Phylo, 0:10/10),include.lowest = T, labels = 1:10))

GAMPredDF[,paste0(rep(WormGroups,each=2),"_",rep(c("Fit","SE"),4))] <- lapply(BAMList, function(a) predict.bam(a, newdata = GAMPredDF, se.fit = T)) %>% 
  bind_cols

FitSlopeTime <- gather(GAMPredDF, key = "Model", value = "Fit", 
                       ends_with("Fit")) %>%
  bind_cols(gather(GAMPredDF, key = "Model", value = "SE", 
                   ends_with("SE"))) %>%
  mutate(Lower = Fit-SE, Upper = Fit+SE) %>%
  mutate_at(c("Fit", "Upper", "Lower"), logistic)

FitSlopeTime %>% 
  filter(Space == dplyr::last(Space),
         Phylo == dplyr::last(Phylo)) %>% 
  ggplot(aes(DietSim, Fit)) + 
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = as.factor(Eaten)), colour = NA, alpha = 0.3) +
  geom_line(aes(colour = as.factor(Eaten), group = as.factor(Eaten))) +
  labs(y = "Predicted Value") +
  lims(x = c(0,1), y = c(0,1)) +
  facet_wrap(~ Model, ncol = 2) +
  scale_color_discrete_sequential(palette = "purp", nmax = 4, order = c(2,4))  +
  scale_fill_discrete_sequential(palette = "purp", nmax = 4, order = c(2,4))  +
  theme(strip.background = element_rect(fill = "white", colour = "dark grey")) +
  ggsave("Figures/Helminth_Diet.jpeg", units = "mm", width = 150, height = 150)

FitSlopeTime %>% 
  filter(DietSim == dplyr::last(DietSim), 
         Eaten == 1,
         Phylo %in% c(0,dplyr::last(PhyloRange), 0.25, 0.5)) %>% 
  ggplot(aes(Space, Fit, colour = as.factor(Phylo), fill = as.factor(Phylo))) + 
  geom_ribbon(aes(ymin = Lower, ymax = Upper), colour = NA, alpha = 0.3) +
  geom_line() +
  labs(y = "Predicted Value") +
  facet_wrap(~ Model, ncol = 2) +
  lims(x = c(0,1), y = c(0,1)) +
  scale_color_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = 5:8)  +
  scale_fill_discrete_sequential(palette = AlberPalettes[[2]], nmax = 8, order = 5:8)  +
  theme(strip.background = element_rect(fill = "white", colour = "dark grey")) +
  ggsave("Figures/Helminth_Space.jpeg", units = "mm", width = 200, height = 200)

FitSlopeTime %>% 
  filter(DietSim == dplyr::last(DietSim), 
         Eaten == 1,
         Space %in% c(0,dplyr::last(SpaceRange), 0.25, 0.5)) %>% 
  ggplot(aes(Phylo, Fit, colour = as.factor(Space), fill = as.factor(Space))) + 
  geom_ribbon(aes(ymin = Lower, ymax = Upper), colour = NA, alpha = 0.3) +
  geom_line() +
  labs(y = "Predicted Value") +
  facet_wrap(~ Model, ncol = 2) +
  scale_color_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = 5:8)  +
  scale_fill_discrete_sequential(palette = AlberPalettes[[1]], nmax = 8, order = 5:8)  +
  ggsave("Figures/Helminth_Phylo.jpeg", units = "mm", width = 200, height = 200)


