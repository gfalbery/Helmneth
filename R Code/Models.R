
# GAM SubModels ####

# Rscript "~/Albersnet/R Code/1_Sharing Models/5z_BRMS GAM SubModels.R" # This is the terminal run code

source("R Code/00_Master Code.R")

library(rstan); library(tidyverse); library(reskew); library(brms)

SubResps <- WormGroups

# Binomial model for viral sharing of viral subtypes ####

SubDataList <- StanDataList <- SubGAMModelList <- list()

# Import data

for(r in 1:length(SubResps)){
  
  SubDataList[[r]] <- FinalHostMatrix[!is.na(FinalHostMatrix[,SubResps[r]]),]
  
  # Generate species-level trait data
  # Get Sp and Sp2 in "d" on the same factor levels
  
  #SubDataList[[r]]$Sp <- factor(as.character(SubDataList[[r]]$Sp),
  #                              levels = union(SubDataList[[r]]$Sp, SubDataList[[r]]$Sp2)
  #)
  
  #SubDataList[[r]]$Sp2 <- factor(as.character(SubDataList[[r]]$Sp2),
  #                               levels = union(SubDataList[[r]]$Sp, SubDataList[[r]]$Sp2)
  #)
  
  SubDataList[[r]]$Sharing <- SubDataList[[r]][,SubResps[r]]
  
  SubModel <- brm(Sharing ~ t2(Space, Phylo) + s(DietSim) + Eaten,
                    #(1 + mmc(domestic, domestic.Sp2) + mmc(d_cites_s, d_cites_s2) | mm(Sp, Sp2)),
                  data = SubDataList[[r]], 
                  family = bernoulli(), 
                  cores = 8, 
                  seed = 17,
                  iter = 1500, 
                  warmup = 500, 
                  thin = 10, 
                  refresh = 0)
  
  SubGAMModelList[[r]] <- SubModel
  
  saveRDS(SubModel, 
          file = paste0(Resps[r],"SubModel.rds"))
  
  remove(SubModel)
  
}

saveRDS(SubGAMModelList, file = "Output Files/SubModelList.rds")
saveRDS(SubDataList, file = "Output Files/SubDataList.rds")

