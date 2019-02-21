

# 0i_Importing Ecological data ####

library(tidyverse); library(ggbiplot)

EltonTraits <- read.delim("Data/MamFuncDat.txt") %>% na.omit
EltonTraits$Scientific <- EltonTraits$Scientific %>% str_replace(" ","_")

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt")
Panth2 <- read.delim("data/PanTHERIA_1-0_WR93_Aug2008.txt")

Panth1$MSW05_Binomial <- 
  Panth1$MSW05_Binomial %>% str_replace(" ", "_")

EcoDistVars <- c(
  "X5.1_AdultBodyMass_g", 
  "X15.1_LitterSize", 
  "X16.1_LittersPerYear",
  "X22.1_HomeRange_km2",
  "X27.2_HuPopDen_Mean_n.km2",
  "X28.1_Precip_Mean_mm",
  "X28.2_Temp_Mean_01degC",
  "X6.2_TrophicLevel",
  "X10.2_SocialGrpSize",
  "X6.1_DietBreadth", 
  "X12.1_HabitatBreadth"
)

kt <- ktab.data.frame(Panth1[,EcoDistVars])

EltonTraits$Carnivore <- ifelse(rowSums(EltonTraits[,c("Diet.Vend","Diet.Vect","Diet.Vfish")])>10,1,0)



DietComp <- EltonTraits %>% select(starts_with("Diet")) %>% select(1:10) %>% as.tibble
Remove <- which(rowSums(DietComp)==0)
DietComp <- DietComp %>% slice(-Remove)

VD <- vegdist(DietComp) %>% as.matrix

colnames(VD) <- rownames(VD) <- EltonTraits %>% slice(-Remove) %>% select(Scientific) %>% unlist

LongDiet <- reshape2::melt(VD) %>%     
  rename(Sp = Var1, Sp2 = Var2, DietSim = value)

FinalHostMatrix <- FinalHostMatrix %>% left_join(LongDiet, by = c("Sp","Sp2"))

FinalHostMatrix <- FinalHostMatrix %>% merge(EltonTraits[,c("Scientific","Carnivore")], by.x = "Sp", by.y = "Scientific")
FinalHostMatrix <- FinalHostMatrix %>% merge(EltonTraits[,c("Scientific","Carnivore")], by.x = "Sp2", by.y = "Scientific",
                                       suffixes = c("",".Sp2"))

FinalHostMatrix$Eaten <- ifelse(FinalHostMatrix$Carnivore==FinalHostMatrix$Carnivore.Sp2,0,1)




