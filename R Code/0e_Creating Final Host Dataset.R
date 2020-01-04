# Creating final dataset

library(tidyverse); library(vegan)

Panth1 <- read.delim("data/PanTHERIA_1-0_WR05_Aug2008.txt") %>%
  dplyr::rename(Sp = MSW05_Binomial, hOrder = MSW05_Order)
Panth1$Sp <- Panth1$Sp %>% str_replace(" ", "_")

NonEutherians <- c("Diprotodontia",
                   "Dasyuromorphia",
                   "Paucituberculata",
                   "Didelphimorphia",
                   "Microbiotheria",
                   "Peramelemorphia", 
                   "Notoryctemorphia",
                   "Monotremata")

NonEutherianSp <- Panth1[Panth1$hOrder%in%NonEutherians,"Sp"]

if(file.exists("Data/FullRangeOverlap.Rdata")) load("Data/FullRangeOverlap.Rdata") else{
  source(paste0(CodeRoot,"/","0_Data Import/0c2_Exhaustive Spatial Data Import.R"))
}

if(!file.exists("Data/intermediate/FullSTMatrix.csv")){
  
  library(geiger);library(ape);library(picante);library(dplyr)
  
  STFull <- read.nexus("data/ele_1307_sm_sa1.tre")[[1]]
  FullSTMatrix <- as.data.frame(cophenetic(STFull)) %>% as.matrix
  
} else{ FullSTMatrix <- as.matrix(read.csv("data/intermediate/FullSTMatrix.csv", header = T)) }

EutherianSp <- colnames(FullSTMatrix) %>% setdiff(NonEutherianSp)

FullSTMatrix <- FullSTMatrix[EutherianSp, EutherianSp]

rownames(Hosts) = Hosts$Sp

tFullSTMatrix <- 1 - (FullSTMatrix - min(FullSTMatrix))/max(FullSTMatrix)

MammalHostNames <- reduce(list(as.character(Hosts$Sp), 
                               rownames(FullRangeAdj1), 
                               #colnames(CytBMatrix),
                               colnames(FullSTMatrix),
                               rownames(HostAdj)), intersect)

FHN <- MammalHostNames; length(FHN)

HostThemselves <- # Removing diagonals, as they're uninformative
  which(upper.tri(HostAdj[FHN,FHN], diag = T)&lower.tri(HostAdj[FHN,FHN], diag  = T))

UpperHosts <- # Removing diagonals, as they're uninformative
  which(upper.tri(HostAdj[FHN,FHN], diag = T))

HostMatrixdf <- data.frame(Helminth = c(HostAdj[FHN, FHN]),
                           Space = c(FullRangeAdj1[FHN, FHN]),
                           #Phylo = c(tCytBMatrix[FHN, FHN]), 
                           Phylo = c(tFullSTMatrix[FHN, FHN]),
                           Sp = as.character(rep(FHN, each = length(FHN))),
                           Sp2 = as.character(rep(FHN, length(FHN)))
)

HostMatrixdf$Sp <- as.character(HostMatrixdf$Sp)
HostMatrixdf$Sp2 <- as.character(HostMatrixdf$Sp2)

UpperHosts <- # Removing diagonals and 
  which(upper.tri(HostAdj[FHN,FHN], diag = T))

FinalHostMatrix <- HostMatrixdf[-UpperHosts,]
FinalHostMatrix$HelminthBinary <- ifelse(FinalHostMatrix$Helminth>0, 1, 0)
FinalHostMatrix$Sp <- factor(FinalHostMatrix$Sp, levels = union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))
FinalHostMatrix$Sp2 <- factor(FinalHostMatrix$Sp2, levels = union(FinalHostMatrix$Sp,FinalHostMatrix$Sp2))

# Establishing all mammal data ####

# Replacing absent names in the full ST matrix ####

AllMammals <- intersect(colnames(FullSTMatrix),colnames(FullRangeAdj1))
AllMammals <- AllMammals[order(AllMammals)]
AbsentHosts <- FHN[which(!FHN%in%AllMammals)]

# LIST HERE

#rownames(FullSTMatrix) <- colnames(FullSTMatrix) <- sapply(rownames(FullSTMatrix), function(a) ifelse(a%in%AbsentHosts, NameReplace[a], a))

for(x in 1:length(WormGroups)){
  
  df <- reshape2::melt(HostAdjList[[x]]) %>%
    dplyr::rename(Sp = Var1, Sp2 = Var2)
  
  names(df)[3] <- WormGroups[x]  
  
  FinalHostMatrix <- FinalHostMatrix %>% left_join(df, by = c("Sp","Sp2"))
  
}

FinalHostMatrix[,paste0(WormGroups,"Binary")] <- apply(FinalHostMatrix[,WormGroups], 2, function(a) ifelse(a>0, 1, 0))


EltonTraits <- read.delim("Data/MamFuncDat.txt") %>% na.omit
EltonTraits$Scientific <- EltonTraits$Scientific %>% str_replace(" ","_")


EltonTraits$Carnivore <- ifelse(rowSums(EltonTraits[,c("Diet.Vend","Diet.Vect","Diet.Vfish")])>10,1,0)

DietComp <- EltonTraits %>% select(starts_with("Diet")) %>% select(1:10) %>% as.tibble
Remove <- which(rowSums(DietComp)==0)
DietComp <- DietComp %>% slice(-Remove)

VD <- vegdist(DietComp) %>% as.matrix

colnames(VD) <- rownames(VD) <- EltonTraits %>% slice(-Remove) %>% select(Scientific) %>% unlist

LongDiet <- reshape2::melt(VD) %>%     
  dplyr::rename(Sp = Var1, Sp2 = Var2, DietSim = value)

FinalHostMatrix <- FinalHostMatrix %>% left_join(LongDiet, by = c("Sp","Sp2")) %>%
  mutate(DietSim = 1 - DietSim)

FinalHostMatrix <- FinalHostMatrix %>% merge(EltonTraits[,c("Scientific","Carnivore")], by.x = "Sp", by.y = "Scientific")
FinalHostMatrix <- FinalHostMatrix %>% merge(EltonTraits[,c("Scientific","Carnivore")], by.x = "Sp2", by.y = "Scientific",
                                             suffixes = c("",".Sp2"))

FinalHostMatrix$Eaten <- ifelse(FinalHostMatrix$Carnivore==FinalHostMatrix$Carnivore.Sp2,0,1)

SlopeTime <- gather(FinalHostMatrix, key = "Group", value = "Shared", paste0(WormGroups,"Binary")) %>%
  filter(!is.na(Shared))

ggplot(SlopeTime, aes(Space, Shared, colour = Group)) + 
  facet_wrap(~Group) + 
  geom_point() + 
  geom_smooth()

ggplot(SlopeTime, aes(Phylo, Shared, colour = Group)) + 
  facet_wrap(~Group) + 
  geom_point() + 
  geom_smooth()

ggplot(SlopeTime, aes(DietSim, Shared, colour = Group, lty = as.factor(Eaten))) + 
  facet_wrap(~Group) + 
  geom_point() + 
  geom_smooth()
