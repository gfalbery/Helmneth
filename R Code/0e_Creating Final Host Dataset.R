# Creating final dataset

library(tidyverse)

if(file.exists("Data/FullRangeOverlap.Rdata")) load("Data/FullRangeOverlap.Rdata") else{
  source(paste0(CodeRoot,"/","0c2_Exhaustive Spatial Data Import.R"))
}

if(!file.exists("Data/intermediate/FullSTMatrix.csv")){
  
  library(geiger);library(ape);library(picante);library(dplyr)
  
  STFull <- read.nexus("data/ele_1307_sm_sa1.tre")[[1]]
  FullSTMatrix <- as.data.frame(cophenetic(STFull)) %>% as.matrix
  
} else{ FullSTMatrix <- as.matrix(read.csv("data/intermediate/FullSTMatrix.csv", header = T)) }

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

#AllMammals <- intersect(colnames(FullSTMatrix),colnames(FullRangeAdj1))
#AllMammals <- AllMammals[order(AllMammals)]
#AbsentHosts <- FHN[which(!FHN%in%AllMammals)]

# LIST HERE

#rownames(FullSTMatrix) <- colnames(FullSTMatrix) <- sapply(rownames(FullSTMatrix), function(a) ifelse(a%in%AbsentHosts, NameReplace[a], a))

for(x in 1:length(WormGroups)){
  
  df <- reshape2::melt(HostAdjList[[x]]) %>%
    rename(Sp = Var1, Sp2 = Var2)
  
  names(df)[3] <- WormGroups[x]  
  
  FinalHostMatrix <- FinalHostMatrix %>% left_join(df, by = c("Sp","Sp2"))
  
}

FinalHostMatrix[,paste0(WormGroups,"Binary")] <- apply(FinalHostMatrix[,WormGroups], 2, function(a) ifelse(a>0, 1, 0))

SlopeTime <- gather(FinalHostMatrix, key = "Group", value = "Shared", paste0(WormGroups,"Binary")) %>%
  filter(!is.na(Shared))

ggplot(SlopeTime, aes(Space, Shared, colour = Group)) + 
  facet_wrap(~Group) + 
  geom_point() + 
  geom_smooth()

