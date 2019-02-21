# Data import for Helmneth ###

#remove(list = ls())

library(igraph); library(tidyverse); require(RCurl); library(readr); library(Matrix)

AssocsBase <- read_csv("Data/HelmAssocs.csv") %>% data.frame()
HostTraits <- read_csv("Data/CLC_database_hosts.csv") %>% data.frame()
HelminthTraits <- read_csv("Data/CLC_database_lifehistory.csv") %>% data.frame()

AssocsBase <- AssocsBase %>%
  rename(Helminth = Parasite) %>%
  mutate(Helminth = as.factor(Helminth), Host = as.factor(Host)) %>%
  mutate(Host = str_replace(Host, " ", "_"),
         Helminth = str_replace(Helminth, " ", "_"))

AssocsBase2 <- AssocsBase %>% filter(!Host == "Homo_sapiens")

# Making bipartite projections ####

AssocsTraits <- AssocsBase2[,c("Helminth","Host")]

m <- table(AssocsTraits)

attributes(m)$class <- "matrix"

M <- m %>%  as("dgCMatrix")

bipgraph <- graph.incidence(M, weighted = T)

Helminthgraph <- bipartite.projection(bipgraph)$proj1
Hostgraph <- bipartite.projection(bipgraph)$proj2

HelminthAdj <- as.matrix(get.adjacency(Helminthgraph, attr = "weight"))
diag(HelminthAdj) <- table(AssocsBase2$Helminth)

HostAdj <- as.matrix(get.adjacency(Hostgraph, attr = "weight"))
diag(HostAdj) <- table(AssocsBase2$Host)

# Deriving metrics from the networks ####

Hosts <- data.frame(Sp = names(V(Hostgraph)),
                    Degree = degree(Hostgraph),
                    Eigenvector = eigen_centrality(Hostgraph)$vector
                    #Kcore = coreness(Hostgraph),
                    #Between = betweenness(Hostgraph)
                    )

Helminths <- data.frame(Sp = names(V(Helminthgraph)),
                        Degree = degree(Helminthgraph),
                        Eigenvector = eigen_centrality(Helminthgraph)$vector
                        #Kcore = coreness(Helminthgraph),
                        #Between = betweenness(Helminthgraph)
                        )

Hosts <- merge(Hosts, HostTraits, by.x = "Sp", by.y = "hHostNameFinal", all.x = T)
Helminths <- merge(Helminths, HelminthTraits, by.x = "Sp", by.y = "vHelminthNameCorrected", all.x = T)

Hosts <- Hosts %>% dplyr::rename(hDom = hWildDomFAO)

Domestics <- Hosts[Hosts$hDom == "domestic", "Sp"]
Wildlife <- Hosts[Hosts$hDom == "wild", "Sp"]

DomesticHelminths <- as.factor(AssocsBase[AssocsBase$Host %in% Domestics, "Helminth"])
WildlifeHelminths <- as.factor(AssocsBase[AssocsBase$Host %in% Wildlife, "Helminth"])
HumanHelminths <- as.factor(AssocsBase[AssocsBase$Host == "Homo_sapiens", "Helminth"])

ZoonoticHelminths <- intersect(HumanHelminths, WildlifeHelminths)

AssocsTraits <- merge(AssocsTraits, HostTraits, by.x = "Host", by.y = "hHostNameFinal", all.x = T)
AssocsTraits <- merge(AssocsTraits, HelminthTraits, by.x = "Helminth", by.y = "vHelminthNameCorrected", all.x = T)

AssocsTraits$Domestic <- ifelse(AssocsTraits$Host%in%Domestics,1,0)
AssocsTraits$Wildlife <- ifelse(AssocsTraits$Host%in%Wildlife,1,0)
AssocsTraits$Zoonosis <- ifelse(AssocsTraits$Helminth%in%ZoonoticHelminths,1,0)

Hosts <- Hosts %>% 
  mutate(
    Domestic = ifelse(Sp %in% Domestics, 1, 0),
    Wildlife = ifelse(Sp %in% Wildlife, 1, 0),
    hZoonosisCount = c(table(AssocsTraits[AssocsTraits$Helminth%in%ZoonoticHelminths,"Host"])),
    Records = c(table(AssocsTraits$Host))
  ) %>%
  mutate(hZoonosisProp = hZoonosisCount/Records)

Helminths <- Helminths %>%
  mutate(
    Human = case_when(
      Sp %in% HumanHelminths ~ 1,
      TRUE ~ 0),
    
    Domestic = case_when(
      Sp %in% DomesticHelminths ~ 1,
      TRUE ~ 0),
    
    Wildlife = case_when(
      Sp %in% WildlifeHelminths ~ 1,
      TRUE ~ 0),
    
    DomesticCount = c(table(AssocsTraits[AssocsTraits$Domestic == 1,"Helminth"])),
    WildlifeCount = c(table(AssocsTraits[AssocsTraits$Wildlife == 1,"Helminth"])),
    ZoonosisCount = c(table(AssocsTraits[AssocsTraits$Zoonosis == 1,"Helminth"])),
    HumanCount = c(table(AssocsBase[AssocsBase$Host == "Homo_sapiens", "Helminth"]))[as.character(Helminths$Sp)],
    
    Records = c(table(AssocsTraits$Helminth))
  ) %>% mutate(
    
    
    PropDomestic = DomesticCount/Records,
    PropWildlife = WildlifeCount/Records,
    PropHuman = HumanCount/Records,
    
    PropZoonosis = ZoonosisCount/Records
    
  )


Helminths$HumDomWild <- factor(with(Helminths, 
                                    paste(ifelse(Human,"Human",""), 
                                          ifelse(Domestic,"Domestic",""), 
                                          ifelse(Wildlife,"Wild",""), sep = "")))

# Loading functions, determining themes ####

#devtools::install_github("gfalbery/ggregplot")
library(ggregplot); library(ggplot2); library(RColorBrewer)

ParasitePalettes<-c("PuRd","PuBu","BuGn","Purples","Oranges")
ParasiteColours<-c("#DD1c77","#2B8CBE","#2CA25F",brewer.pal(5,"Purples")[4],brewer.pal(5,"Oranges")[4])

AlberPalettes <- c("YlGnBu","Reds","BuPu", "PiYG")
AlberColours <- sapply(AlberPalettes, function(a) RColorBrewer::brewer.pal(5, a)[4])
AlberColours[length(AlberColours)+1:2] <- RColorBrewer::brewer.pal(11, AlberPalettes[[4]])[c(2,10)]

AlberTheme <- theme_bw() +
  theme(axis.title.x = element_text(vjust = -0.35, 
                                    size = 12, 
                                    colour = "black"), 
        axis.title.y = element_text(vjust = 1.2, 
                                    size = 12, 
                                    colour = "black"),
        strip.background = element_rect(fill = "white", colour = "dark grey"))

theme_set(AlberTheme)

# Trying out subgraphs ####

WormGroups <- unique(AssocsBase$group)

SubWorms <- lapply(WormGroups, function(a) AssocsBase %>% filter(group == a))

HelminthGraphs <- HostGraphs <- HelminthAdjList <- HostAdjList <-  list()

for(w in 1:length(SubWorms)){
  
  m <- SubWorms[[w]] %>% select(Helminth, Host) %>% table()
  
  attributes(m)$class <- "matrix"
  
  M <- m %>%  as("dgCMatrix")
  
  bipgraph <- graph.incidence(M, weighted = T)
  
  HelminthGraphs[[w]] <- bipartite.projection(bipgraph)$proj1
  HostGraphs[[w]] <- bipartite.projection(bipgraph)$proj2
  
  HelminthAdjList[[w]] <- as.matrix(get.adjacency(HelminthGraphs[[w]], attr = "weight"))
  diag(HelminthAdjList[[w]]) <- SubWorms[[w]] %>% select(Helminth) %>% table()
  
  HostAdjList[[w]] <- as.matrix(get.adjacency(HostGraphs[[w]], attr = "weight"))
  diag(HostAdjList[[w]]) <- SubWorms[[w]] %>% select(Host) %>% table()
  
}

names(SubWorms) <- WormGroups
