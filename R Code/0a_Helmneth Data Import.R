# Data import for Helmneth ###

#remove(list = ls())

library(igraph); library(tidyverse); require(RCurl); library(readr); library(Matrix)

AssocsBase <- read_csv("Data/HelmAssocs.csv") %>% data.frame()
HostTraits <- read_csv("Data/CLC_database_hosts.csv") %>% data.frame()
HelminthTraits <- read_csv("Data/CLC_database_lifehistory.csv") %>% data.frame()

names(AssocsBase)[3] <- 'Helminth'
AssocsBase <- AssocsBase %>%
  #dplyr::rename(Helminth = Parasite) %>%
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
