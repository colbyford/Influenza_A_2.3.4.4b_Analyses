## Load in StrainHub
library(strainhub)
library(dplyr)

## Read in tree, metadata, and geodata
# treedata <- ape::read.nexus("../data/clade-of-interest-around-PQ809550-strainhub/tree0.nexus")
# metadata <- readr::read_csv("../data/clade-of-interest-around-PQ809550-strainhub/H5N1_N2_HA_nuc_djv28_Pulit-Penaloza-blast-LA-NM_noDup.host.rollup1.clade-of-interest.csv", col_names = TRUE)
treedata <- ape::read.tree("../data/clade-of-interest-around-PQ809550-strainhub/Jan8/RAxML_bestTree.H5N1_N2_HA_nuc_djv28_Pulit-Penaloza-clade-of-interest-around-PQ809550.out")
metadata <- readr::read_csv("../data/clade-of-interest-around-PQ809550-strainhub/Jan8/H5N1_N2_HA_nuc_djv28_Pulit-Penaloza-blast-LA-NM_noDup.host.rollup1.clade-of-interest-good-jan8.csv", col_names = TRUE)

# geodata <- readr::read_csv("", col_names = TRUE)


## Check to See Which States are available by which to generate the network
list_states(treedata,
            metadata,
            treeType = "parsimonious")


## Make the Transmission Network
graph <- make_transnet(treedata,
                       metadata,
                       columnSelection = "Source_Clean_FIX",
                       centralityMetric = 6,
                       treeType = "parsimonious")

print(
  graph %>% 
    visNetwork::visPhysics(solver = "forceAtlas2Based", 
               forceAtlas2Based = list(gravitationalConstant = -30))
    )


## Network Customization
library(visNetwork)

nodes <- graph$x$nodes %>%
  mutate(shape = "dot",
         # color = ifelse(label == "Primates", "red", "black"),
         color = case_when(
           label == "Primates" ~ "red",
           label == "Anseriformes" ~ "darkred",
           label == "Galliformes" ~ "darkred",
           TRUE ~ "black"
         ),
         font.size = 20)
  


edges <- graph$x$edges %>%
  mutate(
    arrows = "to",
    smooth = TRUE,
    # color = ifelse(to == 12, "red", "grey"),
    color = case_when(
      to == 12 ~ "red",
      from == 12 ~ "darkred",
      TRUE ~ "lightgrey"
    ),
    # width = ifelse(value == 1, 1, 4),
    # value = NULL
    )

visNetwork(nodes, edges) %>% 
  visPhysics(solver = "forceAtlas2Based", 
             forceAtlas2Based = list(gravitationalConstant = -60))
  
