## Load in StrainHub
library(strainhub)
library(dplyr)

## Read in tree, metadata, and geodata
treedata <- ape::read.tree("../data/clade-of-interest-around-76-human/RAxML_bestTree.H5N1_N2_HA_nuc_djv28_Pulit-Penaloza-blast-LA-NM_noDup.clade-of-interest-around-76-human.out")
metadata <- readr::read_csv("../data/clade-of-interest-around-76-human/H5N1_N2_HA_nuc_djv28_Pulit-Penaloza-blast-LA-NM_noDup.host.rollup1must-fix-Khaled.clade-of-interest-76-human-demerge.csv", col_names = TRUE)

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
           label == "Carnivora" ~ "darkred",
           label == "Artiodactyla" ~ "darkred",
           label == "Galliformes" ~ "darkred",
           TRUE ~ "black"
         ),
         font.size = 40)


edges <- graph$x$edges %>%
  mutate(
    arrows = "to",
    smooth = TRUE,
    # color = ifelse(to == 12, "red", "grey"),
    color = case_when(
      to == 10 ~ "red",
      from == 10 ~ "darkred",
      TRUE ~ "lightgrey"
    ),
    # width = ifelse(value == 1, 1, 4),
    # value = NULL
    )

visNetwork(nodes, edges) %>% 
  # visPhysics(solver = "forceAtlas2Based", 
  #            forceAtlas2Based = list(gravitationalConstant = -60)) %>% 
  visPhysics(solver = "forceAtlas2Based",
             forceAtlas2Based = list(gravitationalConstant = -100))
  
