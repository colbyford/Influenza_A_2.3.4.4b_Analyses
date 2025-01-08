## Load in StrainHub
library(strainhub)

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

print(graph)


## Make Leaflet/Swoopy map
# make_map(graph,
#          geodata = geodata,
#          columnSelection = "Country")
