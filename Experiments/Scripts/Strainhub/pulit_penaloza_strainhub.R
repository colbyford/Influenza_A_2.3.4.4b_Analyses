## Load in StrainHub
library(strainhub)

## Read in tree, metadata, and geodata
treedata <- ape::read.nexus("../data/H5N1_N2_HA_nuc_djv28_Pulit-Penaloza-no-gap-0-is-2024-host-only.nexus")
metadata <- readr::read_csv("../data/H5N1_N2_HA_nuc_djv28_Pulit-Penaloza-no-gap-0-is-2024-host-only_metadata.csv", col_names = TRUE)
# geodata <- readr::read_csv("", col_names = TRUE)


## Check to See Which States are available by which to generate the network
list_states(treedata,
            metadata,
            treeType = "parsimonious")


## Make the Transmission Network
graph <- make_transnet(treedata,
                       metadata,
                       columnSelection = "Source",
                       centralityMetric = 6,
                       treeType = "parsimonious")

print(graph)


## Make Leaflet/Swoopy map
# make_map(graph,
#          geodata = geodata,
#          columnSelection = "Country")
