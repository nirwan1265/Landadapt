library(dplyr)

# Load maize and sorghum landraces
maize_sorghum <- read.csv("data/geoloc/taxa_geoloc_pheno.csv") %>% dplyr::select(c(2,3,6,7))

maize <- maize_sorghum %>% filter(sp == "Zea mays") %>% dplyr::select(-sp)
names(maize) <- c("Lines", "Long","Lat")

sorghum <- maize_sorghum %>% filter(sp == "Sorghum bicolor") %>% dplyr::select(-sp)
names(sorghum) <- c("Lines", "Long","Lat")

# Load arabidopsis
at <- read.csv("data/geoloc/at_1001_genome_geoloc.csv") %>% dplyr::select(c(1,4,5))
names(at) <- c("Lines", "Lat","Long")

# Load Rice
rice <- read.csv("data/geoloc/3000_RG_geoloc.csv") %>% dplyr::select(c(1,2,3))
names(rice) <- c("Lines", "Long","Lat")
