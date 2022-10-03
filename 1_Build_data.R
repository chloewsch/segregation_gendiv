# 1. Build data ------------
library(tidyverse)
library(tigris)
library(sf)
library(rgeos)

# Merge census data with spatial data ------
## Load IPUMS demographic data:
p5 <- read.csv("censusdata_2010_block.csv", header = TRUE)

# Filter, rename column codes
# CBSAA: Metropolitan Statistical Area/Micropolitan Statistical Area Code
# UAA: Urban Area Code

p5_fil <- p5 %>% 
  filter(!STATE %in% c("Hawaii", "Puerto Rico"),
         H7Z001 > 0) %>% # this is total # of people in block; remove blocks with no people.
  select(GISJOIN, CBSAA, UAA, H7Z001, H7Z002,	H7Z003, H7Z004, H7Z005, H7Z006, H7Z007, 
         H7Z008, H7Z009, H7Z010, H7Z011, H7Z012, H7Z013, H7Z014, H7Z015, H7Z016, H7Z017) %>% 
  rename(total = H7Z001, non_hispanic = H7Z002,	
         white = H7Z003, black =	H7Z004, indigenous =	H7Z005, asian = H7Z006,
         pacific_islander = H7Z007, other = H7Z008, mixed =	H7Z009,
         hispanic = H7Z010, white_hisp	= H7Z011, black_hisp =	H7Z012, indi_hisp = H7Z013,
         asian_hisp = H7Z014, pi_hisp =	H7Z015, other_hisp =	H7Z016, mixed_hisp =	H7Z017)

rm(p5)

# Attached demographic data to a map:
## Create a join column with unique IDs for blocks
# In block shapefile (tigris/census): State = 2 digits, County = 3 digits, Tract = 6 digits, Block = 4 digits

# NHGIS data has GISJOIN code, e.g.: G01000100201001000 = G (01) 0 (001) 0 (020100) (1000)
#                                                           STATE  COUNTY   TRACT    BLOCK

# Create the same column in census shapefile 

## USA census blocks:
statenames <- as.list(c("AL","AK","AZ","AR", "CA", "CO", "CT", "DE", "FL", "GA", "ID", "IL", "IN", "IA", "KS", "KY", 
                        "LA", "ME","MD", "MA", "MI", "MN", "MS", "MO", "MT", "NE", "NV", "NH", "NJ", "NM", "NY", "NC", "ND", 
                        "OH", "OK", "OR", "PA", "RI", "SC", "SD", "TN", "TX", "UT", "VT", "VA", "WA", "WV", "WI", "WY", "DC"))

blocks_list <- list()

for (i in statenames) {
  tryCatch(
    expr = {
      b <- blocks(i)
    },
    error = function(e){cat("ERROR :",conditionMessage(e), "\n")})
  blocks_list[[i]] <- b
} 
# Note: large objects, sometimes fails. The last successful state is duplicated
# Repeat for states that fail as needed

usblocks <- do.call(rbind, blocks_list)

# Create join column:
usblocks$GISJOIN <- paste("G", usblocks$STATEFP10, 0, usblocks$COUNTYFP10, 0, usblocks$TRACTCE10, usblocks$BLOCKCE10, sep = "")

# Remove duplicates if any
usblox <- usblocks[,c(18,19)] # only take the Geometry column to decrease size
un <- unique(usblox$GISJOIN)
matches <- match(un, usblox$GISJOIN)
usblocks_uni <- usblox[c(matches),]
#rm(usblocks, usblox, matches, un)

## Merge spatial block data with demographic data -----------
usblocks_merge <-  merge(usblocks, p5_fil, by = "GISJOIN", all.x = FALSE, all.y = TRUE)

rm(p5_fil, usblocks)

# Summarizing racial composition for genetic data sites --------------
# Load genetic data:
gdata <- read.csv("Schmidt_Garroway_PNAS_2022_analysisdata.csv", header = TRUE)

# Convert to sf
sites <- sf::st_as_sf(gdata, coords = c("lon", "lat"), crs = 4326) # WGS84

# Convert to equal area for buffering in m:
AEA <- "+proj=aea +lat_1=20 +lat_2=60 +lat_0=40 +lon_0=-96 +x_0=0 +y_0=0 +ellps=GRS80 
+datum=NAD83 +units=m +no_defs"

sitesAEA <- st_transform(sites, crs = AEA)

# Create buffer around sites (0km, 0.5km, 1km, 5km):
sites_buffAEA.5 <- st_buffer(sitesAEA, dist = 500)
sites_buffAEA1 <- st_buffer(sitesAEA, dist = 1000)
sites_buffAEA5 <- st_buffer(sitesAEA, dist = 5000)

# Convert to NAD83 (same as blocks; EPSG code 4269):
sites_buff0 <- st_transform(sites, crs = st_crs(usblocks_merge)) # no buffer (point feature)
sites_buff.5 <- st_transform(sites_buffAEA.5, crs = st_crs(usblocks_merge))
sites_buff1 <- st_transform(sites_buffAEA1, crs = st_crs(usblocks_merge))
sites_buff5 <- st_transform(sites_buffAEA5, crs = st_crs(usblocks_merge))

# Extract demographic data per site
# This returns a dataframe with 1 row for every block a site overlaps
cendat0 <- st_join(sites_buff0, usblocks_merge)
cendat.5 <- st_join(sites_buff.5, usblocks_merge)
cendat1 <- st_join(sites_buff1, usblocks_merge)
cendat5 <- st_join(sites_buff5, usblocks_merge)

# Average the # of different races from each block:
#### 0.5 km BUFFER
dat.5 <- cendat.5 %>% 
  group_by(pop) %>% 
  summarise(total.5 = mean(total, na.rm = TRUE), # returns NaN when there are no overlaps
            # Non-hispanic or latino:
            non_hispanic.5 = mean (non_hispanic, na.rm = TRUE),
            white.5 = mean (white, na.rm = TRUE),  
            black.5 = mean (black, na.rm = TRUE),
            indigenous.5 = mean (indigenous, na.rm = TRUE),
            asian.5 = mean (asian, na.rm = TRUE),
            pacific_isl.5 = mean (pacific_islander, na.rm = TRUE),
            other.5 = mean (other, na.rm = TRUE),
            mixed.5 = mean (mixed, na.rm = TRUE),
            # Hispanic or latino:
            hispanic.5 = mean (hispanic, na.rm = TRUE),
            white_hisp.5 = mean (white_hisp, na.rm = TRUE),
            black_hisp.5 = mean (black_hisp, na.rm = TRUE),
            indig_hisp.5 = mean (indi_hisp, na.rm = TRUE),
            asian_hisp.5 = mean (asian_hisp, na.rm = TRUE),
            pi_hisp.5 = mean (pi_hisp, na.rm = TRUE),
            other_hisp.5 = mean (other_hisp, na.rm = TRUE),
            mixed_hisp.5 = mean (mixed_hisp, na.rm = TRUE)) %>%
  drop_na() %>%  # remove populations with no block overlaps
  as.data.frame() %>% 
  dplyr::select(-geometry)

#### 1 km BUFFER
dat1 <- cendat1 %>% 
  group_by(pop) %>% 
  summarise(total1 = mean(total, na.rm = TRUE), # returns NaN when there are no overlaps
            # Non-hispanic or latino:
            non_hispanic1 = mean (non_hispanic, na.rm = TRUE),
            white1 = mean (white, na.rm = TRUE),  
            black1 = mean (black, na.rm = TRUE),
            indigenous1 = mean (indigenous, na.rm = TRUE),
            asian1 = mean (asian, na.rm = TRUE),
            pacific_isl1 = mean (pacific_islander, na.rm = TRUE),
            other1 = mean (other, na.rm = TRUE),
            mixed1 = mean (mixed, na.rm = TRUE),
            # Hispanic or latino:
            hispanic1 = mean (hispanic, na.rm = TRUE),
            white_hisp1 = mean (white_hisp, na.rm = TRUE),
            black_hisp1 = mean (black_hisp, na.rm = TRUE),
            indig_hisp1 = mean (indi_hisp, na.rm = TRUE),
            asian_hisp1 = mean (asian_hisp, na.rm = TRUE),
            pi_hisp1 = mean (pi_hisp, na.rm = TRUE),
            other_hisp1 = mean (other_hisp, na.rm = TRUE),
            mixed_hisp1 = mean (mixed_hisp, na.rm = TRUE)) %>%
  drop_na() %>%  # remove populations with no block overlaps
  as.data.frame() %>% 
  dplyr::select(-geometry)

#### 5 km BUFFER
dat5 <- cendat5 %>% 
  group_by(pop) %>% 
  summarise(total5 = mean(total, na.rm = TRUE), # returns NaN when there are no overlaps
            # Non-hispanic or latino:
            non_hispanic5 = mean (non_hispanic, na.rm = TRUE),
            white5 = mean (white, na.rm = TRUE),  
            black5 = mean (black, na.rm = TRUE),
            indigenous5 = mean (indigenous, na.rm = TRUE),
            asian5 = mean (asian, na.rm = TRUE),
            pacific_isl5 = mean (pacific_islander, na.rm = TRUE),
            other5 = mean (other, na.rm = TRUE),
            mixed5 = mean (mixed, na.rm = TRUE),
            # Hispanic or latino:
            hispanic5 = mean (hispanic, na.rm = TRUE),
            white_hisp5 = mean (white_hisp, na.rm = TRUE),
            black_hisp5 = mean (black_hisp, na.rm = TRUE),
            indig_hisp5 = mean (indi_hisp, na.rm = TRUE),
            asian_hisp5 = mean (asian_hisp, na.rm = TRUE),
            pi_hisp5 = mean (pi_hisp, na.rm = TRUE),
            other_hisp5 = mean (other_hisp, na.rm = TRUE),
            mixed_hisp5 = mean (mixed_hisp, na.rm = TRUE)) %>%
  drop_na() %>%  # remove populations with no block overlaps
  as.data.frame() %>% 
  dplyr::select(-geometry)

# Attach the rest of the genetic data:
gdata_cen <- Reduce(function(x,y) merge(x,y, by = "pop", all=TRUE), 
                    list(gdata, dat0, dat.5, dat1, dat5))

# Data for final analyses:
gdatacen <- gdata_cen %>% 
  select(c(1:13),
         white.5, total.5,
         white1, total1,
         white5, total5)

write.csv(gdatacen, "SchmidtGarroway_PNAS_2022_analysisdata.csv", quote = FALSE, row.names = FALSE)