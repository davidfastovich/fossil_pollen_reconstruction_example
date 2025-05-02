# Author: David Fastovich, dfastovi@syr.edu
# Date: 4/16/2025

###############
# LOAD PACKAGES
###############

library(openxlsx)
library(rioja)
library(analogue)
library(tidyverse)

#########################################
# SOURCE NEOTOMA FUNCTION TO COMPILE TAXA
#########################################

source("src/compile_taxa.R")

#####################
# READ IN EXCEL FILES
#####################

# New pollen counts
gds519 <- read.xlsx("pollen_counts/GDS519_for_tilia.xlsx")

# Data needed for analog calculations and climate reconstructions
modern_pollen <- read.xlsx("pollen_counts/whitmoreetal2005_v1-8.xlsx", sheet = "POLLEN DATA")
modern_climate <- read.xlsx("pollen_counts/whitmoreetal2005_v1-8.xlsx", sheet = "CLIMATE+BIOCLIMATE")
modern_pollen_splits <- read.xlsx("pollen_counts/whitmoreetal2005_v1-8.xlsx", sheet = "EAST_WEST") %>% 
  rename(ID2 = ID2_W05v1_6) # immediately rename the columns

####################################
# COMPILE TAXA NAMES TO COMMON NAMES
####################################

# Whitmore uses a non-standard naming convention so I'm going to change it to
# standard names usig a self-made conversion list.

from_whitmore <- read.csv("taxa_conversion_list/from_whitmore.csv")

modern_pollen_renamed <- as.data.frame(
  compile_taxa(modern_pollen[,14:147], alt.table = from_whitmore, list.name = "from_whitmore")
)

############################
# COMPILE TAXA NAMES TO WS64
############################

# Now that Whitmore has standard names we can use my extended WS64 list to
# compile taxa to the recommended list from Williams and Shuman (2008, QSR).

to_ws64 <- read.csv("taxa_conversion_list/ws64_extended.csv")

modern_pollen_ws64 <- as.data.frame(
  compile_taxa(modern_pollen_renamed, alt.table = to_ws64, list.name = "WS64")
) # Errors are okay because those taxa don't belong to WS64

modern_pollen_ws64$ID1 <- modern_pollen$ID1
modern_pollen_ws64$ID2 <- modern_pollen$ID2
modern_pollen_ws64$SITENAME <- modern_pollen$SITENAME

##########################################################################
# REMOVE ALL PICEA WEST SITES AND SPLIT PINUS INTO NORTHEAST AND SOUTHEAST
# POPULATIONS
##########################################################################

# Merge the relevant split data onto the pollen counts using the column
# ID2
modern_pollen_ws64_with_splits <- modern_pollen_splits %>% 
  select( # selecting only certain columns to keep
    ID2,
    PicEast,
    PinNE,
    PinSE
  ) %>% 
  left_join( # this is the merge command
    modern_pollen_ws64,
    by = "ID2"
  ) %>% 
  replace_na( # replacing NAs with 0s
    list(
      PicEast = 0,
      PinNE = 0,
      PinSE = 0
      )
  )

# Grab only the Picea east sites
modern_pollen_ws64_with_splits <- modern_pollen_ws64_with_splits %>% 
  filter(
    PicEast == 1
  )

############################################################
# SPLIT PINUS INTO NORTHEASTERN AND SOUTHEASTERN POPULATIONS
############################################################

# What I'm doing to split northeastern from southeastern pine populations:
# 1. Whitmore is not a perfect dataset and some sites are identified as both
# NE and SE pines. Those sites are removed.
# 2. Create new columns where the entries are 0 for NE and SE pines (undiff, 
# diplo, and haplo). This means six new columns filled with 0's.
# 3. Populate the new columns with the correct data. Taht is NE Pine undiff
# will only be filled with pine abundances from sites that are in the NE.
# 4. Remove the unsplit Pine columns.
# 5. Done. Only split Pine columns remain

# Step 1. Simple logic, add PinSE and PinNE together and remove any rows where
# the sum is two because that means that site is both NE and SE. 41 sites total
modern_pollen_ws64_with_splits <- modern_pollen_ws64_with_splits %>% 
  filter(
    modern_pollen_ws64_with_splits$PinNE + modern_pollen_ws64_with_splits$PinSE <= 1
    )

# Step 2.
modern_pollen_ws64_with_splits$`NE.Pinus diploxylon` <- 0
modern_pollen_ws64_with_splits$`NE.Pinus haploxylon` <- 0
modern_pollen_ws64_with_splits$`NE.Pinus undifferentiated` <- 0

modern_pollen_ws64_with_splits$`SE.Pinus diploxylon` <- 0
modern_pollen_ws64_with_splits$`SE.Pinus haploxylon` <- 0
modern_pollen_ws64_with_splits$`SE.Pinus undifferentiated` <- 0

# Step 3.
modern_pollen_ws64_with_splits$`NE.Pinus diploxylon`[modern_pollen_ws64_with_splits$PinNE == 1] <- modern_pollen_ws64_with_splits$`Pinus diploxylon`[modern_pollen_ws64_with_splits$PinNE == 1]
modern_pollen_ws64_with_splits$`NE.Pinus haploxylon`[modern_pollen_ws64_with_splits$PinNE == 1] <- modern_pollen_ws64_with_splits$`Pinus haploxylon`[modern_pollen_ws64_with_splits$PinNE == 1]
modern_pollen_ws64_with_splits$`NE.Pinus undifferentiated`[modern_pollen_ws64_with_splits$PinNE == 1] <- modern_pollen_ws64_with_splits$`Pinus undifferentiated`[modern_pollen_ws64_with_splits$PinNE == 1]

modern_pollen_ws64_with_splits$`SE.Pinus diploxylon`[modern_pollen_ws64_with_splits$PinSE == 1] <- modern_pollen_ws64_with_splits$`Pinus diploxylon`[modern_pollen_ws64_with_splits$PinSE == 1]
modern_pollen_ws64_with_splits$`SE.Pinus haploxylon`[modern_pollen_ws64_with_splits$PinSE == 1] <- modern_pollen_ws64_with_splits$`Pinus haploxylon`[modern_pollen_ws64_with_splits$PinSE == 1]
modern_pollen_ws64_with_splits$`SE.Pinus undifferentiated`[modern_pollen_ws64_with_splits$PinSE == 1] <- modern_pollen_ws64_with_splits$`Pinus undifferentiated`[modern_pollen_ws64_with_splits$PinSE == 1]

# Step 4.
modern_pollen_ws64_with_splits <- modern_pollen_ws64_with_splits %>% 
  select(
    !c(
      `Pinus diploxylon`,
      `Pinus haploxylon`,
      `Pinus undifferentiated`
    )
  )

#######################
# PREP THE CLIMATE DATA
#######################

# Select the corresponding climate data for our paired down training set
# Labeled "ready" because its the final data frame used for analyses
modern_climate_ready <- modern_climate %>% 
  filter(ID1 %in% modern_pollen_ws64_with_splits$ID1)

# Double check to make sure that the order is the same
all(modern_climate_ready$ID1 == modern_pollen_ws64_with_splits$ID1)

###########################################################
# REMOVE THE METADATA COLUMNS FROM THE MODERN POLLEN COUNTS
###########################################################

# Labeled "ready" because its the final data frame used for analyses
rownames(modern_pollen_ws64_with_splits) <- modern_pollen_ws64_with_splits$ID1

modern_pollen_ready <- modern_pollen_ws64_with_splits %>% 
  select(
    !c(
      ID1,
      ID2,
      SITENAME,
      PicEast,
      PinNE,
      PinSE
    )
  )

###########################
# COMPILE THE FOSSIL COUNTS
###########################

# Need to transpose data frame first to get species as columns
gds519_t <- as.data.frame(
  t(
    gds519[, -1]
  )
)
# Assign species names as columns while preserving the correct order
names(gds519_t) <- gds519[, 1]

# New pollen counts
gds519_ws64 <- as.data.frame(
  compile_taxa(gds519_t, alt.table = to_ws64, list.name = "WS64")
) # Polypodiophtya (zonate_ error okay - not in WS64

########################
# CONVERT TO PROPORTIONS
########################

gds519_ws64_pct <- gds519_ws64/rowSums(gds519_ws64)

modern_pollen_ready <- modern_pollen_ready/rowSums(modern_pollen_ready)

################################
# GET RID OF 'OTHER' AND POACEAE
################################

gds519_ws64_pct <- gds519_ws64_pct[,-grep("Other", colnames(gds519_ws64_pct))]

modern_pollen_ready <- modern_pollen_ready[,-grep("Other", colnames(modern_pollen_ready))]

##################################################################
# SPLIT FOSSIL COUNTS ACCORDING TO WORK FROM MARSICECK ET AL. 2018
##################################################################

# Southeastern from 0-336 cm all else northeastern
depth <- as.numeric(rownames(gds519_ws64_pct))

gds519_ws64_pct$`SE.Pinus haploxylon` <- gds519_ws64_pct$`Pinus haploxylon`
gds519_ws64_pct$`SE.Pinus undifferentiated` <- gds519_ws64_pct$`Pinus undifferentiated`
gds519_ws64_pct$`NE.Pinus haploxylon` <- gds519_ws64_pct$`Pinus haploxylon`
gds519_ws64_pct$`NE.Pinus undifferentiated` <- gds519_ws64_pct$`Pinus undifferentiated`

gds519_ws64_pct$`SE.Pinus haploxylon`[depth > 336] <- 0
gds519_ws64_pct$`SE.Pinus undifferentiated`[depth > 336] <- 0

gds519_ws64_pct$`NE.Pinus haploxylon`[depth <= 336] <- 0
gds519_ws64_pct$`NE.Pinus undifferentiated`[depth <= 336] <- 0

gds519_ws64_pct <- gds519_ws64_pct[,-c(25, 26)] # Get rid of unsplit columns

#################################
# 'JOIN' MODERN AND FOSSIL POLLEN
#################################

# Analog analysis requires that the species are identical in the training and 
# fossil dataset. This makes sure that the same species are present in both
# and any that are not observered as set to 0
gds519_ws64_pct_join <- analogue::join(gds519_ws64_pct, modern_pollen_ready)

# For instsance gds519_ws64_pct does not have cactaceae, but whitmore does.
# After this join GDS519 does hae cactacaea with 0's for all depths.

#################
# BUILD MAT MODEL
#################

# Build models
mat_model_tave <- rioja::MAT(
  y = modern_pollen_ready,
  x = modern_climate_ready$tave,
  dist.method ="sq.chord",
  k = 7,
  lean = FALSE
)

#######################################################
# PREDICT TEMPERATURE FROM THE FOSSIL POLLEN ABUNDANCES
#######################################################

# This also calculates the squared chord distance between the nearest modern
# pollen dataset and has all of the requested results
gds519_prediction <- predict(mat_model_tave, newdata = gds519_ws64_pct_join$gds519_ws64_pct, k = 7, sse = TRUE, n.boot = 100)

# Reconstructions and errors
gds519_prediction$fit.boot
gds519_prediction$SEP.boot

# Squared chord distance to samples
gds519_prediction$dist.n

# The row index for the nearest sample
gds519_prediction$match.name

# Create a Data Frame of the results and save as a csv
results <- data.frame(
  depth = as.numeric(rownames(gds519_t)),
  temperature = gds519_prediction$fit.boot,
  error = gds519_prediction$SEP.boot,
  analog1_ID1 = gds519_prediction$match.name[, 1],
  analog1_distance = gds519_prediction$dist.n[,1],
  analog1_name = sapply(gds519_prediction$match.name[, 1], function(x) modern_pollen$SITENAME[modern_pollen$ID1 == as.numeric(x)]),
  analog2_ID1 = gds519_prediction$match.name[, 2],
  analog2_distance = gds519_prediction$dist.n[,2],
  analog2_name = sapply(gds519_prediction$match.name[, 2], function(x) modern_pollen$SITENAME[modern_pollen$ID1 == as.numeric(x)]),
  analog3_ID1 = gds519_prediction$match.name[, 3],
  analog3_distance = gds519_prediction$dist.n[,3],
  analog3_name = sapply(gds519_prediction$match.name[, 3], function(x) modern_pollen$SITENAME[modern_pollen$ID1 == as.numeric(x)]),
  analog4_ID1 = gds519_prediction$match.name[, 4],
  analog4_distance = gds519_prediction$dist.n[,4],
  analog4_name = sapply(gds519_prediction$match.name[, 4], function(x) modern_pollen$SITENAME[modern_pollen$ID1 == as.numeric(x)]),
  analog5_ID1 = gds519_prediction$match.name[, 5],
  analog5_distance = gds519_prediction$dist.n[,5],
  analog5_name = sapply(gds519_prediction$match.name[, 5], function(x) modern_pollen$SITENAME[modern_pollen$ID1 == as.numeric(x)]),
  analog6_ID1 = gds519_prediction$match.name[, 6],
  analog6_distance = gds519_prediction$dist.n[,6],
  analog6_name = sapply(gds519_prediction$match.name[, 6], function(x) modern_pollen$SITENAME[modern_pollen$ID1 == as.numeric(x)]),
  analog7_ID1 = gds519_prediction$match.name[, 7],
  analog7_distance = gds519_prediction$dist.n[,7],
  analog7_name = sapply(gds519_prediction$match.name[, 7], function(x) modern_pollen$SITENAME[modern_pollen$ID1 == as.numeric(x)])
) %>% 
  write_csv("results.csv")

# Quickly plot the results
plot(
  x = results$depth,
  y = results$temperature.MAT,
  type = "l"
)

