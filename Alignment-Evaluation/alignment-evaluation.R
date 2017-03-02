
# Library -----------------------------------------------------------------

# Data manipulation
library(magrittr)
library(purrr)
library(dplyr)
library(tidyr)

# Plotting
library(ggplot2)

# Shape-based clustering
library(dtwclust)


# Import ------------------------------------------------------------------

dat <- NA #### IMPORT YOUR DATA HERE


# Data Checks -------------------------------------------------------------

# For our data, the columns are:
# metaboliteID: ID representing a unique measured metabolite.
# plateWellID: ID representing a sample in which metabolites were measured.
# plate, well: The plate and well that a sample was measured in. Every plate containts multiple wells
# value: The peak abundance measured in a sample. There are no NAs in this column.

# THIS SCRIPT IS ONLY GUARANTEED TO WORK IF YOUR DATA HAS THESE EXACT COLUMNS
# WITH THE EXACT SAME NAMES WITH THE EXACT SAME INFORMATION IN THE EXACT DATA
# FORMAT.

# Below are some checks that attempt to help ensure that the data is formatted
# correctly before proceeding. However, passing the checks does *not* guarentee
# that the script will work for your data

# Check for column Names
stopifnot(all(c('metaboliteID', 'plateWellID',
                'plate', 'well', 'value', 'method') %in% names(dat)))

# Check for data format
stopifnot(is.numeric(dat$metaboliteID))
stopifnot(is.character(dat$plateWellID))
stopifnot(is.integer(dat$plate))
stopifnot(is.integer(dat$well))
stopifnot(is.numeric(dat$value))


# Misalignment ------------------------------------------------------------


misalign <-
  dat %>%
  group_by(metaboliteID, plate) %>%
  summarise(misaligned = all(is.na(value))) %>%
  ungroup() %>%
  group_by(metaboliteID) %>%
  mutate(misaligned = mean(misaligned)) %>%
  ungroup()


dat %<>% left_join(misalign, by = c('metaboliteID', 'plate'))


# Standard Plotting -------------------------------------------------------


metabsToPlot <-
  dat %>%
  distinct(metaboliteID) %>%
  sample_n(10) # random 100 metabolites

blankThm <-
  theme(
    text = element_blank(),
    strip.background = element_blank(),
    axis.ticks = element_blank()
  )

# All Metabolites PLotting
dat %>%
  ggplot() +
  geom_point(aes(x = plateWellID, y = value, colour = misaligned)) +
  facet_wrap( ~ metaboliteID, scales = 'free_y') +
  scale_color_gradient(low = 'blue', high = 'red') +
  blankThm

# Subset PLotting
dat %>%
  inner_join(metabsToPlot, by = 'metaboliteID') %>%
  ggplot() +
  geom_point(aes(x = plateWellID, y = value, colour = misaligned)) +
  facet_wrap( ~ metaboliteID, scales = 'free_y') +
  scale_color_gradient(low = 'blue', high = 'red') +
  blankThm


# Shape Based Clustering --------------------------------------------------

# Format the data for the clustering function
tsdat <-
  dat %>%
  select(metaboliteID, value) %>%
  split(.$metaboliteID) %>%
  map(function(x) {
    x[, 'value']
  })

nms <- names(tsdat)

tsdat %<>%
  bind_cols %>%
  set_names(nms) %>%
  as.data.frame %>%
  t()

# Here, we replace NA's with zero. Recognizing this step can be important for understanding why certain metabolites may end up clustered together, especially when trying to understand why a mostly misalgned metabolite might be assigned to a specific cluster.
tsdat[is.na(tsdat)] <- 0


# Compute Hierarchical clusters
clusts <-
  dtwclust(tsdat,
           k = 10, # Number of clusters
           type = 'hierarchical', # Clustering technique
           method = 'ward.D2', # linkage method
           distance = 'sbd', # Shape based distance
           preproc = zscore)

# Get the cluster labels
clusterLabels <-
  data.frame(metaboliteID = nms,
             label =  cutree(clusts, k = 10))


# Cluster Plot Backgrounds ------------------------------------------------
# For coloured backgrounds
rectangles <-
  dat %>%
  group_by(metaboliteID) %>%
  arrange(plateWellID) %>%
  summarise(xmin = plateWellID[1],
            xmax = plateWellID[n()],
            ymin = min(value, na.rm = T),
            ymax = max(value, na.rm = T)) %>%
  left_join(clusterLabels) %>%
  mutate(label = as.factor(label))


# Cluster Plots -----------------------------------------------------------
#Full Plot
dat %>%
  left_join(clusterLabels, by  = 'metaboliteID') %>%
  ggplot() +
  geom_point(aes(x = plateWellID, y = value, colour = misaligned)) +
  facet_wrap( ~ label + metaboliteID, scales = 'free_y') +
  scale_color_gradient(low = 'blue', high = 'red') +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = label),
            alpha = 0.3,
            data = rectangles) +
  blankThm


# Subset Plot
dat %>%
  inner_join(metabsToPlot, by = 'metaboliteID') %>%
  left_join(clusterLabels, by  = 'metaboliteID') %>%
  # mutate(plateWellID = as.numeric(as.factor(plateWellID))) %>%
  ggplot() +
  geom_point(aes(x = plateWellID, y = value, colour = misaligned)) +
  facet_wrap( ~ label + metaboliteID, scales = 'free_y') +
  scale_color_gradient(low = 'blue', high = 'red') +
  geom_rect(aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = label),
            alpha = 0.3,
            data = inner_join(rectangles, metabsToPlot, by = 'metaboliteID')) +
  blankThm

