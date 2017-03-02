
# Library -----------------------------------------------------------------

library(magrittr)
library(purrr)
library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)


# Import ------------------------------------------------------------------

dat <- NA #### IMPORT YOUR DATA HERE

# Data Checks -------------------------------------------------------------

# For our data, the columns are:
# metaboliteID: ID representing a unique measured metabolite.
# plateWellID: ID representing a sample in which metabolites were measured.
# plate, well: The plate and well that a sample was measured in. Every plate containts multiple wells
# value: The peak abundance measured in a sample. There are no NAs in this column.
# method: We generated a dataset with the above columns for the uncorrected data
  # and data after applying a correction algorithm. Those datasets were stacked
  # togehter, with this column identifying the dataset each row belongs to.

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
stopifnot(all(!is.na(dat$value)))
stopifnot(is.character(dat$method))


# Perc Misaligned CDF -----------------------------------------------------
# Only include metabolites that are 50% missing or less
maxPerc <- 0.5

misalign <-
  dat %>%
  # label any plates where a metaboliteID is totally missing
  group_by(method, metaboliteID, plate) %>%
  summarise(misaligned = all(value == 0)) %>%
  ungroup() %>%
  # Calculate % of plates totally missing
  group_by(method, metaboliteID) %>%
  summarise(misaligned = mean(misaligned)) %>%
  ungroup() %>%
  # Count the number of metaboliteIDs that are at any given misalignment %
  group_by(method, misaligned) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  # filter metabolites that are too missing
  filter(misaligned <= maxPerc) %>%
  # Calculate the cumulative count of metabolites that are 5% misaligned or better.
  group_by(method) %>%
  mutate(cumMisaligned = cumsum(n)) %>%
  ungroup()


# AUC ---------------------------------------------------------------------

auc <-
  misalign %>%
  group_by(method) %>%
  mutate(width = c(diff(misaligned), maxPerc - max(misaligned))) %>%
  mutate(auc = width * cumMisaligned) %>%
  summarise(auc = sum(auc),
            y = max(cumMisaligned)) %>%
  mutate(label = paste('AUCDF', round(auc, 1)),
         x = maxPerc) %>%
  select(method, label, x, y)


# No. Missing -------------------------------------------------------------

nMetabolites <-
  dat %>%
  group_by(method) %>%
  summarise(nMetabolite = length(unique(metaboliteID)))

perfaligned <-
  misalign %>%
  left_join(nMetabolites, by = 'method') %>%
  filter(misaligned == 0) %>%
  group_by(method) %>%
  summarise(label = cumMisaligned / nMetabolite,
            naligned = label * nMetabolite,
            x = -0.00,
            y = naligned) %>%
  mutate() %>%
  mutate(label = paste('No Misalignment Percentage:', round(label, 2))) %>%
  ungroup() %>%
  select(method, label, x, y)

cdftails <-
  misalign %>%
  group_by(method) %>%
  summarise(cumMisaligned = max(cumMisaligned)) %>%
  mutate(misaligned = maxPerc)

misalign %<>% full_join(cdftails, by = c("method", "misaligned", "cumMisaligned"))


# Plot --------------------------------------------------------------------


plotText <-
  rbind(perfaligned, auc)

p <-
  misalign %>%
  ggplot(aes(colour = method)) +
  geom_step(aes(x = misaligned, y = cumMisaligned), direction = 'hv', size = 1) +
  ylab('Cumulative Metabolite Count') +
  xlab('Percent Misalignment') +
  scale_colour_discrete('')

p +
  ggrepel::geom_label_repel(data = plotText, aes(label = label, x = x, y = y))
