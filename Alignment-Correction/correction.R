
# Library -----------------------------------------------------------------

### easy-read syntax and data-manipulation ###
library(magrittr)
library(dplyr)
library(purrr)

### Data import ###
library(readr)

### splines ###
library(gam)

### plotting ###
library(ggplot2)

### xml file manipulation ###
library(rvest)
library(xml2)



# Import ------------------------------------------------------------------

dat <- NA #### IMPORT YOUR DATA HERE

# Data Checks -------------------------------------------------------------
# For our data, the columns are:
# metaboliteID: ID representing a unique measured metabolite.
# plateWellID: ID representing a sample in which metabolites were measured.
# retentionTimeMean: The mean retention time of metabolite across all samples.
# The mean is the referent value that we attempt to adjust values to.
# retentionTime: The measured retention time of a metabolite in a given sample.
# diff: The difference between the referent value and the measured value. These
# differences are what we attempt to model as a function of retention time.


# THIS SCRIPT IS ONLY GUARANTEED TO WORK IF YOUR DATA HAS THESE EXACT COLUMNS
# WITH THE EXACT SAME NAMES WITH THE EXACT SAME INFORMATION IN THE EXACT DATA
# FORMAT.

# Below are some checks that attempt to help ensure that the data is formatted
# correctly before proceeding. However, passing the checks does *not* guarentee
# that the script will work for your data

# Check for column Names
stopifnot(all(c('metaboliteID', 'plateWellID',
                'retentionTimeMean', 'retentionTime', 'diff') %in% names(dat)))

# Check for data format
stopifnot(is.numeric(dat$metaboliteID))
stopifnot(is.character(dat$plateWellID))
stopifnot(is.numeric(dat$retentionTimeMean))
stopifnot(is.numeric(dat$retentionTime))
stopifnot(is.numeric(dat$diff))


# Plots of Difference from Referent ---------------------------------------

dat %>%
  ggplot(aes(x = retentionTimeMean, y = diff, colour = plateWellID)) +
  geom_point() +
  scale_color_discrete(guide = FALSE) +
  stat_smooth(colour = 'black') +
  facet_wrap(~ plateWellID) +
  ggtitle('Difference from Referent for 6 Samples') +
  ylab('Difference from Referent (minutes)') +
  xlab('Referent Retention Time (minutes)')


dat %>%
  ggplot(aes(x = retentionTimeMean, y = diff, colour = plateWellID)) +
  geom_point() +
  stat_smooth(colour = 'black') +
  ggtitle('Difference from Referent for 6 Samples') +
  ylab('Difference from Referent (minutes)') +
  xlab('Referent Retention Time (minutes)')


# Add time based domain knowledge -----------------------------------------

time_region <- function(x) {
  ifelse(x <= 1, 1,
         ifelse(x >= 6, 3, 2))
}

dat %<>%
  mutate(timeRegion = time_region(retentionTime)) %>%
  mutate(timeRegion = as.factor(timeRegion))


# Modeling ----------------------------------------------------------------

models <-
  dat %>%
  split(.$plateWellID) %>%
  map(function(x) {
    gam(
      diff ~ s(retentionTime, 16) + timeRegion + timeRegion * s(retentionTime, 16) ,
      data = x
    )
  })



# mzXML Processing --------------------------------------------------------

# Used with mzMine v3.2
# YOU CANNONT USE THIS AS IT IS. YOU MUST MAKE ADJUSTMENTS. SEARCH for '#**'

dir <-
  'data/folderWithMZMINExmls/' #**

outputpath <-
  'path/to/output/folder' #**

suffix <-
  'stringToTagOutputfileWith' #**

mzXMLdocs <- list.files(dir)

for (d in mzXMLdocs) {

  # Get the platewell identifier from the fileName
  pwid <-
    stringr::str_extract(d, 'plateWellID-regex') #**

  doc <-
    read_xml(paste0(dir, d))

  # Extract Scan Elements and retention times from XML
  scanNodes <-
    xml_children(doc)[[1]] %>%
    xml_children() %>%
    .[xml_name(.) == "scan"]

  # cleaning and extraction
  retentionTimes <-
    xml_attr(scanNodes, 'retentionTime') %>%
    {
      gsub('^PT', '', .)
    } %>%
    {
      gsub('S$', '', .)
    } %>%
    as.numeric() %>%
    {
      . / 60 #  minute/second conversion
    }

  #creating a dataframe of values to predict on
  predDat <-
    data.frame(retentionTimes,
               time_region(retentionTimes)) %>%
    setNames(c('retentionTime', 'timeRegion'))

  predDat$timeRegion %<>% factor(levels = levels(dat$timeRegion))

  # Grab the model corresponding to the plate and well combo of the data
  model <-
    models[[pwid]]

  predDiffs <-
    predict(model, newdata = predDat) %>%
    as.vector()

  correctedRetentionTimes <-
    retentionTimes %>%
    `-`(predDiffs) %>%
    `*`(60) # Minute/second conversion

  # Adjust for those that were over corrected
  # Corrected final values must be constantly increasing (software constraint)
  overcorrected <-
    c(FALSE, diff(correctedRetentionTimes) <= 0)

  while (any(overcorrected)) {
    # index of first overcorrected time
    oc <-
      which.max(overcorrected)

    # Overcorrected time, and all times after it, are shifted such that the
    # overcorrected time is delta greater than the time previous
    delta <- 0.0001

    correction <-
      (correctedRetentionTimes[(oc - 1)] + delta) - correctedRetentionTimes[(oc)]

    correctedRetentionTimes[oc:length(correctedRetentionTimes)] <-
      correctedRetentionTimes[oc:length(correctedRetentionTimes)] + correction

    overcorrected <-
      c(FALSE, diff(correctedRetentionTimes) <= 0)
  }

  stopifnot(all(correctedRetentionTimes > 0))
  stopifnot(all(diff(correctedRetentionTimes) > 0))

  correctedFormattedRetentionTimes <-
    correctedRetentionTimes %>%
    {
      paste0("PT", ., "S")
    }

  for (i in 1:length(correctedFormattedRetentionTimes)) {
    xml_attr(scanNodes[i], 'retentionTime') <-
      correctedFormattedRetentionTimes[i]
  }

  stopifnot(all(xml_attr(scanNodes, 'retentionTime') == correctedFormattedRetentionTimes))

  write_xml(doc, paste0(outputpath, suffix, '-', d))
}
