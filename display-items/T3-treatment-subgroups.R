library(dplyr)

setwd('UCSF/ipd-ma-cd2-2')

source('helper/subgroup.R') # subgroup.cohorts()

# load preprocessed data
drug_ordering <- read.csv('data/treatment_preferences-n10000.csv')

# group patients by similar treatment preferences
subgroup.cohorts(drug_ordering)
#   drug1 drug2 drug3 p12_ohe p23_ohe     n
# 1 drug1 drug2 drug3       0       0  3142
# 2 drug1 drug2 il12        0       1     4
# 3 drug1 drug2 intg        0       1   354
# 4 il12  drug2 drug3       1       0   138
# 5 il12  tnfi  intg        1       1     1
# 6 tnfi  drug2 drug3       1       0  2021
# 7 tnfi  intg  il12        1       1    43