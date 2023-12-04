# read in ABBA-BABA results
overall_result <- read.csv("result.Observed.txt", header = T, sep = "\t", na.strings = "")
window_result10k <- read.csv("window_10k_result.Observed.txt", header = T, sep = "\t", na.strings = "")
window_result5k <- read.csv("window_5k_result.Observed.txt", header = T, sep = "\t", na.strings = "")
window_result1k <- read.csv("window_1k_result.Observed.txt", header = T, sep = "\t", na.strings = "")
