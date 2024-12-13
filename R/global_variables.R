# This file declares global variables to suppress R CMD check warnings
# about undefined global variables in the package code.

utils::globalVariables(c(
  "Reads", "Group", "Sample", "CAP1", "CAP2",
  "subplot", "value", "Taxa", "Abundance",
  "Prevalence", "Genus", "Phylum", "qval"
))
