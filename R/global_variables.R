# This file declares global variables to suppress R CMD check warnings
# about undefined global variables in the package code.

utils::globalVariables(c(
  "Reads", "Group", "Sample", "CAP1", "CAP2",
  "subplot", "value", "Taxa", "Abundance",
  "Prevalence", "Genus", "Phylum", "qval", "ASV",
  ".", "Species", "cluster", "hue", "nChrHue", "shade", "prare2",
  "tax_fix", "phyloseq_validate", "ps_calc_dominant"
))
