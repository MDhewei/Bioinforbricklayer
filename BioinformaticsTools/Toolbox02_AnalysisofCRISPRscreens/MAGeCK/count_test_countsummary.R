Sweave("count_test_countsummary.Rnw");
library(tools);

texi2dvi("count_test_countsummary.tex",pdf=TRUE);

