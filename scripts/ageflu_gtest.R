# G-test of independence

# Get 2x2 contingency table inputs 
myArgs <- as.numeric(commandArgs(trailingOnly = TRUE))
input <- matrix(myArgs, nrow=2, ncol=2, byrow = TRUE)

output <- DescTools::GTest(input, correct="yates")
output$p.value