# Banard exact tests 

myArgs <- as.numeric(commandArgs(trailingOnly = TRUE))
input <- matrix(myArgs, nrow=2, ncol=2, byrow = TRUE)

banard.output <- Barnard::barnard.test(input[1,1],input[2,1],input[1,2],input[2,2])
banard.output$p.value