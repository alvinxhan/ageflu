library(PropCIs)
myArgs <- as.numeric(commandArgs(trailingOnly = TRUE))
PropCIs::orscoreci(myArgs[1], myArgs[1]+myArgs[2], myArgs[3], myArgs[3]+myArgs[4], myArgs[5])
