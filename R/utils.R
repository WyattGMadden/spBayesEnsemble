
logit <- function(x) log(x / (1 - x))

ilogit <- function(x) 1 / (1 + exp(-x))

