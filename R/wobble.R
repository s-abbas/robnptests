## ----------------------------------------------------------------------------
## Add random noise if too many values in the samples are equal
## ----------------------------------------------------------------------------

wobble <- function(x, y, partially = FALSE) {
  ## Determine number of different values in both samples and in joint sample
  no.values.x <- length(unique(x))
  no.values.y <- length(unique(y))
  no.values.xy <- length(unique(c(x, y)))


  if (no.values.x == length(x) & no.values.y == length(y) & no.values.xy == length(c(x, y))) {
    ## If all values are distinct, return original observations
    return(list(x = x, y = y))
  } else if (!partially) {
    z <- c(x, y)

    ## Maximal number of digits after decimal point:
    # this expression always returns two (part before and part after the decimal
    # point)
    # max.digits <- max(sapply(strsplit(as.character(z), "\\."), length)) - 1
    # Instead:
    digits <- sapply(as.character(z),
                     function(x) nchar(unlist(strsplit(x, "\\."))[2]))
    digits[is.na(digits)] <- 0 # if the values are discrete we will get NA from the strsplit

    # I propose we use the full range of values between the observations:
    ## Add random noise from a uniform distribution
    z.wobble <- z + runif(length(z),
                          -0.5*10^(-max(digits)), 0.5*10^(-max(digits)))

    #z.wobble <- z + runif(length(z), min = -10^(-(max.digits + 1)), max = 10^(-(max.digits + 1)))
  } else {
    ind <- which(duplicated(c(x, y)))
    z.wobble <- c(x, y)
    z1 <- z.wobble[ind]

    digits <- sapply(as.character(z1), function(x) nchar(unlist(strsplit(x, "\\."))[2]))
    digits[is.na(digits)] <- 0

    z1.wobble <- z1 + runif(length(z1),
                          -0.5*10^(-max(digits)), 0.5*10^(-max(digits)))

    z.wobble[ind] <- z1.wobble
  }
    return(list(x = z.wobble[1:length(x)], y = z.wobble[(length(x) + 1):length(c(x, y))]))
}




