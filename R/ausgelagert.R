win_var <- function(x, gamma = 0, type = c("trim", "q2", "sk2", "sk5"), na.rm = FALSE) {
  type <- match.arg(type)
  
  ## Error handling
  if (gamma < 0 || gamma > 0.5) {
    stop ("gamma has to be in [0, 0.5]")
  }
  
  ## NA handling
  if (!na.rm & any(is.na(x))) {
    return(NA)
  } else if (na.rm & any(is.na(x))) {
    x <- stats::na.omit(x)
  }
  
  ## Calculate winsorized variance
  n <- length(x)
  
  switch (type,
          trim = {
            r <- floor(gamma * n)
            
            x.sort <- sort(x)
            x.lower <- x.sort[r + 1]
            x.upper <- x.sort[n - r]
            x.sort[x.sort < x.lower] <- x.lower
            x.sort[x.sort > x.upper] <- x.upper
            
            res <- 1 / (n - 1) * sum((x.sort - mean(x.sort)) ^ 2)
            h <- n - 2 * r
          },
          q2 = {
            n <- length(x)
            if (floor(0.05 * n) == 0) {
              stats::var(x)
            } else {
              x.sort <- sort(x)
              
              U.005 <-
                mean(x.sort[(n - floor(0.05 * n) + 1):n])
              L.005 <- mean(x.sort[1:floor(0.05 * n)])
              U.05 <-
                mean(x.sort[(n - floor(0.5 * n) + 1):n])
              L.05 <- mean(x.sort[1:floor(0.5 * n)])
              
              gamma.lower <-
                0.1 * (U.005 - L.005) / (U.005 - L.005 + U.05 - L.05)
              gamma.upper <- 0.1 - gamma.lower
              
              r.lower <- 1 + floor(gamma.lower * n + 0.5)
              r.upper <- n - floor(gamma.upper * n + 0.5)
              
              x.lower <- x.sort[r.lower]
              x.upper <- x.sort[r.upper]
              x.sort[x.sort < x.lower] <- x.lower
              x.sort[x.sort > x.upper] <- x.upper
              
              res <- 1 / (r.upper - r.lower + 1) * sum((x.sort - trim_mean(x, type = "q2")) ^
                                                         2)
              h <- r.upper - r.lower + 1
            }
          },
          sk2 = {
            x.sort <- sort(x)
            
            gamma.lower <-
              0.1 * (min(x) - stats::median(x)) / (min(x) - max(x))
            gamma.upper <- 0.1 - gamma.lower
            
            r.lower <- 1 + floor(gamma.lower * n + 0.5)
            r.upper <- n - floor(gamma.upper * n + 0.5)
            
            x.lower <- x.sort[r.lower]
            x.upper <- x.sort[r.upper]
            x.sort[x.sort < x.lower] <- x.lower
            x.sort[x.sort > x.upper] <- x.upper
            
            res <- 1 / (r.upper - r.lower + 1) * sum((x.sort - trim_mean(x, type = "sk2")) ^
                                                       2)
            h <- r.upper - r.lower + 1
          },
          sk5 = {
            gamma.lower <- 0.1 * (min(x) - stats::median(x)) / (min(x) - max(x))
            gamma.upper <- 0.1 - gamma.lower
            
            r.lower <- 1 + floor(gamma.lower * n + 0.5)
            r.upper <- n - floor(gamma.upper * n + 0.5)
            
            x.sort <- sort(x)
            x.lower <- x.sort[r.lower]
            x.upper <- x.sort[r.upper]
            x.sort[x.sort < x.lower] <- x.lower
            x.sort[x.sort > x.upper] <- x.upper
            
            res <- 1 / (r.upper - r.lower + 1) * sum((x.sort - trim_mean(x, type = "sk5")) ^
                                                       2)
            h <- r.upper - r.lower + 1
          }
  )
  
  return(list(var = res, h = h))
}


trim_mean <- function(x, gamma = 0.2, type = c("trim", "q2", "sk2", "sk5"), na.rm = FALSE) {
  type <- match.arg(type)
  
  ## Error handling
  if (gamma < 0 || gamma > 0.5) {
    stop ("gamma has to be in [0, 0.5]")
  }
  
  ## NA handling
  if (!na.rm & any(is.na(x))) {
    return(NA)
  } else if (na.rm & any(is.na(x))) {
    x <- stats::na.omit(x)
  }
  
  ## Calculate trimmed mean
  res <- switch (type,
                 trim = {
                   mean(x, trim = gamma) ## BB in R integrierte Funktion die einfach das arithmetische getrimmte Mittel berechnet
                 },
                 q2 = {
                   n <- length(x)
                   
                   if (floor(0.05 * n) == 0) {
                     ## No observations trimmed
                     mean(x)
                   } else {
                     x.sort <- sort(x)
                     
                     U.005 <- mean(x.sort[(n - floor(0.05 * n) + 1):n]) ## BB Woher kommen hier die 0.05? ist das 0.10 durch zwei? 
                     L.005 <- mean(x.sort[1:floor(0.05 * n)])           ## also der entprechende Trimmanteil? Nein, einfach das 
                     U.05 <- mean(x.sort[(n - floor(0.5 * n) + 1):n])   ## ensprechende Quantil!
                     L.05 <- mean(x.sort[1:floor(0.5 * n)])
                     
                     gamma.lower <- 0.1 * (U.005 - L.005)/(U.005 - L.005 + U.05 - L.05)
                     gamma.upper <- 0.1 - gamma.lower
                     
                     r.lower <- floor(n * gamma.lower) + 1
                     r.upper <- n - floor(n * gamma.upper)
                     
                     mean(x.sort[r.lower:r.upper])
                   }
                 },
                 sk2 = {
                   gamma.lower <- 0.1 * (min(x) - stats::median(x))/(min(x) - max(x))
                   gamma.upper <- 0.1 - gamma.lower
                   
                   n <- length(x)
                   r.lower <- floor(n * gamma.lower) + 1
                   r.upper <- n - floor(n * gamma.upper)
                   
                   x.sort <- sort(x)
                   mean(x.sort[r.lower:r.upper])
                 },
                 sk5 = {
                   gamma.lower <- 0.25 * (min(x) - mean(x))/(min(x) - max(x))
                   gamma.upper <- 0.25 - gamma.lower
                   
                   n <- length(x)
                   r.lower <- floor(n * gamma.lower) + 1
                   r.upper <- n - floor(n * gamma.upper)
                   
                   x.sort <- sort(x)
                   res <- mean(x.sort[r.lower:r.upper])
                 }
  )
  
  return(res)
}
