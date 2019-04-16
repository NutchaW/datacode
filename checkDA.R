checkDA <- function(x) {
  N <- 500000
  
  ug <- max(x[1, ]) # Corrects if S>N
  if (ug > N) { x[1, x[1, ] > N] <- N - 1 }
  
  ug <- max(x[2, ]) # Corrects if I>N
  if (ug > N) { x[2, x[2, ] > N] <- median(x[2, ]) }
  
  ug <- min(min(x)) # Corrects if any state or parameter nudges negative
  if (ug <= 0) {
    for (i in 1:7) { x[i, x[i, ] < 0] <- mean(x[i, ]) }
  }
  
  ug <- min(x[6, ]) # Corrects if L < 200 days
  if (ug < 200) { x[6, x[6, ] < 200] <- median(x[6, ]) }
  
  ug <- min(x[7, ]) # Corrects if D < .5
  if (ug <= .5) { x[7, x[7, ] < 0.5] <- median(x[7, ]) }
  
  ug <- min(x[4, ] - x[5, ]) # Corrects if R0mx <= R0mn
  if (ug <= 0) { x[4, x[4, ] < x[5, ]] <- x[5, x[4, ] < x[5, ]] + 0.01 }
  
  return(x)
}