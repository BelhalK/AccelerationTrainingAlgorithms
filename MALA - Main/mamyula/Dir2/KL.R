KL <- function(p, q, seuil=1e-10 ) {
    step = 1/sum(p) 
    use <- (p>seuil) & (q>seuil) #pour éviter des pbms numériques 
    k <- sum( p[use] * log(p[use]/q[use]) )
    return( k * step )
    }