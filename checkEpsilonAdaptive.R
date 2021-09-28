matplot(t(burned$samples$epsilon), type = "l", ylab = "epsilon")
samples <- burned$samples$alpha[1, ,]
accept <- apply(samples, 1, diff) != 0

meanacc <- matrix(0, nrow = nrow(accept), ncol = ncol(accept))
for (i in 1:nrow(accept)) {
  for(j in 1:ncol(accept)){
    meanacc[i,j]=     mean(accept[round(i/2) :i,j])
  }
}

matplot(meanacc, type = "l", ylab = "mean acceptance rate")
