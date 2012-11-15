goodness.fit <- function(fdp, fda, starts, datas,
                method="SANN"){
  
  resultado = parametros = dados_ordenados = n = y = u <- NULL
  W_temp = A_temp = A_2 = W_2 <- NULL
  
  vero = function(par,x){
    -sum(log(fdp(par,x)))
  }
  
  if(method!="nlminb"){
    resultado = optim(starts,fn = vero, x=datas, control = list(maxit = 50000),
                      method=as.character(method))
  }
  if(method=="nlminb"){
    resultado = nlminb(starts,objective = vero, x=datas)
  }
  
  parametros = resultado$par
  dados_ordenados = sort(datas)
  v = fda(as.vector(parametros),dados_ordenados) # Dados ordenados.
  n = length(datas) # Tamanho da amostra.
  y = qnorm(v) # Inversa da acumulada da normal.
  u = pnorm((y-mean(y))/sqrt(var(y)))
  
  W_temp <- vector()
  A_temp <- vector()
  
  for(i in 1:n){
    W_temp[i] = (u[i] - (2*i-1)/(2*n))^2
    A_temp[i] = (2*i-1)*log(u[i]) + (2*n+1-2*i)*log(1-u[i])
  }
  
  A_2 = -n - mean(A_temp)
  W_2 = sum(W_temp) + 1/(12*n)
  W_estrela = W_2*(1+0.5/n)
  A_estrela = A_2*(1+0.75/n + 2.25/n^2)
  
  p = length(parametros)
  log.likelihood = -1*vero(parametros,datas) 
  AICc = -2*log.likelihood + 2*p + 2*(p*(p+1))/(n-p-1)
  AIC  = -2*log.likelihood + 2*p
  BIC  = -2*log.likelihood + p*log(n)
  resultado = (list("W" = W_estrela,"A" = A_estrela,
                    "EMV" = parametros, "AIC" = AIC ,"CAIC " = AICc, "BIC" = BIC))
  class(resultado) <- "list" 
  return(resultado)
  
} 
