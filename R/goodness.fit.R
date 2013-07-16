goodness.fit <- function(fdp, fda, starts, datas, method = "L-BFGS-B", domain = c(0,Inf),
                         emv = NULL){
  
  if(missingArg(fda)==TRUE) stop("Unknown cumulative distribution function. The function needs to be informed.")
  if(missingArg(fdp)==TRUE) stop("Unknown probability density function. The function needs to be informed.")
  if(class(fdp)!="function") stop("The argument fdp must be a function. See the example in the documentation!")
  if(class(fda)!="function") stop("The argument fda must be a function. See the example in the documentation!")
  if(missingArg(datas)==TRUE) stop("Database missing!")
  if(TRUE%in%is.nan(datas)==TRUE) warning("The data have missing information!")
  if(length(domain)!=2) stop("The domain must have two arguments!")
  
  if(is.null(emv)==TRUE){
    if(missingArg(starts)==TRUE) stop("The initial shots were not informed!")
  }else{
    starts = emv 
  }   
  
  # Verifying properties of cumulative distribution function.
  
  if(fda(par=starts, x = domain[2])!=1) warning("The fda function informed is not a cumulative distribution function! The function no takes value 1 in Inf.")
  if(fda(par=starts, x = domain[1])!=0) warning("Check if the cumulative distribution informed is actually a distribution function.")
  
  myintegrate = function(...) tryCatch(integrate(...), error=function(e) NA)
  
  value_int = as.numeric(myintegrate(f=fdp,par=starts,lower=domain[1],
                                     upper=domain[2])[1])
  if(isTRUE(is.na(value_int))==TRUE) warning("Make sure that fdp is a probability density function. The integral in the domain specified is not 					     convergent.")
  
  if(isTRUE(is.na(value_int))!=TRUE){
     #Verifying properties of probability density function.
    if(value_int<0.99){
      warning("The integral from ", domain[1], " to ", domain[2]," of the probability density function has different from 1. Make sure the option 		      domain is correct.")
    }
    
    if(round(value_int)!=1){
      warning("fdp is not a probability density function.")
    }
  }
  
  if(is.null(emv)==TRUE){    
    vero = function(par,x){
      -sum(log(fdp(par,x)))
    }
    
    if(method == "Nelder-Mead" || method == "N"){
      resultado = optim(par = starts, fn = vero, x = datas,
                        method = "Nelder-Mead", hessian = TRUE)
    }
    
    if(method == "CG" || method == "C"){
      resultado = optim(par = starts, fn = vero, x = datas,
                        method = "CG", hessian = TRUE)
    }  
    
    if(method == "L-BFGS-B" || method == "L"){
      resultado = optim(par=starts, fn = vero,method="L-BFGS-B", x = datas,
                        lower=c(1e-10,1e-10,1e-10,1e-10,1e-10), upper=c(Inf,Inf,Inf,Inf,Inf), hessian=TRUE)
    }
    
    if(method == "SANN" || method == "S"){
      resultado = optim(par = starts, fn = vero, x = datas,
                        method = "SANN", hessian = TRUE)
    }  
    
    if(method == "BFGS" || method == "B"){
      resultado = optim(par = starts, fn = vero, x = datas,
                        method = "BFGS", hessian = TRUE)
    }
    
    if((FALSE %in% (method != c("L-BFGS-B", "L", "BFGS", "B",
                                "Nelder-Mead", "N", "SANN", "S", "CG", "C")))==FALSE){
      stop("Valid options are: L-BFGS-B or L, BFGS or B, Nelder-Mead or N, SANN or S, CG or C.")
    }
    
    parametros = resultado$par
    hessiana = resultado$hessian # matriz hessiana.
    
    dados_ordenados = sort(datas)
    v = fda(as.vector(parametros), dados_ordenados) # Dados ordenados.
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
    ks.testg = function(...) tryCatch(ks.test(...),
                                      warning = function(war) NA)
    KS = ks.testg(x = datas, y= "fda", par = as.vector(parametros))
    
    resultado = (list("W" = W_estrela,"A" = A_estrela, "KS" = KS,
                      "EMV" = parametros, "AIC" = AIC ,"CAIC " = AICc,
                      "BIC" = BIC, "Erro" = sqrt(diag(solve(hessiana))),
                      "Value" = resultado$value, "Convergence" = resultado$convergence))
    class(resultado) <- "list" 
    return(resultado)
  }
  
  if(class(domain)=="numeric"){
    
    vero = function(par,x){
      -sum(log(fdp(par,x)))
    }
    
    parametros = emv
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
    ks.testg = function(...) tryCatch(ks.test(...),
                                      warning = function(war) NA)
    KS = ks.testg(x = datas, y= "fda", par = as.vector(parametros))
    
    resultado = (list("W" = W_estrela,"A" = A_estrela, "KS" = KS, "AIC" = AIC,
                      "CAIC " = AICc, "BIC" = BIC))
    class(resultado) <- "list" 
    return(resultado)
  }
}
