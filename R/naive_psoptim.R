#' General-purpose optimization a naive Particle Swarm Optimization algorithm. Supports box-constrained optimization.
#'
#' @param string input character vector
#' @return numeric vector giving number of characters in each element of the
#'   character vector.  Missing strings have missing length.
#' @seealso \code{\link{nchar}} which this function wraps
#' @export
#' @examples
#' str_length(letters)
#' str_length(c("i", "like", "programming", NA))

naive_psoptim <- function(func,number_par,S,b_lo,b_up,lim_sup,lim_inf,N,data=NULL){
  
  dimension = number_par
  swarm_xi = swarm_pi = swarm_vi = matrix(NA,nrow=S,ncol=dimension)
  
  # Melhor posição do exame.
  g = runif(n=dimension,min=b_lo,max=b_up)
  
  # Função objetivo calculada em g.
  f_g = func(par=as.vector(g),x=as.vector(data))  

  if(as.character(f_g)=="NaN"){
    while(as.character(f_g)=="NaN"){
      g = runif(n=dimension,min=b_lo,max=b_up)
      f_g = func(par=g,x=as.vector(data))
    }
  }
  
  # AQUI COMECA A INICIALIZACAO DO ALGORITMO
  
  for(i in 1:S){
    # Inicializando a posição inicial de cada particula.
    x_i = runif(n=dimension,min=b_lo,max=b_up)
    
    # Inicializando a melhor posição da particula i à posição inicial. 
    swarm_pi[i,] = swarm_xi[i,] = x_i
      
    f_pi = func(par=x_i,x=as.vector(data))
    
    if(as.character(f_pi)=="NaN"){
      while(as.character(f_pi)=="NaN"){
        x_i = runif(n=dimension,min=b_lo,max=b_up)
        swarm_pi[i,] = swarm_xi[i,] = x_i
        f_pi = func(par=x_i,x=as.vector(data))
      }
    }
    
    if(as.character(f_pi)=="NaN") stop("f_pi igual à NaN para i igual à ",i)
    
    if(f_pi < f_g) g = as.vector(swarm_pi[i,])
    
    # Inicializando as valocidades das particulas
    swarm_vi[i,] = runif(n=dimension,min=-abs(b_up-b_lo),max=abs(b_up-b_lo)) 
  }
  
  # AQUI TERMINA A INICIALIZACAO DO ALGORITMO
    
  omega = 0.5
  phi_p = 0.5
  phi_g = 0.5
  
  is.integer0 <- function(x){
    is.integer(x) && length(x)==0L
  }
  
  m=1
  while(m<N){  
    for(i in 1:S){
      # r_p e r_g são numeros pseudo-aleatorios em (0,1).
      r_p = runif(n=dimension,min=0,max=1)
      r_g = runif(n=dimension,min=0,max=1)
      
      # Atualizando o vetor velocidade.
      swarm_vi[i,] = omega*swarm_vi[i,]+phi_p*r_p*(swarm_pi[i,]-swarm_xi[i,])+phi_g*r_g*(g-swarm_xi[i,])
      
      # Atualizando a posição de cada partícula.
      swarm_xi[i,] = swarm_xi[i,]+swarm_vi[i,]
      
      f_xi = func(par=as.vector(swarm_xi[i,]),x=as.vector(data))
      f_pi = func(par=as.vector(swarm_pi[i,]),x=as.vector(data))
      f_g =  func(par=g,x=as.vector(data))
      
      if(as.character(f_xi)=="NaN" || as.character(f_pi)=="NaN"){
        while(as.character(f_xi)=="NaN"){
          swarm_xi[i,] = runif(n=dimension,min=b_lo,max=b_up)
          f_xi = func(par=as.vector(swarm_xi[i,]),x=as.vector(data))
        }
        while(as.character(f_pi)=="NaN"){
          swarm_pi[i,] = runif(n=dimension,min=b_lo,max=b_up)
          f_pi = func(par=as.vector(swarm_pi[i,]),x=as.vector(data))
        }
      }

      # Existem valores abaixo do limite inferior de restrições?
      id_test_inf = which(swarm_xi[i,]<lim_inf)
      id_test_sup = which(swarm_xi[i,]>lim_sup)
      
      if(is.integer0(id_test_inf)!=TRUE){
        for(k in 1:length(id_test_inf)){
          new_value_inf = runif(1,min=lim_inf[id_test_inf[k]],max=lim_sup[id_test_inf[k]])
          swarm_xi[i,id_test_inf[k]] = new_value_inf
          swarm_pi[i,id_test_inf[k]] = new_value_inf
        }
      }
      
      if(is.integer0(id_test_sup)!=TRUE){
        for(k in 1:length(id_test_sup)){
          new_value_sup = runif(1,min=lim_inf[id_test_sup[k]],max=lim_sup[id_test_sup[k]])
          swarm_xi[i,id_test_sup[k]] = new_value_sup
          swarm_pi[i,id_test_sup[k]] = new_value_sup
        }
      }
   
      if(as.character(f_xi)=="NaN") stop("Fora da inicialização: f_xi igual à NaN para i igual à ",i)
      if(as.character(f_pi)=="NaN") stop("Fora da inicialização: f_pi igual à NaN para i igual à ",i)
      if(as.character(f_g)=="NaN") stop("Fora da inicialização: g igual à NaN para i igual à ",i)
      
      if(f_xi<f_pi){
        swarm_pi[i,] = swarm_xi[i,]
        if(f_pi<f_g) g = swarm_pi[i,] 
      } # Aqui termina o bloco if.
    } # Aqui termina o bloco for.
    m = m+1
  } # Aqui termina o bloco while.
  g
} # Aqui termina a função.