# Likelihood function
ll.fun <- function(a, b, tau){
  ll <- sum(dnorm(x = y, mean = a + b * x, sd = sqrt(1 / tau),log = TRUE))
  return(ll)
}

# prior function (uninformative)
prior.fun <- function(a, b, tau){
  lp <- dnorm(x = a, mean = 0, sd = sqrt(1/10^-6), log = TRUE) +
    dnorm(x = b, mean = 0, sd = sqrt(1/10^-6), log = TRUE) +
    dgamma(x = tau, shape = 0.001, rate = 0.001, log = TRUE)
  return(lp)
}

# MCMC - metropolis algorithm

# N = lengte van de MCMC
# sigma = grootte van de stappen in de markov chain
# init_a en init_b = initiële waarden voor a en b

MCMC_metro <- function(x, y, N = 20, sigma = 1, init_a, init_b){

  set.seed(56)
  # We maken een lege dataframe om de output in te bewaren
  mcmc <- data.frame(iter = 1:N,
                     a = NA,
                     b = NA,
                     tau = NA,
                     posterior = NA,
                     accept = NA)

  mcmc$a[1] <- init_a
  mcmc$b[1] <- init_b
  mcmc$tau[1] <- 0.228  # init_tau

  # Berekening van de posterior voor de eerste lijn (met de startwaarde)
  mcmc$posterior[1] <- ll.fun(a = mcmc$a[1], b = mcmc$b[1], tau = mcmc$tau[1]) +
    prior.fun(a = mcmc$a[1], b = mcmc$b[1], tau = mcmc$tau[1])
  mcmc$accept[1] <- 1   # De initiële waarde wordt altijd geaccepteerd

  for(i in 2:N){

    # Kies een nieuwe waarde uit een normale verdeling
    # met gemiddelde = vorige waarde en standaard afwijking = sigma (stapgrootte)
    new_a <- rnorm(1, mean = mcmc$a[i - 1], sd = sigma)
    new_b <- rnorm(1, mean = mcmc$b[i - 1], sd = sigma)
    #new_tau <- rnorm(1, mean = mcmc$tau[i - 1], sd = sigma)
    new_tau <- mcmc$tau[i -1]

    # Bereken de posterior
    new_post <- ll.fun(a = new_a, b = new_b, tau = new_tau) +
      prior.fun(a = new_a, b = new_b, tau = new_tau)

    # Bereken de log ratio  (in de log schaal is dit het verschil)
    log.ratio <- new_post - mcmc$post[i-1]

    # Als de log.ratio > 0 => accepteren
    if (log.ratio > 0) {
      mcmc$accept[i] <- 1
      # Als log.ratio < 0 is de kans op accepteren afhankelijk van de log.ratio
    }else{
      if (rbinom(n = 1, size = 1, prob = exp(log.ratio))) {
        mcmc$accept[i] <- 1
      }else{
        mcmc$accept[i] <- 0
      }
    }

    # De nieuwe waarde wordt geaccepteerd
    if(mcmc$accept[i]) {
      mcmc$a[i] <- new_a
      mcmc$b[i] <- new_b
      mcmc$tau[i] <- new_tau
      mcmc$posterior[i] <- new_post
    }else{
      # De nieuwe waarde wordt niet geaccepteerd (de oude blijft behouden)
      mcmc$a[i] <- mcmc$a[i-1]
      mcmc$b[i] <- mcmc$b[i-1]
      mcmc$tau[i] <- mcmc$tau[i-1]
      mcmc$posterior[i] <- mcmc$posterior[i-1]
    }
  }
  return(mcmc = mcmc)
}
