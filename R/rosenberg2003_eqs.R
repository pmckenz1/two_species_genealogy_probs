# equation 1 in the Rosenberg 2003:
# (note the dependency on two short functions
# apk and abk, which are defined below)
g <- function(n, j, T_) {
  sum_gprobs <- 0
  
  # start of the summation
  for (k in j:n) {
    prod1 <- exp(-k*(k-1)*T_/2)
    prod2 <- (2*k-1)
    prod3 <- (-1)**(k-j)
    prod4 <- apk(j,k-1)
    prod5 <- abk(n,k)
    prod6 <- factorial(j)
    prod7 <- factorial(k-j)
    prod8 <- apk(n,k)
    
    val <- prod1*prod2*prod3*prod4*prod5/prod6/prod7/prod8
    
    sum_gprobs <- sum_gprobs + val
  }
  return(sum_gprobs)
}


# described under equation 1, and used in it:
apk <- function(a,k) {
  if (k == 0) {
    return(1)
  }
  mults <- prod(a+0:(k-1))
  return(mults)
}


# described under equation 1, and used in it:
abk <- function(a,k) {
  if (k == 0) {
    return(1)
  }
  mults <- prod(a-0:(k-1))
  return(mults)
}


# equation 14, prob of reciprocal monophyly given
# time (coalescent units) and sample sizes
PrC1 <- function(ra, rb, coalT_A, coalT_B) {
  sum_prob <- 0
  for (qa in 1:ra) {
    for (qb in 1:rb) {
      prod1 <-  g(ra,qa,coalT_A)*g(rb,qb,coalT_B)
      prod2 <- 2 / (choose(qa+qb, qa) * (qa+qb-1))
      sum_prob <- sum_prob + prod1*prod2
    }
  }
  return(sum_prob)
}


# equation 15, prob of monophyletic A and paraphyletic B given
# time (coalescent units) and sample sizes
PrC2 <- function(ra, rb, coalT_A, coalT_B) {
  sum_prob <- 0
  for (qa in 1:ra) {
    for (qb in 1:rb) {
      prod1 <-  g(ra,qa,coalT_A)*g(rb,qb,coalT_B)
      prod2 <- 2 / choose(qa+qb, qa)
      prod3 <- ((2*qa + qb) * (qb - 1)) / (qa * (qa + 1) * (qa + qb - 1))
      sum_prob <- sum_prob + prod1*prod2*prod3
    }
  }
  return(sum_prob)
}


# equation 16, prob of paraphyletic A and monophyletic B given
# time (coalescent units) and sample sizes
PrC3 <- function(ra, rb, coalT_A, coalT_B) {
  sum_prob <- 0
  for (qa in 1:ra) {
    for (qb in 1:rb) {
      prod1 <-  g(ra,qa,coalT_A)*g(rb,qb,coalT_B)
      prod2 <- 2 / choose(qa+qb, qa)
      prod3 <- ((qa + 2*qb) * (qa - 1)) / (qb * (qb + 1) * (qa + qb -1))
      sum_prob <- sum_prob + prod1*prod2*prod3
    }
  }
  return(sum_prob)
}

# equation 17, prob of polyphyly given
# time (coalescent units) and sample sizes
PrC4 <- function(ra, rb, coalT_A, coalT_B) {
  sum_prob <- 0
  for (qa in 1:ra) {
    for (qb in 1:rb) {
      prod1 <-  g(ra,qa,coalT_A)*g(rb,qb,coalT_B)
      prod2 <- 2 / choose(qa+qb, qa)
      prod3 <- ((qa + qb) / (qa * (qa + 1))) + ((qa + qb) / (qb * (qb + 1))) - (1 / (qa + qb - 1))
      sum_prob <- sum_prob + prod1*(1-prod2*prod3)
    }
  }
  return(sum_prob)
}

