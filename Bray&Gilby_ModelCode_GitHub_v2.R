#Model Code

#Title: Social relationships among adult male chimpanzees (Pan troglodytes schweinfurthii): Variation in the strength and quality of social bonds
#Authors: Joel Bray & Ian C. Gilby

############################################################
#BOND STRENGTH MODELS
############################################################

#
#Dyadic Association Index (DAI)
#

m.strength.dai  <- map2stan(
  alist(
    
    DAI ~ dzagamma2( p, mu , scale ),
    
    logit(p) ~ ap + ap_id[id1] + ap_id[id2] + ap_dyad[dyad] + ap_year[year] +
      (bp_rank + bp_rank_id[id1])*rank_id1 + (bp_rank + bp_rank_id[id2])*rank_id2 +
      (bp_age + bp_age_id[id1])*age_id1 + (bp_age + bp_age_id[id2])*age_id2 +
      (bp_rankDiff + bp_rankDiff_dyad[dyad])*rank_diff +
      bp_ageDiff*age_diff +
      bp_kin*maternal_kin + 
      bp_kinXageDiff*maternal_kin*age_diff +
      bp_kinXrankDiff*maternal_kin*rank_diff,
    
    log(mu) ~ am + am_id[id1] + am_id[id2] + am_dyad[dyad] + am_year[year] +
      (bm_rank + bp_rank_id[id1])*rank_id1 + (bm_rank + bm_rank_id[id2])*rank_id2 +
      (bm_age + bp_age_id[id1])*age_id1 + (bm_age + bm_age_id[id2])*age_id2 +
      (bm_rankDiff + bm_rankDiff_dyad[dyad])*rank_diff +
      bm_ageDiff*age_diff +
      bm_kin*maternal_kin + 
      bm_kinXageDiff*maternal_kin*age_diff +
      bm_kinXrankDiff*maternal_kin*rank_diff,
    
    c(ap_year, am_year)[year] ~ dmvnormNC( sigma_year, Rho_year ),
    c(ap_id, am_id, bp_rank_id, bm_rank_id, bp_age_id, bm_age_id)[id1] ~ dmvnormNC( sigma_id, Rho_id ),
    c(ap_dyad, am_dyad, bp_rankDiff_dyad, bm_rankDiff_dyad)[dyad] ~ dmvnormNC( sigma_dyad, Rho_dyad ),
    ap ~ dnorm(0,2),
    am ~ dnorm(0,2),
    c(bp_rank, bm_rank, bp_age, bm_age, bp_kin, bm_kin, bp_rankDiff, bm_rankDiff, bp_ageDiff, bm_ageDiff, bp_kinXageDiff, bm_kinXageDiff, bp_kinXrankDiff, bm_kinXrankDiff) ~ dnorm(0,2),
    c(sigma_year, sigma_id, sigma_dyad) ~ dexp(1), 
    c(Rho_year, Rho_id, Rho_dyad) ~ dlkjcorr(3),
    scale ~ dexp(1)
    
  ),
  
  data=data.strength.dai, cores=3 , chains=3 , warmup=1000, iter=4000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE
  
)



#
#Dyadic Grooming Index (DGI)
#

m.strength.dgi  <- map2stan(
  alist(
    
    DGI ~ dzagamma2( p, mu , scale ),
    
    logit(p) ~ ap + ap_id[id1] + ap_id[id2] + ap_dyad[dyad] + ap_year[year] +
      (bp_rank + bp_rank_id[id1])*rank_id1 + (bp_rank + bp_rank_id[id2])*rank_id2 +
      (bp_age + bp_age_id[id1])*age_id1 + (bp_age + bp_age_id[id2])*age_id2 +
      (bp_rankDiff + bp_rankDiff_dyad[dyad])*rank_diff +
      bp_ageDiff*age_diff +
      bp_kin*maternal_kin + 
      bp_kinXageDiff*maternal_kin*age_diff +
      bp_kinXrankDiff*maternal_kin*rank_diff,
    
    log(mu) ~ am + am_id[id1] + am_id[id2] + am_dyad[dyad] + am_year[year] +
      (bm_rank + bp_rank_id[id1])*rank_id1 + (bm_rank + bm_rank_id[id2])*rank_id2 +
      (bm_age + bp_age_id[id1])*age_id1 + (bm_age + bm_age_id[id2])*age_id2 +
      (bm_rankDiff + bm_rankDiff_dyad[dyad])*rank_diff +
      bm_ageDiff*age_diff +
      bm_kin*maternal_kin + 
      bm_kinXageDiff*maternal_kin*age_diff +
      bm_kinXrankDiff*maternal_kin*rank_diff,
    
    c(ap_year, am_year)[year] ~ dmvnormNC( sigma_year, Rho_year ),
    c(ap_id, am_id, bp_rank_id, bm_rank_id, bp_age_id, bm_age_id)[id1] ~ dmvnormNC( sigma_id, Rho_id ),
    c(ap_dyad, am_dyad, bp_rankDiff_dyad, bm_rankDiff_dyad)[dyad] ~ dmvnormNC( sigma_dyad, Rho_dyad ),
    ap ~ dnorm(0,2),
    am ~ dnorm(0,2),
    c(bp_rank, bm_rank, bp_age, bm_age, bp_kin, bm_kin, bp_rankDiff, bm_rankDiff, bp_ageDiff, bm_ageDiff, bp_kinXageDiff, bm_kinXageDiff, bp_kinXrankDiff, bm_kinXrankDiff) ~ dnorm(0,2),
    c(sigma_year, sigma_id, sigma_dyad) ~ dexp(1), 
    c(Rho_year, Rho_id, Rho_dyad) ~ dlkjcorr(3),
    scale ~ dexp(1)
    
  ),
  
  data=data.strength.dgi, cores=3 , chains=3 , warmup=1000, iter=4000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE
  
)





############################################################
#GROOMING EQUALITY MODELS
############################################################

#
#Model with dyadic association index
#

m.eq.index.dai <- map2stan(
  alist(
    
    eq_index ~ dbeta2( p , theta ),
    
    logit(p) ~ ap + ap_id[id1] + ap_id[id2] + ap_dyad[dyad] + tau*z[year] +
      (bp_rank + bp_rank_id[id1])*rank_id1 + (bp_rank + bp_rank_id[id2])*rank_id2 + 
      (bp_age + bp_age_id[id1])*age_id1 + (bp_age + bp_age_id[id2])*age_id2 + 
      (bp_rankDiff + bp_rankDiff_dyad[dyad])*rank_diff +
      bp_ageDiff*age_diff +
      bp_kin*maternal_kin + 
      bp_dai*DAI +
      bp_male_comm_size*male_comm_size,
    
    c(ap_id, bp_age_id, bp_rank_id )[id1] ~ dmvnormNC( sigma_id, Rho_id ),
    c(ap_dyad, bp_rankDiff_dyad )[dyad] ~ dmvnormNC( sigma_dyad, Rho_dyad ),
    c(ap) ~ dnorm(0,2), 
    c(bp_rank, bp_age, bp_rankDiff, bp_ageDiff, bp_kin, bp_dai, bp_male_comm_size) ~ dnorm(0,2), 
    c(sigma_id, sigma_dyad) ~ dexp(1),
    c(Rho_id, Rho_dyad) ~ dlkjcorr(3),
    z[year] ~ dnorm(0,1),
    tau ~ dcauchy(0,1),
    theta ~ dexp(1)
    
  ),
  
  data=data.equality, cores=3 , chains=3 , warmup=1000, iter=4000, control=list(adapt_delta=0.99), start=list(theta = 3), constraints=list(tau="lower=0"), WAIC=TRUE
  
)


#
#Model with dyadic grooming index
#

m.eq.index.dgi <- map2stan(
  alist(
    
    eq_index ~ dbeta2( p , theta ),
    
    logit(p) ~ ap + ap_id[id1] + ap_id[id2] + ap_dyad[dyad] + tau*z[year] +
      (bp_rank + bp_rank_id[id1])*rank_id1 + (bp_rank + bp_rank_id[id2])*rank_id2 + 
      (bp_age + bp_age_id[id1])*age_id1 + (bp_age + bp_age_id[id2])*age_id2 + 
      (bp_rankDiff + bp_rankDiff_dyad[dyad])*rank_diff +
      bp_ageDiff*age_diff +
      bp_kin*maternal_kin + 
      bp_dgi*DGI +
      bp_male_comm_size*male_comm_size,
    
    c(ap_id, bp_age_id, bp_rank_id )[id1] ~ dmvnormNC( sigma_id, Rho_id ),
    c(ap_dyad, bp_rankDiff_dyad )[dyad] ~ dmvnormNC( sigma_dyad, Rho_dyad ),
    c(ap) ~ dnorm(0,2), 
    c(bp_rank, bp_age, bp_rankDiff, bp_ageDiff, bp_kin, bp_dgi, bp_male_comm_size) ~ dnorm(0,2), 
    c(sigma_id, sigma_dyad) ~ dexp(1),
    c(Rho_id, Rho_dyad) ~ dlkjcorr(3),
    z[year] ~ dnorm(0,1),
    tau ~ dcauchy(0,1),
    theta ~ dexp(1)
    
  ),
  
  data=data.equality, cores=3 , chains=3 , warmup=1000, iter=4000, control=list(adapt_delta=0.99), start=list(theta = 3), constraints=list(tau="lower=0"), WAIC=TRUE
  
)


#
#Model with no bond strength index
#

m.eq.index.slim <- map2stan(
  alist(
    
    eq_index ~ dbeta2( p , theta ),
    
    logit(p) ~ ap + ap_id[id1] + ap_id[id2] + ap_dyad[dyad] + tau*z[year] +
      (bp_rank + bp_rank_id[id1])*rank_id1 + (bp_rank + bp_rank_id[id2])*rank_id2 + 
      (bp_age + bp_age_id[id1])*age_id1 + (bp_age + bp_age_id[id2])*age_id2 + 
      (bp_rankDiff + bp_rankDiff_dyad[dyad])*rank_diff +
      bp_ageDiff*age_diff +
      bp_kin*maternal_kin + 
      bp_male_comm_size*male_comm_size,
    
    c(ap_id, bp_age_id, bp_rank_id )[id1] ~ dmvnormNC( sigma_id, Rho_id ),
    c(ap_dyad, bp_rankDiff_dyad )[dyad] ~ dmvnormNC( sigma_dyad, Rho_dyad ),
    c(ap) ~ dnorm(0,2), 
    c(bp_rank, bp_age, bp_rankDiff, bp_ageDiff, bp_kin, bp_male_comm_size) ~ dnorm(0,2), 
    c(sigma_id, sigma_dyad) ~ dexp(1),
    c(Rho_id, Rho_dyad) ~ dlkjcorr(3),
    z[year] ~ dnorm(0,1),
    tau ~ dcauchy(0,1),
    theta ~ dexp(1)
    
  ),
  
  data=data.equality, cores=1, chains=1 , warmup=1000, iter=4000, control=list(adapt_delta=0.99), start=list(theta = 3), constraints=list(tau="lower=0"), WAIC=TRUE
  
)





############################################################
#BOND STABILITY MODELS
############################################################

#
#Dyadic Association Index
#

m.prior1.dai <- map2stan(
  alist(
    
    DAI ~ dzagamma2( p, mu , scale ),
    
    logit(p) ~ ap + ap_id[id1] + ap_id[id2] + ap_dyad[dyad] + ap_year[year] +
      (bp_rank + bp_rank_id[id1])*rank_id1 + (bp_rank + bp_rank_id[id2])*rank_id2 +
      (bp_age + bp_age_id[id1])*age_id1 + (bp_age + bp_age_id[id2])*age_id2 +
      (bp_rankDiff + bp_rankDiff_dyad[dyad])*rank_diff +
      bp_ageDiff*age_diff +
      bp_kin*maternal_kin + 
      bp_kinXageDiff*maternal_kin*age_diff +
      bp_kinXrankDiff*maternal_kin*rank_diff +
      bp_prior1*dai_prior1,
    
    log(mu) ~ am + am_id[id1] + am_id[id2] + am_dyad[dyad] + am_year[year] +
      (bm_rank + bp_rank_id[id1])*rank_id1 + (bm_rank + bm_rank_id[id2])*rank_id2 +
      (bm_age + bp_age_id[id1])*age_id1 + (bm_age + bm_age_id[id2])*age_id2 +
      (bm_rankDiff + bm_rankDiff_dyad[dyad])*rank_diff +
      bm_ageDiff*age_diff +
      bm_kin*maternal_kin + 
      bm_kinXageDiff*maternal_kin*age_diff +
      bm_kinXrankDiff*maternal_kin*rank_diff +
      bm_prior1*dai_prior1,
    
    c(ap_year, am_year)[year] ~ dmvnormNC( sigma_year, Rho_year ),
    c(ap_id, am_id, bp_rank_id, bm_rank_id, bp_age_id, bm_age_id)[id1] ~ dmvnormNC( sigma_id, Rho_id ),
    c(ap_dyad, am_dyad, bp_rankDiff_dyad, bm_rankDiff_dyad)[dyad] ~ dmvnormNC( sigma_dyad, Rho_dyad ),
    ap ~ dnorm(0,2),
    am ~ dnorm(0,2),
    c(bp_rank, bm_rank, bp_age, bm_age, bp_kin, bm_kin, bp_rankDiff, bm_rankDiff, bp_ageDiff, bm_ageDiff, bp_kinXageDiff, bm_kinXageDiff, bp_kinXrankDiff, bm_kinXrankDiff, bp_prior1, bm_prior1) ~ dnorm(0,2),
    c(sigma_year, sigma_id, sigma_dyad) ~ dexp(1), 
    c(Rho_year, Rho_id, Rho_dyad) ~ dlkjcorr(3),
    scale ~ dexp(1)
    
  ),
  
  data=data.prior1.dai, cores=3 , chains=3 , warmup=1000, iter=4000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE
  
)


m.prior2.dai <- map2stan(
  alist(
    
    DAI ~ dzagamma2( p, mu , scale ),
    
    logit(p) ~ ap + ap_id[id1] + ap_id[id2] + ap_dyad[dyad] + ap_year[year] +
      (bp_rank + bp_rank_id[id1])*rank_id1 + (bp_rank + bp_rank_id[id2])*rank_id2 +
      (bp_age + bp_age_id[id1])*age_id1 + (bp_age + bp_age_id[id2])*age_id2 +
      (bp_rankDiff + bp_rankDiff_dyad[dyad])*rank_diff +
      bp_ageDiff*age_diff +
      bp_kin*maternal_kin + 
      bp_kinXageDiff*maternal_kin*age_diff +
      bp_kinXrankDiff*maternal_kin*rank_diff +
      bp_prior2*dai_prior2,
    
    log(mu) ~ am + am_id[id1] + am_id[id2] + am_dyad[dyad] + am_year[year] +
      (bm_rank + bp_rank_id[id1])*rank_id1 + (bm_rank + bm_rank_id[id2])*rank_id2 +
      (bm_age + bp_age_id[id1])*age_id1 + (bm_age + bm_age_id[id2])*age_id2 +
      (bm_rankDiff + bm_rankDiff_dyad[dyad])*rank_diff +
      bm_ageDiff*age_diff +
      bm_kin*maternal_kin + 
      bm_kinXageDiff*maternal_kin*age_diff +
      bm_kinXrankDiff*maternal_kin*rank_diff +
      bm_prior2*dai_prior2,
    
    c(ap_year, am_year)[year] ~ dmvnormNC( sigma_year, Rho_year ),
    c(ap_id, am_id, bp_rank_id, bm_rank_id, bp_age_id, bm_age_id)[id1] ~ dmvnormNC( sigma_id, Rho_id ),
    c(ap_dyad, am_dyad, bp_rankDiff_dyad, bm_rankDiff_dyad)[dyad] ~ dmvnormNC( sigma_dyad, Rho_dyad ),
    ap ~ dnorm(0,2),
    am ~ dnorm(0,2),
    c(bp_rank, bm_rank, bp_age, bm_age, bp_kin, bm_kin, bp_rankDiff, bm_rankDiff, bp_ageDiff, bm_ageDiff, bp_kinXageDiff, bm_kinXageDiff, bp_kinXrankDiff, bm_kinXrankDiff, bp_prior2, bm_prior2) ~ dnorm(0,2),
    c(sigma_year, sigma_id, sigma_dyad) ~ dexp(1), 
    c(Rho_year, Rho_id, Rho_dyad) ~ dlkjcorr(3),
    scale ~ dexp(1)
    
  ),
  
  data=data.prior2.dai, cores=3 , chains=3 , warmup=1000, iter=4000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE
  
)


m.prior3.dai <- map2stan(
  alist(
    
    DAI ~ dzagamma2( p, mu , scale ),
    
    logit(p) ~ ap + ap_id[id1] + ap_id[id2] + ap_dyad[dyad] + ap_year[year] +
      (bp_rank + bp_rank_id[id1])*rank_id1 + (bp_rank + bp_rank_id[id2])*rank_id2 +
      (bp_age + bp_age_id[id1])*age_id1 + (bp_age + bp_age_id[id2])*age_id2 +
      (bp_rankDiff + bp_rankDiff_dyad[dyad])*rank_diff +
      bp_ageDiff*age_diff +
      bp_kin*maternal_kin + 
      bp_kinXageDiff*maternal_kin*age_diff +
      bp_kinXrankDiff*maternal_kin*rank_diff +
      bp_prior3*dai_prior3,
    
    log(mu) ~ am + am_id[id1] + am_id[id2] + am_dyad[dyad] + am_year[year] +
      (bm_rank + bp_rank_id[id1])*rank_id1 + (bm_rank + bm_rank_id[id2])*rank_id2 +
      (bm_age + bp_age_id[id1])*age_id1 + (bm_age + bm_age_id[id2])*age_id2 +
      (bm_rankDiff + bm_rankDiff_dyad[dyad])*rank_diff +
      bm_ageDiff*age_diff +
      bm_kin*maternal_kin + 
      bm_kinXageDiff*maternal_kin*age_diff +
      bm_kinXrankDiff*maternal_kin*rank_diff +
      bm_prior3*dai_prior3,
    
    c(ap_year, am_year)[year] ~ dmvnormNC( sigma_year, Rho_year ),
    c(ap_id, am_id, bp_rank_id, bm_rank_id, bp_age_id, bm_age_id)[id1] ~ dmvnormNC( sigma_id, Rho_id ),
    c(ap_dyad, am_dyad, bp_rankDiff_dyad, bm_rankDiff_dyad)[dyad] ~ dmvnormNC( sigma_dyad, Rho_dyad ),
    ap ~ dnorm(0,2),
    am ~ dnorm(0,2),
    c(bp_rank, bm_rank, bp_age, bm_age, bp_kin, bm_kin, bp_rankDiff, bm_rankDiff, bp_ageDiff, bm_ageDiff, bp_kinXageDiff, bm_kinXageDiff, bp_kinXrankDiff, bm_kinXrankDiff, bp_prior3, bm_prior3) ~ dnorm(0,2),
    c(sigma_year, sigma_id, sigma_dyad) ~ dexp(1), 
    c(Rho_year, Rho_id, Rho_dyad) ~ dlkjcorr(3),
    scale ~ dexp(1)
    
  ),
  
  data=data.prior3.dai, cores=3 , chains=3 , warmup=1000, iter=4000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE
  
)


m.prior4.dai <- map2stan(
  alist(
    
    DAI ~ dzagamma2( p, mu , scale ),
    
    logit(p) ~ ap + ap_id[id1] + ap_id[id2] + ap_dyad[dyad] + ap_year[year] +
      (bp_rank + bp_rank_id[id1])*rank_id1 + (bp_rank + bp_rank_id[id2])*rank_id2 +
      (bp_age + bp_age_id[id1])*age_id1 + (bp_age + bp_age_id[id2])*age_id2 +
      (bp_rankDiff + bp_rankDiff_dyad[dyad])*rank_diff +
      bp_ageDiff*age_diff +
      bp_kin*maternal_kin + 
      bp_kinXageDiff*maternal_kin*age_diff +
      bp_kinXrankDiff*maternal_kin*rank_diff +
      bp_prior4*dai_prior4,
    
    log(mu) ~ am + am_id[id1] + am_id[id2] + am_dyad[dyad] + am_year[year] +
      (bm_rank + bp_rank_id[id1])*rank_id1 + (bm_rank + bm_rank_id[id2])*rank_id2 +
      (bm_age + bp_age_id[id1])*age_id1 + (bm_age + bm_age_id[id2])*age_id2 +
      (bm_rankDiff + bm_rankDiff_dyad[dyad])*rank_diff +
      bm_ageDiff*age_diff +
      bm_kin*maternal_kin + 
      bm_kinXageDiff*maternal_kin*age_diff +
      bm_kinXrankDiff*maternal_kin*rank_diff +
      bm_prior4*dai_prior4,
    
    c(ap_year, am_year)[year] ~ dmvnormNC( sigma_year, Rho_year ),
    c(ap_id, am_id, bp_rank_id, bm_rank_id, bp_age_id, bm_age_id)[id1] ~ dmvnormNC( sigma_id, Rho_id ),
    c(ap_dyad, am_dyad, bp_rankDiff_dyad, bm_rankDiff_dyad)[dyad] ~ dmvnormNC( sigma_dyad, Rho_dyad ),
    ap ~ dnorm(0,2),
    am ~ dnorm(0,2),
    c(bp_rank, bm_rank, bp_age, bm_age, bp_kin, bm_kin, bp_rankDiff, bm_rankDiff, bp_ageDiff, bm_ageDiff, bp_kinXageDiff, bm_kinXageDiff, bp_kinXrankDiff, bm_kinXrankDiff, bp_prior4, bm_prior4) ~ dnorm(0,2),
    c(sigma_year, sigma_id, sigma_dyad) ~ dexp(1), 
    c(Rho_year, Rho_id, Rho_dyad) ~ dlkjcorr(3),
    scale ~ dexp(1)
    
  ),
  
  data=data.prior4.dai, cores=3 , chains=3 , warmup=1000, iter=4000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE
  
)


m.prior5.dai <- map2stan(
  alist(
    
    DAI ~ dzagamma2( p, mu , scale ),
    
    logit(p) ~ ap + ap_id[id1] + ap_id[id2] + ap_dyad[dyad] + ap_year[year] +
      (bp_rank + bp_rank_id[id1])*rank_id1 + (bp_rank + bp_rank_id[id2])*rank_id2 +
      (bp_age + bp_age_id[id1])*age_id1 + (bp_age + bp_age_id[id2])*age_id2 +
      (bp_rankDiff + bp_rankDiff_dyad[dyad])*rank_diff +
      bp_ageDiff*age_diff +
      bp_kin*maternal_kin + 
      bp_kinXageDiff*maternal_kin*age_diff +
      bp_kinXrankDiff*maternal_kin*rank_diff +
      bp_prior5*dai_prior5,
    
    log(mu) ~ am + am_id[id1] + am_id[id2] + am_dyad[dyad] + am_year[year] +
      (bm_rank + bp_rank_id[id1])*rank_id1 + (bm_rank + bm_rank_id[id2])*rank_id2 +
      (bm_age + bp_age_id[id1])*age_id1 + (bm_age + bm_age_id[id2])*age_id2 +
      (bm_rankDiff + bm_rankDiff_dyad[dyad])*rank_diff +
      bm_ageDiff*age_diff +
      bm_kin*maternal_kin + 
      bm_kinXageDiff*maternal_kin*age_diff +
      bm_kinXrankDiff*maternal_kin*rank_diff +
      bm_prior5*dai_prior5,
    
    c(ap_year, am_year)[year] ~ dmvnormNC( sigma_year, Rho_year ),
    c(ap_id, am_id, bp_rank_id, bm_rank_id, bp_age_id, bm_age_id)[id1] ~ dmvnormNC( sigma_id, Rho_id ),
    c(ap_dyad, am_dyad, bp_rankDiff_dyad, bm_rankDiff_dyad)[dyad] ~ dmvnormNC( sigma_dyad, Rho_dyad ),
    ap ~ dnorm(0,2),
    am ~ dnorm(0,2),
    c(bp_rank, bm_rank, bp_age, bm_age, bp_kin, bm_kin, bp_rankDiff, bm_rankDiff, bp_ageDiff, bm_ageDiff, bp_kinXageDiff, bm_kinXageDiff, bp_kinXrankDiff, bm_kinXrankDiff, bp_prior5, bm_prior5) ~ dnorm(0,2),
    c(sigma_year, sigma_id, sigma_dyad) ~ dexp(1), 
    c(Rho_year, Rho_id, Rho_dyad) ~ dlkjcorr(3),
    scale ~ dexp(1)
    
  ),
  
  data=data.prior5.dai, cores=3 , chains=3 , warmup=1000, iter=4000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE
  
)



#
#Dyadic Grooming Index
#

m.prior1.dgi <- map2stan(
  alist(
    
    DGI ~ dzagamma2( p, mu , scale ),
    
    logit(p) ~ ap + ap_id[id1] + ap_id[id2] + ap_dyad[dyad] + ap_year[year] +
      (bp_rank + bp_rank_id[id1])*rank_id1 + (bp_rank + bp_rank_id[id2])*rank_id2 +
      (bp_age + bp_age_id[id1])*age_id1 + (bp_age + bp_age_id[id2])*age_id2 +
      (bp_rankDiff + bp_rankDiff_dyad[dyad])*rank_diff +
      bp_ageDiff*age_diff +
      bp_kin*maternal_kin + 
      bp_kinXageDiff*maternal_kin*age_diff +
      bp_kinXrankDiff*maternal_kin*rank_diff +
      bp_prior1*dgi_prior1,
    
    log(mu) ~ am + am_id[id1] + am_id[id2] + am_dyad[dyad] + am_year[year] +
      (bm_rank + bp_rank_id[id1])*rank_id1 + (bm_rank + bm_rank_id[id2])*rank_id2 +
      (bm_age + bp_age_id[id1])*age_id1 + (bm_age + bm_age_id[id2])*age_id2 +
      (bm_rankDiff + bm_rankDiff_dyad[dyad])*rank_diff +
      bm_ageDiff*age_diff +
      bm_kin*maternal_kin + 
      bm_kinXageDiff*maternal_kin*age_diff +
      bm_kinXrankDiff*maternal_kin*rank_diff +
      bm_prior1*dgi_prior1,
    
    c(ap_year, am_year)[year] ~ dmvnormNC( sigma_year, Rho_year ),
    c(ap_id, am_id, bp_rank_id, bm_rank_id, bp_age_id, bm_age_id)[id1] ~ dmvnormNC( sigma_id, Rho_id ),
    c(ap_dyad, am_dyad, bp_rankDiff_dyad, bm_rankDiff_dyad)[dyad] ~ dmvnormNC( sigma_dyad, Rho_dyad ),
    ap ~ dnorm(0,2),
    am ~ dnorm(0,2),
    c(bp_rank, bm_rank, bp_age, bm_age, bp_kin, bm_kin, bp_rankDiff, bm_rankDiff, bp_ageDiff, bm_ageDiff, bp_kinXageDiff, bm_kinXageDiff, bp_kinXrankDiff, bm_kinXrankDiff, bp_prior1, bm_prior1) ~ dnorm(0,2),
    c(sigma_year, sigma_id, sigma_dyad) ~ dexp(1), 
    c(Rho_year, Rho_id, Rho_dyad) ~ dlkjcorr(3),
    scale ~ dexp(1)
    
  ),
  
  data=data.prior1.dgi, cores=3 , chains=3 , warmup=1000, iter=4000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE
  
)


m.prior2.dgi <- map2stan(
  alist(
    
    DGI ~ dzagamma2( p, mu , scale ),
    
    logit(p) ~ ap + ap_id[id1] + ap_id[id2] + ap_dyad[dyad] + ap_year[year] +
      (bp_rank + bp_rank_id[id1])*rank_id1 + (bp_rank + bp_rank_id[id2])*rank_id2 +
      (bp_age + bp_age_id[id1])*age_id1 + (bp_age + bp_age_id[id2])*age_id2 +
      (bp_rankDiff + bp_rankDiff_dyad[dyad])*rank_diff +
      bp_ageDiff*age_diff +
      bp_kin*maternal_kin + 
      bp_kinXageDiff*maternal_kin*age_diff +
      bp_kinXrankDiff*maternal_kin*rank_diff +
      bp_prior2*dgi_prior2,
    
    log(mu) ~ am + am_id[id1] + am_id[id2] + am_dyad[dyad] + am_year[year] +
      (bm_rank + bp_rank_id[id1])*rank_id1 + (bm_rank + bm_rank_id[id2])*rank_id2 +
      (bm_age + bp_age_id[id1])*age_id1 + (bm_age + bm_age_id[id2])*age_id2 +
      (bm_rankDiff + bm_rankDiff_dyad[dyad])*rank_diff +
      bm_ageDiff*age_diff +
      bm_kin*maternal_kin + 
      bm_kinXageDiff*maternal_kin*age_diff +
      bm_kinXrankDiff*maternal_kin*rank_diff +
      bm_prior2*dgi_prior2,
    
    c(ap_year, am_year)[year] ~ dmvnormNC( sigma_year, Rho_year ),
    c(ap_id, am_id, bp_rank_id, bm_rank_id, bp_age_id, bm_age_id)[id1] ~ dmvnormNC( sigma_id, Rho_id ),
    c(ap_dyad, am_dyad, bp_rankDiff_dyad, bm_rankDiff_dyad)[dyad] ~ dmvnormNC( sigma_dyad, Rho_dyad ),
    ap ~ dnorm(0,2),
    am ~ dnorm(0,2),
    c(bp_rank, bm_rank, bp_age, bm_age, bp_kin, bm_kin, bp_rankDiff, bm_rankDiff, bp_ageDiff, bm_ageDiff, bp_kinXageDiff, bm_kinXageDiff, bp_kinXrankDiff, bm_kinXrankDiff, bp_prior2, bm_prior2) ~ dnorm(0,2),
    c(sigma_year, sigma_id, sigma_dyad) ~ dexp(1), 
    c(Rho_year, Rho_id, Rho_dyad) ~ dlkjcorr(3),
    scale ~ dexp(1)
    
  ),
  
  data=data.prior2.dgi, cores=3 , chains=3 , warmup=1000, iter=4000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE
  
)


m.prior3.dgi <- map2stan(
  alist(
    
    DGI ~ dzagamma2( p, mu , scale ),
    
    logit(p) ~ ap + ap_id[id1] + ap_id[id2] + ap_dyad[dyad] + ap_year[year] +
      (bp_rank + bp_rank_id[id1])*rank_id1 + (bp_rank + bp_rank_id[id2])*rank_id2 +
      (bp_age + bp_age_id[id1])*age_id1 + (bp_age + bp_age_id[id2])*age_id2 +
      (bp_rankDiff + bp_rankDiff_dyad[dyad])*rank_diff +
      bp_ageDiff*age_diff +
      bp_kin*maternal_kin + 
      bp_kinXageDiff*maternal_kin*age_diff +
      bp_kinXrankDiff*maternal_kin*rank_diff +
      bp_prior3*dgi_prior3,
    
    log(mu) ~ am + am_id[id1] + am_id[id2] + am_dyad[dyad] + am_year[year] +
      (bm_rank + bp_rank_id[id1])*rank_id1 + (bm_rank + bm_rank_id[id2])*rank_id2 +
      (bm_age + bp_age_id[id1])*age_id1 + (bm_age + bm_age_id[id2])*age_id2 +
      (bm_rankDiff + bm_rankDiff_dyad[dyad])*rank_diff +
      bm_ageDiff*age_diff +
      bm_kin*maternal_kin + 
      bm_kinXageDiff*maternal_kin*age_diff +
      bm_kinXrankDiff*maternal_kin*rank_diff +
      bm_prior3*dgi_prior3,
    
    c(ap_year, am_year)[year] ~ dmvnormNC( sigma_year, Rho_year ),
    c(ap_id, am_id, bp_rank_id, bm_rank_id, bp_age_id, bm_age_id)[id1] ~ dmvnormNC( sigma_id, Rho_id ),
    c(ap_dyad, am_dyad, bp_rankDiff_dyad, bm_rankDiff_dyad)[dyad] ~ dmvnormNC( sigma_dyad, Rho_dyad ),
    ap ~ dnorm(0,2),
    am ~ dnorm(0,2),
    c(bp_rank, bm_rank, bp_age, bm_age, bp_kin, bm_kin, bp_rankDiff, bm_rankDiff, bp_ageDiff, bm_ageDiff, bp_kinXageDiff, bm_kinXageDiff, bp_kinXrankDiff, bm_kinXrankDiff, bp_prior3, bm_prior3) ~ dnorm(0,2),
    c(sigma_year, sigma_id, sigma_dyad) ~ dexp(1), 
    c(Rho_year, Rho_id, Rho_dyad) ~ dlkjcorr(3),
    scale ~ dexp(1)
    
  ),
  
  data=data.prior3.dgi, cores=3 , chains=3 , warmup=1000, iter=4000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE
  
)


m.prior4.dgi <- map2stan(
  alist(
    
    DGI ~ dzagamma2( p, mu , scale ),
    
    logit(p) ~ ap + ap_id[id1] + ap_id[id2] + ap_dyad[dyad] + ap_year[year] +
      (bp_rank + bp_rank_id[id1])*rank_id1 + (bp_rank + bp_rank_id[id2])*rank_id2 +
      (bp_age + bp_age_id[id1])*age_id1 + (bp_age + bp_age_id[id2])*age_id2 +
      (bp_rankDiff + bp_rankDiff_dyad[dyad])*rank_diff +
      bp_ageDiff*age_diff +
      bp_kin*maternal_kin + 
      bp_kinXageDiff*maternal_kin*age_diff +
      bp_kinXrankDiff*maternal_kin*rank_diff +
      bp_prior4*dgi_prior4,
    
    log(mu) ~ am + am_id[id1] + am_id[id2] + am_dyad[dyad] + am_year[year] +
      (bm_rank + bp_rank_id[id1])*rank_id1 + (bm_rank + bm_rank_id[id2])*rank_id2 +
      (bm_age + bp_age_id[id1])*age_id1 + (bm_age + bm_age_id[id2])*age_id2 +
      (bm_rankDiff + bm_rankDiff_dyad[dyad])*rank_diff +
      bm_ageDiff*age_diff +
      bm_kin*maternal_kin + 
      bm_kinXageDiff*maternal_kin*age_diff +
      bm_kinXrankDiff*maternal_kin*rank_diff +
      bm_prior4*dgi_prior4,
    
    c(ap_year, am_year)[year] ~ dmvnormNC( sigma_year, Rho_year ),
    c(ap_id, am_id, bp_rank_id, bm_rank_id, bp_age_id, bm_age_id)[id1] ~ dmvnormNC( sigma_id, Rho_id ),
    c(ap_dyad, am_dyad, bp_rankDiff_dyad, bm_rankDiff_dyad)[dyad] ~ dmvnormNC( sigma_dyad, Rho_dyad ),
    ap ~ dnorm(0,2),
    am ~ dnorm(0,2),
    c(bp_rank, bm_rank, bp_age, bm_age, bp_kin, bm_kin, bp_rankDiff, bm_rankDiff, bp_ageDiff, bm_ageDiff, bp_kinXageDiff, bm_kinXageDiff, bp_kinXrankDiff, bm_kinXrankDiff, bp_prior4, bm_prior4) ~ dnorm(0,2),
    c(sigma_year, sigma_id, sigma_dyad) ~ dexp(1), 
    c(Rho_year, Rho_id, Rho_dyad) ~ dlkjcorr(3),
    scale ~ dexp(1)
    
  ),
  
  data=data.prior4.dgi, cores=3 , chains=3 , warmup=1000, iter=4000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE
  
)


m.prior5.dgi <- map2stan(
  alist(
    
    DGI ~ dzagamma2( p, mu , scale ),
    
    logit(p) ~ ap + ap_id[id1] + ap_id[id2] + ap_dyad[dyad] + ap_year[year] +
      (bp_rank + bp_rank_id[id1])*rank_id1 + (bp_rank + bp_rank_id[id2])*rank_id2 +
      (bp_age + bp_age_id[id1])*age_id1 + (bp_age + bp_age_id[id2])*age_id2 +
      (bp_rankDiff + bp_rankDiff_dyad[dyad])*rank_diff +
      bp_ageDiff*age_diff +
      bp_kin*maternal_kin + 
      bp_kinXageDiff*maternal_kin*age_diff +
      bp_kinXrankDiff*maternal_kin*rank_diff +
      bp_prior5*dgi_prior5,
    
    log(mu) ~ am + am_id[id1] + am_id[id2] + am_dyad[dyad] + am_year[year] +
      (bm_rank + bp_rank_id[id1])*rank_id1 + (bm_rank + bm_rank_id[id2])*rank_id2 +
      (bm_age + bp_age_id[id1])*age_id1 + (bm_age + bm_age_id[id2])*age_id2 +
      (bm_rankDiff + bm_rankDiff_dyad[dyad])*rank_diff +
      bm_ageDiff*age_diff +
      bm_kin*maternal_kin + 
      bm_kinXageDiff*maternal_kin*age_diff +
      bm_kinXrankDiff*maternal_kin*rank_diff +
      bm_prior5*dgi_prior5,
    
    c(ap_year, am_year)[year] ~ dmvnormNC( sigma_year, Rho_year ),
    c(ap_id, am_id, bp_rank_id, bm_rank_id, bp_age_id, bm_age_id)[id1] ~ dmvnormNC( sigma_id, Rho_id ),
    c(ap_dyad, am_dyad, bp_rankDiff_dyad, bm_rankDiff_dyad)[dyad] ~ dmvnormNC( sigma_dyad, Rho_dyad ),
    ap ~ dnorm(0,2),
    am ~ dnorm(0,2),
    c(bp_rank, bm_rank, bp_age, bm_age, bp_kin, bm_kin, bp_rankDiff, bm_rankDiff, bp_ageDiff, bm_ageDiff, bp_kinXageDiff, bm_kinXageDiff, bp_kinXrankDiff, bm_kinXrankDiff, bp_prior5, bm_prior5) ~ dnorm(0,2),
    c(sigma_year, sigma_id, sigma_dyad) ~ dexp(1), 
    c(Rho_year, Rho_id, Rho_dyad) ~ dlkjcorr(3),
    scale ~ dexp(1)
    
  ),
  
  data=data.prior5.dgi, cores=3 , chains=3 , warmup=1000, iter=4000, constraints=list(scale="lower=0") , control=list(adapt_delta=0.99), WAIC=TRUE
  
)



# END