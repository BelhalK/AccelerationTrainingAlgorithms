model2 <- inlineModel("
                      [COVARIATE]
                      input={p_F}

                      DEFINITION:
                      gender = { type        = categorical, 
                                 categories  = {0,1},
                                 P(gender=0) = p_F }
                      
                      [INDIVIDUAL]
                      input={e0_pop,o_e0,emax_pop,o_emax,gender, beta_F,e50_pop,o_e50}
                      gender={type=categorical,categories={0,1}}

                      
                      DEFINITION:
                      e0  ={distribution=lognormal, 
                            prediction=e0_pop,
                            sd=o_e0}
                      emax   ={distribution=lognormal, prediction=emax_pop,   sd=o_emax}
                      e50  ={distribution=lognormal, 
                              reference=e50_pop,  
                              covariate = gender,
                              coefficient  = {beta_F,0},
                              sd=o_e50}

                       [LONGITUDINAL]
                      input = {e0, emax, e50,a}

                      
                      EQUATION:
                      Cc = e0+emax*t/(e50+t)
                      
                      DEFINITION:
                      y1 ={distribution=normal, prediction=Cc, sd=a}

                      ")


adm  <- list(amount=1, time=seq(0,50,by=50))


p <- c(p_F=0.3,  beta_F=0.6,e0_pop=e0_true, o_e0=o_e0_true,
       emax_pop=emax_true, o_emax=o_emax_true, 
       e50_pop=e50_true, o_e50=o_e50_true,  
       a=0.1)
y1 <- list(name='y1', time=seq(0,to=50,by=10))
ind <- list(name=c("gender"))

res2a2 <- simulx(model = model2,
                 treatment = adm,
                 parameter = p,
                 group = list(size=100, level="covariate"),
                 output = list(ind, y1))

  writeDatamlx(res2a2, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/incr_pd.csv")
  pd.data <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/ISAEM/saemixB/data/incr_pd.csv", header=T, sep=",")
  pd.data <- pd.data[pd.data[,4]!=1 ,c(1,2,3,5)]
  pd.data[,3] <- res2a2$y1[,3]
  colnames(pd.data) <- c("subject","dose","response","gender")
