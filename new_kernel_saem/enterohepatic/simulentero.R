library(mlxR)


model2 <- inlineModel("
    [LONGITUDINAL] 
    input = {Ntt, ka, kEhc, Fent, Cl, V, tg, a, b}

    PK:
    compartment(cmt=1, amount=Ehc)
    compartment(cmt=2, amount=Ak)
    iv(cmt=1, p=-1/amtDose)
    iv(cmt=1, Tlag=tg, p=1/amtDose)
    oral(cmt=2, Mtt=Ntt/ka, Ktr=ka, ka)

    EQUATION:
    kb = Cl/V
    t0 = 0
    Ehc_0 = 1
    ddt_Ehc=0 
    ddt_Ak = Ehc*kEhc*Agb - ka*Ak
    ddt_Acc = ka*Ak - kb*Acc 
    ddt_Agb = Fent*kb*Acc - Ehc*kEhc*Agb
    Cc = Acc/V
    g = a + b*Cc

    DEFINITION:
    y = {distribution=normal, prediction=Cc, sd=g}

    [INDIVIDUAL]
    input = {Cl_pop, omega_Cl, V_pop, omega_V, kEhc_pop, omega_kEhc, 
    Fent_pop, omega_Fent, ka_pop, omega_ka, tg_pop, omega_tg}

    DEFINITION:
    tg = {distribution=lognormal, prediction=tg_pop, sd=omega_tg}
    Cl = {distribution=lognormal, prediction=Cl_pop, sd=omega_Cl}
    V = {distribution=lognormal, prediction=V_pop, sd=omega_V}
    ka = {distribution=lognormal, prediction=ka_pop, sd=omega_ka}
    kEhc = {distribution=lognormal, prediction=kEhc_pop, sd=omega_kEhc}
    Fent = {distribution=logitnormal, prediction=Fent_pop, sd=omega_Fent}
                      ")

adm  <- list(amount=1, time=seq(0,50,by=50))
p <- c(base_pop=1, o_base=1,
       slope_pop=1.5, o_slope=1.5,  
       a=1)
y1 <- list(name='y1', time=seq(1,to=50,by=5))


res2a2 <- simulx(model = model2,
                 treatment = adm,
                 parameter = p,
                 group = list(size=10, level="individual"),
                 output = y1)


writeDatamlx(res2a2, result.file = "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/enterohepatic/entero_synth.csv")
table <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/enterohepatic/entero_synth.csv", header=T, sep=",")
head(table)
table[1:45,]


#modification for mlxsaem dataread function
obj <- read.table("/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/enterohepatic/entero_synth.csv", header=T, sep=",")
obj <- obj[obj$amount !=1,]
obj[,4] <- 1
write.table(obj, "/Users/karimimohammedbelhal/Documents/GitHub/saem/new_kernel_saem/enterohepatic/entero_synth.csv", sep=",", row.names=FALSE,quote = FALSE, col.names=TRUE)


