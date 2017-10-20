library(mlxR)
library(gridExtra)
theme_set(theme_bw())

model_macros <- inlineModel("
[LONGITUDINAL] 
input = {Ntt, ka, kEhc, Fent, kb, Cl, V, tg}

PK:
compartment(cmt=1, amount=Ehc)
compartment(cmt=2, amount=Ak)
iv(cmt=1, p=-1/amtDose)
iv(cmt=1, Tlag=tg, p=1/amtDose)
oral(cmt=2, Mtt=Ntt/ka, Ktr=ka, ka)

EQUATION:
t0 = 0
Ehc_0 = 1
ddt_Ehc=0 
ddt_Ak = Ehc*kEhc*Agb - ka*Ak
ddt_Acc = ka*Ak - Fent*kb*Acc - (1-Fent)*Cl/V*Acc
ddt_Agb = Fent*kb*Acc - Ehc*kEhc*Agb
Cc = Acc/V
")

p <- c(Cl=20, V=200, Ntt=3, kEhc=0.85, Fent=0.5, kb=0.04, ka=2.5, tg=6)
adm  <- list(amount=1, time=seq(0,50,by=12))
out  <- list(name = c('Ehc', 'Cc'), time = seq(0, 100, by=0.1))

res <- simulx(model=model_macros, 
              output=out, 
              treatment=list(adm), 
              parameter=p)

pl1 <- ggplot(res$Ehc) + geom_line(aes(time,Ehc))
pl2 <- ggplot(res$Cc)  + geom_line(aes(x=time,y=Cc))
grid.arrange(pl1,pl2)
## adm.w <- list(
##   tfd    = list(widget="slider", value=0,  min=0, max=24, step=2),
##   nd     = list(widget="slider", value=4,  min=0, max=12, step=1),
##   ii     = list(widget="slider", value=12, min=3, max=24, step=1),
##   amount = list(widget="slider", value=1,  min=0, max=20, step=1)
## )
## 
## p1 <- c(Cl=20, V=200, Ntt=3, ka=2.5)
## p2 <- c(tg=6, kEhc=0.85, Fent=0.5, kb=0.04)
## Ehc  <- list(name = 'Ehc', time = seq(0, 100, by=0.1))
## Cc  <- list(name = 'Cc' , time = seq(0, 100, by=0.1))
## out <- list(Ehc, Cc)
## 
## shiny.app <- shinymlx(model     = model_macros,
##                       output    = list(Ehc, Cc),
##                       treatment = adm.w,
##                       parameter = list(p1,p2),
##                       style     = "navbar1")
## shiny::runApp(shiny.app)
enterohepatic_model <- inlineModel("
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

p <- c(Cl_pop   = 20,     omega_Cl = 0.2, 
       V_pop    = 200,     omega_V = 0.2, 
       kEhc_pop = 0.85, omega_kEhc = 0.2,
       Fent_pop = 0.5,  omega_Fent = 0.2, 
       ka_pop   = 2.5,    omega_ka = 0.2, 
       tg_pop   = 6,      omega_tg = 0.2, 
       Ntt      = 3, 
       a        = 0.3, 
       b        = 0.02)

adm  <- list(amount=1000, time=seq(0,50,by=12))
Cc  <- list(name = 'Cc', time = seq(0, 100, by=0.1))
y   <- list(name = 'y',  time = c(1,2,3,4,5,6,7,8,9,12,24,27,30,33,36,48,54,60,72))

res <- simulx(model     = enterohepatic_model, 
              output    = y, 
              treatment = list(adm), 
              group     = list(size=100, level='individual'),
              parameter = p)

ggplot() +  geom_point(data=res$y, aes(x=time, y=y, color=id)) + 
  geom_line(data=res$y, aes(time,y, color=id)) +  theme(legend.position="none")

writeDatamlx(res, result.file="entero_data.csv")
