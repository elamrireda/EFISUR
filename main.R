##############  TEST

Objective_function <- fct1
Constraint_function <- list(fct2)

d <- 1 # dim of the control space
m <- 1 # dim of the uncertain space
bornesControl <- rbind(c(-1,1)) #X-space
bornesUncert <- rbind(c(0.35,0.65)) # U-space
bornes <- rbind(bornesControl,bornesUncert) #joint-space
seuil <- -2 #c(-57,-32,0.59,0.57,1.42) # threshold
sign <- "<" #c(">",">",">","<","<")
alpha <- 0.05 # parameter that controls the probaility of constraint satisfaction
N <- 200 # realizations to estimate the probability of A
M_CRN <- 50 # CRN (commmun random numbers)
alea <- matrix(randtoolbox::sobol(M_CRN,m),ncol=m)
for(i in 1:m) alea[,i] <- bornes[(i+d),1] + (bornes[(i+d),2] - bornes[(i+d),1])*alea[,i]

############################
##############  INITIAL DoE
############################

n = 100*(d+m) # initial design size
DoE = lhsDesign(n,d+m)$design # initial LHS design in the joint space
DoE  <- maximinSA_LHS(DoE)$design # maximin-LHS design in the joint space
for(i in 1:(d+m)) DoE[,i] <- bornes[i,1] + (bornes[i,2] - bornes[i,1])*DoE[,i]

############################
##############  INITIAL GP models
############################
 
d_c <- length(Constraint_function) # constraints models
listCon <- NULL
response_f <- apply(DoE, 1, Objective_function) # objective function responses f(x,u)
model_O <- km(formula= ~1, design = data.frame(x=DoE), response= response_f,estim.method="LOO",
              control = list(trace=F),lower=rep(1e-10,d+m)) # Objective function metamodel (F)
for(l in 1:d_c){
  response_const <- apply(DoE, 1, Constraint_function[[l]]) # constraint function responses g(x,u)
  model_Const <- km(formula= ~1, design = data.frame(x=DoE), response= response_const,estim.method="LOO",
                    control = list(trace=F),lower=rep(1e-10,d+m)) # Constraint function metamodel (G)
  listCon[[l]] <- model_Const
}

############################
##############  EFISUR
############################
solver <- "bobyqa"
res <- EFISUR(model_O,listCon,bornes,d,m,alea,N,alpha,sign,seuil,solver)
