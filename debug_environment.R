######################################################################################################
#### Fonction explore.envir pour accéder au nom d'un objet à l'intérieur d'une méthode S4

setGeneric(name="explore.envir",
           def=function(object,type=c("ppred"),...){standardGeneric("explore.envir")}
)

setMethod(f="explore.envir",
          signature="numeric",
          def=function(object,type=c("ppred"),...) {
            namObj<-deparse(substitute(object))
            print(namObj)
            print(object)
            assign(namObj,4,envir=parent.frame()) # update object invisibly with predictions
          })

setMethod(f="explore.envir",
          signature="SaemixObject",
          def=function(object,type=c("ppred"),...) {
            namObj<-deparse(substitute(object))
            type<-match.arg(type)
            print(namObj)
#            print(object)
          })

cat("Comportement attendu, affiche 6\n")
explore.envir(6)

cat("Comportement attendu, affiche 'x' et modifie la valeur de x à 4\n")
x<-6
explore.envir(x)
cat("Dans une fonction cat:    ",explore.envir(x),"=",x,"\n")

#### Fonction predict dans saemix
progDir<-"/Users/karimimohammedbelhal/Desktop/Phd/R_Package/contributions/git-saemixextension"
datDir<-"/Users/karimimohammedbelhal/Desktop/Phd/R_Package/contributions/git-saemixextension/data"

cat("Library files\n")
{
  source(file.path(progDir,"R","aaa_generics.R"))
  source(file.path(progDir,"R","SaemixData.R"))
  source(file.path(progDir,"R","SaemixRes.R"))
  source(file.path(progDir,"R","SaemixModel.R"))
  source(file.path(progDir,"R","SaemixObject.R"))
  source(file.path(progDir,"R","main.R"))
  source(file.path(progDir,"R","func_aux.R"))
  source(file.path(progDir,"R","main_initialiseMainAlgo.R"))
  source(file.path(progDir,"R","main_estep.R"))
  source(file.path(progDir,"R","main_mstep.R"))
  source(file.path(progDir,"R","func_FIM.R"))
  source(file.path(progDir,"R","func_plots.R"))
  source(file.path(progDir,"R","func_distcond.R"))
  source(file.path(progDir,"R","func_simulations.R"))
  source(file.path(progDir,"R","compute_LL.R"))
  }

theo.saemix<-read.table(file.path(datDir,"theo.saemix.tab"),header=T)
saemix.data<-saemixData(name.data=theo.saemix,header=TRUE,sep=" ",na=NA, 
                        name.group=c("Id"),name.predictors=c("Dose","Time"),
                        name.response=c("Concentration"),name.covariates=c("Weight","Sex"),
                        units=list(x="hr",y="mg/L",covariates=c("kg","-")), name.X="Time")

model1cpt<-function(psi,id,xidep) { 
  dose<-xidep[,1]
  tim<-xidep[,2]  
  ka<-psi[id,1]
  V<-psi[id,2]
  CL<-psi[id,3]
  k<-CL/V
  ypred<-dose*ka/(V*(ka-k))*(exp(-k*tim)-exp(-ka*tim))
  return(ypred)
}
# Default model, no covariate
saemix.model<-saemixModel(model=model1cpt,type="structural",
                          description="One-compartment model with first-order absorption",
                          psi0=matrix(c(1.,20,0.5,0.1,0,-0.01),ncol=3,byrow=TRUE, 
                                      dimnames=list(NULL, c("ka","V","CL"))),transform.par=c(1,1,1))

# Note: remove the options save=FALSE and save.graphs=FALSE to save the results and graphs
saemix.options<-list(seed=632545,save=FALSE,save.graphs=FALSE)
saemix.fit<-saemix(saemix.model,saemix.data,saemix.options)

source('/Users/karimimohammedbelhal/Desktop/Phd/R_Package/contributions/git-saemixextension/SaemixObject-predict.R')
explore.envir(saemix.fit)
cat("Eco: affiche saemix.fit comme prevu\n")

mypred<-predict(saemix.fit)
cat("Eco: affiche object, et saemix.fit n'est pass updaté donc mypred est vide \n")
head(mypred)

##############################s########################################################################
# Un autre truc qui ne va pas
print.foo1=function(x){ 
  print(deparse(substitute(x))) 
}
test1 <- list(a=1, b=2)
test1
#this shows "test" as expected
print(test1)
print.foo1(test1)

#this should show
#"structure(list(a = 1, b = 2), .Names = c(\"a\", \"b\"), class = \"foo\")"
class(test1)="foo1"
test1
# but it shows "x" 
# why ???

######################################################################################################
# Understanding assignements

x <- 0
y <- 10
f <- function() {
  x <- 1
  g()
}
g <- function() {
  x <- 2
  h()
}
h <- function() {
  x <- 3
  x + y
}
f()
h()

######################################################################################################
