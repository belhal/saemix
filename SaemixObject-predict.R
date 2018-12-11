
setMethod(f="predict",
          signature="SaemixObject",
          def=function(object,newdata=NULL,type=c("ipred", "ypred", "ppred", "icpred"), se.fit=FALSE, ...) {
            namObj<-deparse(substitute(object))
            type<-match.arg(type)
            print(namObj)
            saemix.data<-object["data"]
            saemix.model<-object["model"]
            #    se.fit<-match.arg(se.fit) # doesn't work with logical type, change
            #    if(se.fit) cat("Currently predict() does not handle argument se.fit=TRUE.\n")
            if(missing(newdata)) { # Return predictions from fitted object
              xpred<-fitted(object,type,...)
              if(length(xpred)==0) {
                opred<-saemix.predict(object)
                assign(namObj,opred,envir=parent.frame()) # update object invisibly with predictions
                print(parent.frame())
                print(namObj)
                xpred<-fitted(object,type,...)
              }
            } else {# Ignore type - when newdata is given, predictions correspond to the population predictions using the final estimates
              Mcov<-object["model"]["Mcovariates"]
              valnames<-unique(c(names(Mcov)[names(Mcov)!="id"],object["data"]["name.group"],object["data"]["name.predictors"]))
              ival<-match(valnames,names(newdata))
              if(sum(is.na(ival))>0) {
                cat("The new dataframe must contain the predictors, patient id and covariates in the model, with the same names as the original dataset.")
                cat("The following items were not found:",valnames[is.na(ival)],"\n")
                return()
              }
              if(type!="ppred") {
                ival2<-match(object["data"]["name.response"],names(newdata))
                cat("Currently only the population predictions are provided for a new dataset, option type is ignored.\n")
#                if(is.na(ival2)) cat("For individual predictions, a column must be provided with the responses measured for the new individuals.\n")
              }
              ndat<-dim(newdata)[1]
              ndat.ind<-c(tapply(unlist(newdata[,object["data"]["name.group"]]),unlist(newdata[,object["data"]["name.group"]]),length))
              nsuj<-length(unique(newdata[,object["data"]["name.group"]]))
              for(varname in c(object["data"]["name.mdv"],object["data"]["name.cens"])) {
                if(is.na(match(varname,names(newdata)))) {
                  newdata<-cbind(newdata,rep(0,ndat))
                  colnames(newdata)[dim(newdata)[2]]<-varname
                }
              }
              for(varname in c(object["data"]["name.occ"],object["data"]["name.ytype"])) {
                if(is.na(match(varname,names(newdata)))) {
                  newdata<-cbind(newdata,rep(1,ndat))
                  colnames(newdata)[dim(newdata)[2]]<-varname
                }
              }
              if(is.na(match(object["data"]["name.response"],names(newdata)))) {
                newdata<-cbind(newdata,rep(NA,ndat))
                colnames(newdata)[dim(newdata)[2]]<-object["data"]["name.response"]
              }
              newdata<-cbind(newdata,index=rep(1:nsuj,each=ndat.ind))
              object["data"]["data"]<-newdata
              object["data"]["N"]<-nsuj
              object["data"]["ntot.obs"]<-ndat
              object["data"]["nind.obs"]<-ndat.ind
              Mcovariates<-data.frame(id=rep(1,object["data"]["N"]))
              if(dim(Mcov)[2]>1) { # Recover covariates
                name.cov<-names(Mcov)[names(Mcov)!="id"]
                idx<-which(!duplicated(newdata[,object["data"]["name.group"]]))
                Mcovariates<-cbind(Mcovariates,newdata[idx,name.cov,drop=FALSE])
              }
              for(icol in dim(Mcovariates)[2])
                if(is.factor(Mcovariates[,icol])) Mcovariates[,icol]<-as.numeric(Mcovariates[,icol])-1
              # COV: design matrix
              COV<-matrix(nrow=dim(Mcovariates)[1],ncol=0)
              for(j in 1:object["model"]["nb.parameters"]) {
                jcov<-which(object["model"]["betaest.model"][,j]==1)
                aj<-as.matrix(Mcovariates[,jcov])
                COV<-cbind(COV,aj)
              }
              # population parameters given by phi=Ci*mu (COV %*% Lcovariates) and psi=h(phi)
              pop.phi<-COV %*% object["results"]["MCOV"]
              pop.psi<-transphi(pop.phi,object["model"]["transform.par"])
              structural.model<-object["model"]["model"]
              chdat<-new(Class="SaemixRepData",data=object["data"], nb.chains=1)
              IdM<-chdat["dataM"]$IdM
              XM<-chdat["dataM"][,c(saemix.data["name.predictors"],saemix.data["name.cens"],saemix.data["name.mdv"],saemix.data["name.ytype"]),drop=FALSE]
              # Model predictions with pop.psi
              fpred<-structural.model(pop.psi, IdM, XM)
              ind.exp<-which(saemix.model["error.model"]=="exponential")
              for(ityp in ind.exp) fpred[XM$ytype==ityp]<-log(cutoff(fpred[XM$ytype==ityp]))
              xpred<-fpred
            }
            return(xpred)
          }
)
