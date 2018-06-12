#!/usr/bin/Rscript

#title : 'Split and Filter the Paralogous sites'
#author: " Script from Sylvain Glémin (sylvain.glemin@univ-montp2.fr) / C++ version included in the read2snp software, J. David. 1/03/2012, adapted to our analyses by S. Arnoux"
#date : "January 2018"

### DEFINITION OF ARGUMENTS 
library("optparse")
rm(list=ls())
option_list = list(
make_option(c("-d", "--directory"), type="character", default=NULL, 
              help="directory location", metavar="character"),
make_option(c("-i", "--input"), type="character", default=NULL, 
              help="input file type file.out", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$directory)){
  print_help(opt_parser)
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

Pathdir <- opt$directory
FILE <- opt$input

setwd(Pathdir) 

##############################################################################################################
##############################################################################################################
#########################################                            #########################################
#########################################    SPLIT .out to FILES     #########################################
#########################################                            #########################################
##############################################################################################################
##############################################################################################################

Namesss<-"NewVcf_CSF_SPLIT_"
Predata<-read.table(FILE,sep=",") # ATTENTION au séparateur: peut être modifié au besoin
rownames(Predata) <- paste(Predata$V1,Predata$V2,sep = "_")
head(Predata)
data <- Predata[1:nrow(Predata),3:ncol(Predata)]
head(data)
arrondi_down <- function (x) ifelse(x-trunc(x)<=0,trunc(x)+1,trunc(x)) 
number <- arrondi_down(nrow(data)/20000) #get the number to split

SPLIT<- rep(1:number, each = 20000) # here the nrow(data) = 1811466
data$sub=c(SPLIT,rep((number+1), nrow(data)-length(SPLIT)))
spt1 <- split(data, data$sub) 
lapply(names(spt1), function(x){write.table(spt1[[x]], file = paste(Namesss, x, ".subfile.out", sep = ""))})

##############################################################################################################
##############################################################################################################
#########################################                            #########################################
#########################################    FILTER splitted FILES   #########################################
#########################################                            #########################################
##############################################################################################################
##############################################################################################################
library(gtools)

Files2filter<-list.files(path = Pathdir, pattern = ".subfile.out$", full.names = TRUE )

files2filter <- mixedsort(Files2filter)

FiltERCS <- function(FILE){
  pdata<-read.table(FILE,sep=" ") # ATTENTION au séparateur: peut être modifié au besoin
  data<-pdata[,1:(ncol(pdata)-1)]
  #head(data)

  #OUT<-paste(repertoire,"/",FILE,".res",sep="")
  OUT<-paste(FILE,".res",sep="")

  err<-0.001 # error threshold
  fis<-0.5 # FIS threshold # 0.3
  PREC<-0.001 # optimisation threshold

  POLYPARA<-TRUE # Booleen to know if we want the polymorphism details of the paralogous sites.
  SEUIL<-0.01 # p.value threshold for the difference between the model with or without paralogs above which we do not test the paralog polymorphisms

  #Fréquences des génotypes
  Fhom <-function(p,f) {p*p*(1-f)+f*p}
  Fhet <-function(p,f) {2*p*(1-p)*(1-f)}

  #Proportion de reads attendues
  Phom1 <-function(e) {c(1-3*e,e,e,e)}
  Phet <-function(e) {c((1-2*e)/2,(1-2*e)/2,e,e)}
  Phom2 <-function(e) {c(e,1-3*e,e,e)}

  #Proba pour un génotype possible
  pgeno<- function(vect,proba,freq){
    freq*dmultinom(vect,prob=proba)
  }

  #Log Vraissemblance pour un individu: Model0: heterozygote
  lnL0 <-function(vect,e,f,p1){
    A<-pgeno(vect,Phom1(e),Fhom(p1,f))
    A<-A+pgeno(vect,Phet(e),Fhet(p1,f))
    A<-A+pgeno(vect,Phom2(e),Fhom(1-p1,f))
    log(A)
  }


  #Log Vraissemblance pour un individu: Model1: paralogue
  lnL1 <-function(vect,e,f,p1,p2,x){
    A<-pgeno(vect,x*Phom1(e)+(1-x)*Phom1(e),Fhom(p1,f)*Fhom(p2,f))
    A<-A+pgeno(vect,x*Phom1(e)+(1-x)*Phet(e),Fhom(p1,f)*Fhet(p2,f))
    A<-A+pgeno(vect,x*Phom1(e)+(1-x)*Phom2(e),Fhom(p1,f)*Fhom(1-p2,f))
    A<-A+pgeno(vect,x*Phet(e)+(1-x)*Phom1(e),Fhet(p1,f)*Fhom(p2,f))
    A<-A+pgeno(vect,x*Phet(e)+(1-x)*Phet(e),Fhet(p1,f)*Fhet(p2,f))
    A<-A+pgeno(vect,x*Phet(e)+(1-x)*Phom2(e),Fhet(p1,f)*Fhom(1-p2,f))
    A<-A+pgeno(vect,x*Phom2(e)+(1-x)*Phom1(e),Fhom(1-p1,f)*Fhom(p2,f))
    A<-A+pgeno(vect,x*Phom2(e)+(1-x)*Phet(e),Fhom(1-p1,f)*Fhet(p2,f))
    A<-A+pgeno(vect,x*Phom2(e)+(1-x)*Phom2(e),Fhom(1-p1,f)*Fhom(1-p2,f))
    log(A)
  }

  assignation <- function(geno,e,f) {
    #geno<-as.matrix(matrix(as.numeric(x),ncol=4,byrow=TRUE))
    Nind<-length(geno)/4
    #Sort the 4 bases (A,C,G,T) in decreasing order of the number of reads
    ordre<-order(c(sum(geno[,1]),sum(geno[,2]),sum(geno[,3]),sum(geno[,4])),decreasing=TRUE)
    tab<-geno[,ordre]
    
    #Total loglikelihood function, model0: heterozygotes
    logvrais0 <- function(p) {
      LL<-0
      for(i in 1:Nind){
        LL<-LL+lnL0(as.numeric(tab[i,]),e,f,p)
      }
      -LL
    }
    
    #Total loglikelihood function, model1: paralogues
    logvrais1 <- function(par) {
      LL<-0
      for(i in 1:Nind){
        LL<-LL+lnL1(as.numeric(tab[i,]),e,f,par[1],par[2],par[3])
      }
      -LL
    }
    
    #Total loglikelihood function, model2A: sub-model of the model 1 where p1 is fixed either to 1 or 0 (p2 and x are estimated)
    #To test wether p1=1
    logvrais2Asup <- function(par) {
      LL<-0
      for(i in 1:Nind){
        LL<-LL+lnL1(as.numeric(tab[i,]),e,f,1-PREC,par[1],par[2])
      }
      -LL
    }
    #To test wether p1=0
    logvrais2Ainf <- function(par) {
      LL<-0
      for(i in 1:Nind){
        LL<-LL+lnL1(as.numeric(tab[i,]),e,f,PREC,par[1],par[2])
      }
      -LL
    }
    
    #Total loglikelihood function, model2B: sub-model of the model 1 where p2 is fixed either to 1 or 0 (p1 and x are estimated)
    #To test wether p2=1
    logvrais2Bsup <- function(par) {
      LL<-0
      for(i in 1:Nind){
        LL<-LL+lnL1(as.numeric(tab[i,]),e,f,par[1],1-PREC,par[2])
      }
      -LL
    }	
    #To test wether p2=1
    logvrais2Binf <- function(par) {
      LL<-0
      for(i in 1:Nind){
        LL<-LL+lnL1(as.numeric(tab[i,]),e,f,par[1],PREC,par[2])
      }
      -LL
    }
    
    #starting values for optimization
    n1<-sum(tab[,1])
    n2<-sum(tab[,2])
    xi<-n1/(n1+n2)
    
    #Optimization of model 0
    if (is.finite(logvrais0(xi))) {
      ml0<- tryCatch(optim(par = xi,fn = logvrais0,lower=PREC,upper=1-PREC,method="L-BFGS-B"),error=function(e) NA)
      if (any(is.na(ml0))) {
        return(rep("NA",6))
      } else {
        #Optimization of model1
        #Boundaries for model 1
        min<-c(PREC,PREC,0.5) #By convention optimization is constrained to x > 1/2, hence the focal gene is assumed to be the most expressed one
        max<-c(1-PREC,1-PREC,1-PREC) 
        #3 starting points can be used used	
        init<-c(1-PREC,PREC,xi)
        ml1a<-tryCatch(optim(par=init,fn=logvrais1,lower=min,upper=max,method="L-BFGS-B"),error=function(e) NA)
        if (any(is.na(ml1a))) {
          return(rep("NA",6))
        } else {
          #Check of the convergence: the likelihood must be lower or equal than logvrais(p_estim0,p_estim0,1-PREC)
          #where p_estim0 is the frequency estimated in model 0
          #This comparison is done instead of the direct comparison with the likelihood of model 0
          #because of the boundaries that are not 0 and 1
          ml1<-ml1a
          maxLL<-logvrais1(c(ml0$par,ml0$par,1-PREC))
          if(ml1a$value>maxLL){
            init<-c(PREC,1-PREC,xi)
            ml1b<-optim(par=init,fn=logvrais1,lower=min,upper=max,method="L-BFGS-B")
            ml1<-ml1b
            if(ml1b$value>maxLL){
              init<-c(xi,xi,0.5)
              ml1c<-optim(par=init,fn=logvrais1,lower=min,upper=max,method="L-BFGS-B")
              #In case of model 0 has the best likelihood because of approximation problem the best ml1 is chosen
              ml1<-ml1c
              if(ml1b$value>maxLL){
                if(ml1b$value<ml1a$value & ml1b$value<=ml1c$value) ml1<-ml1b
                if(ml1c$value<ml1a$value & ml1c$value<ml1b$value) ml1<-ml1c
              }
            }	
          }
          
          
          #Computation of the p.values: Model1 vs model0
          p.value_para<-pchisq(max(0,2*(ml0$value-ml1$value)),2,lower.tail=FALSE)
          
          if(POLYPARA==FALSE) return(p.value_para)
          
          #Optional procedures
          
          if(p.value_para > SEUIL) return(c(ml1$par,p.value_para,"NA","NA"))
          
          #Optimization of model2A and model 2B
          #Boundaries for model 2A
          min<-c(PREC,0.5)
          max<-c(1-PREC,1-PREC) 
          if(ml1$par[1]>=0.5){
            init<-c(PREC,xi)
            ml2A<-tryCatch(optim(par=init,fn=logvrais2Asup,lower=min,upper=max,method="L-BFGS-B"),error=function(e) NA)
            init<-c(1-PREC,xi)
            ml2B<-tryCatch(optim(par=init,fn=logvrais2Binf,lower=min,upper=max,method="L-BFGS-B"),error=function(e) NA)
          }
          if(ml1$par[1]<0.5){
            init<-c(1-PREC,xi)
            ml2A<-tryCatch(optim(par=init,fn=logvrais2Ainf,lower=min,upper=max,method="L-BFGS-B"),error=function(e) NA)
            init<-c(PREC,xi)
            ml2B<-tryCatch(optim(par=init,fn=logvrais2Bsup,lower=min,upper=max,method="L-BFGS-B"),error=function(e) NA)
          }
          
          #Computation of the p.values
          #Model1 vs model2A
          if (any(is.na(ml2A))) {
            p.value_mono1 <- "NA"
          } else {
            p.value_mono1<-pchisq(max(0,2*(ml2A$value-ml1$value)),2,lower.tail=FALSE)
          }
          #Model1 vs model2B
          if (any(is.na(ml2B))) {
            p.value_mono2 <- "NA"
          } else {
            p.value_mono2<-pchisq(max(0,2*(ml2B$value-ml1$value)),2,lower.tail=FALSE)
          }
        
        #Return p1, p2, x and the 3 p.values
        return(c(ml1$par,p.value_para,p.value_mono1,p.value_mono2))
        }
      }
    } else {
      return(rep("NA",6))
    }  
  }

  #Boucle générale
  write(c("p1","p2","x","p.value_para","p.value_mono1","p.value_mono2"),OUT,ncolumn=6)

  generalAssignation <- function(x){
    genotype<-as.matrix(matrix(as.numeric(x),ncol=4,byrow=TRUE))
    resultat<-assignation(genotype,err,fis)
    write(resultat,OUT,append=TRUE,ncolumn=6)
  }

  res <- apply(data, 1 , generalAssignation)
}


library(parallel)
# Calculate the number of cores
no_cores <- detectCores() - 1
 
# Initiate cluster
cl <- makeCluster(no_cores)

#parSapply allow the sapply to be parallelized and to go faster :D
resfil <- parSapply(cl, files2filter, FiltERCS)
stopCluster(cl)
##############################################################################################################
##############################################################################################################
#########################################                            #########################################
#########################################    MERGING .RES FILES      #########################################
#########################################                            #########################################
##############################################################################################################
##############################################################################################################

cl<-makeCluster(no_cores)

RESOUT<-paste(FILE,".Merged.res",sep="")

#Get the files that have been resulting from firt steps of filtering file(x).out -> file(x).out.res

Files2Merge<-list.files(path = Pathdir, pattern = ".out.res$", full.names = TRUE )

files <- mixedsort(Files2Merge)

clusterExport(cl, "RESOUT")

mergeAllRes <- function(x){
  NAMES<-paste(sub(".out.res",".out",x))
  RES<-paste(x)
  a_joindre <- read.table(RES,sep=" ",header=T) 
  NameS <- read.table(NAMES,sep=" ",header=T)
  Name_Contigs <- unlist(strsplit(rownames(NameS),"_"))[ c(TRUE,FALSE) ]
  Name_Sites <- unlist(strsplit(rownames(NameS),"_"))[ c(FALSE,TRUE) ]
  c_joindre <- data.frame(Contig = Name_Contigs[1:(nrow(a_joindre))], Site = as.numeric(Name_Sites[1:(nrow(a_joindre))]), a_joindre) 
  write.table(c_joindre, RESOUT, append = TRUE, col.names = FALSE, row.names = FALSE)
}

resmer <- parSapply(cl, files, mergeAllRes)

stopCluster(cl)
