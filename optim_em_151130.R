# updated: 11/30/2015
# add fs() functions to use Fisher scoring to update parameters
# updated: 9/23/2015
# replaced tr() in two positions. corrected score_bint_mat() at dat[,i]. EM for either positive or negative bint to get a better starting value.
# updated: 8/20/2015
# use ginv() instead of solve() to find the starting values
# MAJOR updated: 11/7/2014
# use EM algorithm instead of nloptr to find good start
# updated: 11/6/2014
# optimized score functions
# updated: 10/13/2014
# add the function of only specifying the intercept
# updated: 04/24/2014
# use nloptr to find the starting point when there is GxE interaction. Also modified grid_loglike accordingly.
# loglike multiplied by -1 since nloptr only minimizes.
# updated: 03/13/2014
# add grid function to get a good starting value for parameters
# also return hessian so can estimate se of the parameters



tr <- function(a,b) crossprod(as.vector(a),as.vector(b))        #trace of product of 2 square matrices
score_beta <- function(X,invV_resid) as.vector(t(X)%*%invV_resid)

score_sig_or_bint <- function(invV,invV_resid,dV_dpar){
    #dV_dpar: d(V)/d(parameter), parameter could be either sigma or bint
    (-0.5*tr(invV,dV_dpar)+0.5*t(invV_resid)%*%dV_dpar%*%invV_resid)[1]
}

score_bint_mat <- function (dat,beta,famid,ch_var){        #d(V)/d(bint), by famid
    #return a list of d(V)/d(bint), one element per bint
    vect = as.vector(1+as.matrix(dat)%*%matrix(beta,ncol=1))
    score <- list()
    for(i in 1:length(beta)){
        score[[i]] = by(cbind(vect,dat[,i]),INDICES=famid,function(v) {a=tcrossprod(v[,1],v[,2])*ch_var;a+t(a)})
    }
    
    score
}

info_sig_or_bint <- function(invV,dV_dpar_1,dV_dpar_2){
    #calculate fisher information
    #dV_dpar_1,dV_dpar_2: two dV_dpar matrix for different sig bint
    0.5*tr(invV%*%dV_dpar_1,invV%*%dV_dpar_2)
}

info_beta <- function(X,invV) t(X)%*%invV%*%X

ch_gxe_mat <- function(dat,beta,famid){        #generate CHE matrix when ther is GxE
    #dat: data.frame that contain the SNP data
    #beta: vector of interaction effect size(s)
    
    vect = as.vector(1+as.matrix(dat)%*%matrix(beta,ncol=1))
    by(vect,INDICES=famid,function(v) tcrossprod(v))
}

loglike_perFam <- function(theta,y,X,n_beta,n_var,...){
    beta = theta[1:n_beta]
    sig = theta[(n_beta+1):(n_beta+n_var)] 
    
    var = list(...)
    
    V = array(0,dim=dim(var[[1]]))    #initialize omega matrix
    for (i in 1:length(sig)){
        V = V+sig[i]*var[[i]]
    }
    
    n = dim(X)[1]
    L = try(t(chol(V)))
    
    if(class(L)=="try-error"){
        return(NaN)
    }
    
    invV_right = solve(L)%*%(y-as.matrix(X)%*%beta)
    
    loglike = -t(invV_right)%*%invV_right/2 - sum(log(diag(L))) - n*log(2*pi)/2
    
    loglike
}

loglike <- function(theta,y,X,varlist,gxe_int=FALSE,intX,famid){
    #intX: data frame of SNPs that have GxE
    n_beta = dim(X[[1]])[2]    
    n_var = length(varlist)
    if(gxe_int){
        n_bint = dim(intX)[2]
        bint = theta[(n_beta+n_var+1):length(theta)]
        varlist[[n_var-1]] = ch_gxe_mat(intX,bint,famid)        #re-write the CHE matrix, by default the next to the last var is CHE
    }else{
        n_bint = 0
    }
    
    if(length(theta)-n_beta-n_bint!=n_var){
        stop("Problem with the number of variance matrices!\n")
    }
    
    ll_call = paste0(paste("mapply(loglike_perFam,y=y,X=X",paste("varlist[[",1:n_var,"]]",collapse=","),"MoreArgs=list(theta=theta,n_beta=n_beta,n_var=n_var)",sep=","),")")        
    #"paste varlist": list varlist element one by one, for the sake of mapply
    ll_call = parse(text=ll_call)
    loglike = sum(eval(ll_call))
    
    loglike
}

e_step_perFam <- function(theta,y,X,n_beta,n_var,...){
    beta = theta[1:n_beta]
    sig = theta[(n_beta+1):(n_beta+n_var)]
    
    X = as.matrix(X)
    par = list(...)
    var = par[1:n_var]
    n_bint = length(par)-n_var
    
    V = array(0,dim=dim(var[[1]]))    #initialize omega matrix
    for (i in 1:length(sig)){
        V = V+sig[i]*var[[i]]
    }
    
    invV = solve(V)
    xb = X%*%beta
    resid = y-xb
    invV_resid = invV%*%resid
    
    t = numeric(length(sig))
    for (i in 1:length(sig)){          
        t[i] = sig[i]^2*t(invV_resid)%*%var[[i]]%*%invV_resid + sig[i]*dim(var[[i]])[1] - sig[i]^2*tr(invV,var[[i]])
    }
    
    s<-xb+sig[n_var]*invV_resid
    
    if(n_bint>0){
        g = numeric(n_bint)    #grad of bint
        info = matrix(0,nrow=n_bint,ncol=n_bint)
        s_bint = par[(n_var+1):(n_var+n_bint)]
        for (i in 1:n_bint){
            g[i] = score_sig_or_bint(invV,invV_resid,s_bint[[i]])
            #browser()
            for (j in 1:i){
                info[i,j] = info_sig_or_bint(invV,s_bint[[i]],s_bint[[j]])
            }
        }
        info[upper.tri(info)] = info[lower.tri(info)]
        
        return(list(t=t,s=s,g=g,info=info))
    }else{
        return(list(t=t,s=s))
    }
}

em <- function(theta,y,X_spl,varlist,gxe_int=FALSE,intX,famid,inv_xtx_xt,max_iter=10,verbose=FALSE){
    #inv_xtx_xt: ginv(t(X)%*%X)%*%t(X)
    n_beta = dim(X_spl[[1]])[2]
    n_var = length(varlist)
    n_ind = length(unlist(y))
    
    if(gxe_int==FALSE){
        n_bint = 0
        em_call = paste0("mapply(e_step_perFam,y=y,X=X_spl,",paste("varlist[[",1:n_var,"]]",collapse=","),
                               ",MoreArgs=list(theta=theta,n_beta=n_beta,n_var=n_var),SIMPLIFY = F)")       
        #"paste varlist": list varlist element one by one, for the sake of mapply
    }else{
        n_bint = dim(intX)[2]
        bint = theta[(n_beta+n_var+1):length(theta)]
        em_call = paste0("mapply(e_step_perFam,y=y,X=X_spl,",paste("varlist[[",1:n_var,"]]",collapse=","),",",
                               paste("dV_dbint[[",1:n_bint,"]]",collapse=","),",MoreArgs=list(theta=theta,n_beta=n_beta,n_var=n_var),SIMPLIFY = F)")
    }
    em_call = parse(text=em_call)
    
    i=1
    repeat{
        if(i>max_iter){
            conv<-0
            break
        }
        
        if(gxe_int){
            varlist[[n_var-1]] = ch_gxe_mat(intX,bint,famid)        #re-write the CHE matrix, by default the next to the last var is CHE
            ch_var = theta[n_beta+n_var-1]
            dV_dbint = score_bint_mat(intX,bint,famid,ch_var)        #d(V)/d(bint)
        }
        
        em_res = eval(em_call)
        
        suff = rowSums(sapply(em_res,"[[",1))/n_ind
        
        s = as.vector(sapply(em_res,"[[",2))
        beta_new<-inv_xtx_xt%*%s
        
        if(gxe_int){
            g = rowSums(matrix(sapply(em_res,"[[",3),nrow=n_bint))
            info = Reduce("+",lapply(em_res,"[[",4))
            bint = bint + solve(info)%*%g
            theta = c(beta_new,suff,bint)
        }else{
            theta = c(beta_new,suff)
        }
        theta_ll = c(beta_new,suff)
        
        if(verbose){
            ll = loglike(theta_ll,y=y,X=X_spl,varlist=varlist)
            cat(i,"\t",theta,"\t",ll,"\n")
        }
        
        i = i+1
    }
    
    #browser()
    theta
}

fs_perFam <- function(theta,y,X,n_beta,n_var,...){
    beta = theta[1:n_beta]
    sig = theta[(n_beta+1):(n_beta+n_var)]
    
    X = as.matrix(X)
    par = list(...)
    n_bint = length(par)-n_var
    
    V = array(0,dim=dim(par[[1]]))    #initialize omega matrix
    for (i in 1:n_var){
        V = V+sig[i]*par[[i]]
    }
    
#     invL = try(solve(t(chol(V))))
#     
#     if(class(invL)=="try-error"){
#         return(NaN)
#     }
#     
#     invV = t(invL)%*%invL
    invV = solve(V)
    xb = X%*%beta
    resid = y-xb
    invV_resid = invV%*%resid
    
    #score
    g <- numeric(length(theta))
    g[1:n_beta] = score_beta(X,invV_resid)
    for (i in 1:(n_var+n_bint)){
        g[n_beta+i] = score_sig_or_bint(invV,invV_resid,par[[i]])
    }
    
    #information
    i_beta = info_beta(X,invV)
    
    i_sig_and_bint = matrix(0,nrow=n_var+n_bint,ncol=n_var+n_bint)
    for (i in 1:(n_var+n_bint)){
        for (j in 1:i){
            i_sig_and_bint[i,j] = info_sig_or_bint(invV,par[[i]],par[[j]])
        }
    }
    i_sig_and_bint[upper.tri(i_sig_and_bint)] = i_sig_and_bint[lower.tri(i_sig_and_bint)]
    info_mat = as.matrix(Matrix::bdiag(i_beta,i_sig_and_bint))
    

#     i_sig = matrix(0,nrow=n_var,ncol=n_var)
#     for (i in 1:n_var){
#         for (j in 1:i){
#             i_sig[i,j] = info_sig_or_bint(invV,par[[i]],par[[j]])
#         }
#     }
#     i_sig[upper.tri(i_sig)] = i_sig[lower.tri(i_sig)]
# 
#     if(n_bint > 0){
#         i_bint = matrix(0,nrow=n_bint,ncol=n_bint)
#         for (i in 1:n_bint){
#             for (j in 1:i){
#                 i_bint[i,j] = info_sig_or_bint(invV,par[[n_var+i]],par[[n_var+j]])
#             }
#         }
#         i_bint[upper.tri(i_bint)] = i_bint[lower.tri(i_bint)]
#         info_mat = as.matrix(Matrix::bdiag(i_beta,i_sig,i_bint))
#     }else{
#         info_mat = as.matrix(Matrix::bdiag(i_beta,i_sig))
#     }

    return(list( grad = g, info = info_mat ))
}

fs <- function(theta,y,X_spl,varlist,gxe_int=FALSE,intX,famid,max_iter=50,reltol=1e-5,verbose=F){
    n_beta = dim(X_spl[[1]])[2]
    n_var = length(varlist)
    n_ind = length(unlist(y))
    
    if(gxe_int==FALSE){
        n_bint = 0
        fs_call = paste0("mapply(fs_perFam,y=y,X=X_spl,",paste("varlist[[",1:n_var,"]]",collapse=","),
                               ",MoreArgs=list(theta=theta,n_beta=n_beta,n_var=n_var),SIMPLIFY = F)")        
        #"paste varlist": list varlist element one by one, for the sake of mapply
    }else{
        n_bint = dim(intX)[2]
        bint = theta[(n_beta+n_var+1):length(theta)]
        fs_call = paste0("mapply(fs_perFam,y=y,X=X_spl,",paste("varlist[[",1:n_var,"]]",collapse=","),",",
                         paste("dV_dbint[[",1:n_bint,"]]",collapse=","),",MoreArgs=list(theta=theta,n_beta=n_beta,n_var=n_var),SIMPLIFY = F)") 
    }
    fs_call = parse(text=fs_call)
    
    i=1
    repeat{
        if(i>max_iter){
            conv <- -1
            break
        }
        
        theta_old = theta
        
        if(gxe_int){
            varlist[[n_var-1]] = ch_gxe_mat(intX,bint,famid)        #re-write the CHE matrix, by default the next to the last var is CHE
            ch_var = theta[n_beta+n_var-1]
            dV_dbint = score_bint_mat(intX,bint,famid,ch_var)        #d(V)/d(bint)
        }
        
        fs_res = eval(fs_call)
        
        g = rowSums(matrix(sapply(fs_res,"[[",1),nrow=length(theta)))
        info = Reduce("+",lapply(fs_res,"[[",2))
        
        if(i == 1){
            eig = eigen(info)
            if(any(eig$values<0)){
                mod_info = T
                mu = sqrt(max(abs(eig$values))/min(abs(eig$values)))
            }else{
                mod_info = F
            }
        }
        
        if(mod_info){
            theta = theta + solve(info+mu*diag(length(theta)))%*%g
        }else{
            theta = theta + solve(info)%*%g
        }
  
        #browser()
        #parameter space boundary for sigmas
        idx = which(theta[(n_beta+1):(n_beta+n_var)] <= 0)
        if(length(idx)>0){
            print('ha')
            print(theta)
            tot_var = sum(theta[(n_beta+1):(n_beta+n_var)])
            theta[(n_beta+1):(n_beta+n_var)][idx] <- tot_var/100
        }
        
        if(verbose){
            ll = loglike(theta,y=y,X=X_spl,varlist=varlist,gxe_int=gxe_int,intX,famid)
            cat(i,"\t",theta,"\t",ll,"\n")
        }
        
        if(all( abs((theta-theta_old)/theta_old) < reltol )){
            conv <- 0
            break
        }
        
        i = i+1
    }
    
    ll = loglike(theta,y=y,X=X_spl,varlist=varlist,gxe_int=gxe_int,intX,famid)
    
    return(list(convergence = conv, theta = as.vector(theta), loglike = ll, grad = g, info = info, iter = i-1))
}


myOptim <- function(formula,gdata,covar,kin=TRUE,ident=TRUE,ch=FALSE,gxe_int=FALSE,intSNP,max_em_iter=10,verbose_em=F){
    #gdata must have columns names "famid","id","fa","mo","sex" for family ID, individual ID, father ID,
    #    mother ID, sex, respectively.
    #covariates such as age and sex can be specified as covar=c("age","sex"), using the same column name
    #    as in gdata.
    
    if(!all(c("famid","id","fa","mo","sex") %in% names(gdata))){
        stop('Please check if you have all required columns ("famid","id","fa","mo","sex").')
    }
    
    require("kinship2")
    
    vars = all.vars(formula)        #all variables in formula
    rhs = strsplit(deparse(formula)," ~ ")[[1]][2]    #right hand side of formula
    
    y = gdata[,vars[1]]        #response
    
    if(rhs==1){
        X = as.matrix(rep(1,dim(gdata)[1]),ncol=1)
    }else{
        if(!missing(covar)){
            if(!all(covar %in% names(gdata))){
                stop("Any covariates not in the dataset?")
            }
            X = as.matrix(cbind(1,gdata[,covar],gdata[,vars[-1]]))        #as.matrix: for following initial est of beta
        }else{
            X = as.matrix(cbind(1,gdata[,vars[-1]]))
        }
    }
    
    famid = gdata$famid
    n_bint = 0
    bint = NULL
    if(kin==T){
        kinM = by(gdata,INDICES=famid,function(dat) 2*kinship(id=dat[[2]],dadid=dat[[3]],momid=dat[[4]],sex=dat[[5]]))
    }
    if(ch==T){
        chM = by(gdata,INDICES=famid,function(dat) matrix(1,nrow=dim(dat)[1],ncol=dim(dat)[1]))    
        if(gxe_int){        #having interaction
            if(missing(intSNP)){
                stop('Please correctly specify the SNP that interacts with CHE!')
            }else if(! intSNP %in% names(gdata)){
                stop('Please correctly specify the SNP that interacts with CHE!')
            }else{
                intX = gdata[,intSNP,drop=F]        #data frame of GxE interaction SNPs
                n_bint = dim(intX)[2]
            }
        }
    }else{        #no CHE
        if(gxe_int){
            stop("GxE interaction cannot be specified with 'ch=FALSE'!")
        }
    }
    if(ident==T){
        identM = by(gdata,INDICES=famid,function(dat) diag(dim(dat)[1]))
    }
    
    y_spl = split(y,famid)
    X_spl = split(as.data.frame(X),famid)
    
    require('MASS')

    n_var = sum(kin,ident,ch)        #number of variance components
    inv_xtx_xt = tcrossprod(ginv(t(X)%*%X),X)
    beta = inv_xtx_xt%*%y        #initial value

    totvar = var(y)
    
    mats = c("kinM","chM","identM")
    mats = mats[c(kin,ch,ident)]
    mats = paste(mats,collapse=',')
    
    if(gxe_int){
        #inital guess for bint either positive or negative
        time1 <- system.time( init1 <- em(c(beta,rep(totvar/n_var,n_var),1),y=y_spl,X_spl=X_spl,varlist=eval(parse(text=paste0('list(',mats,')'))),
                                         gxe_int = gxe_int,intX = intX,famid=famid,inv_xtx_xt=inv_xtx_xt,max_iter=max_em_iter,verbose = verbose_em) )
        time2 <- system.time( init2 <- em(c(beta,rep(totvar/n_var,n_var),-1),y=y_spl,X_spl=X_spl,varlist=eval(parse(text=paste0('list(',mats,')'))),
                                         gxe_int = gxe_int,intX = intX,famid=famid,inv_xtx_xt=inv_xtx_xt,max_iter=max_em_iter,verbose = verbose_em) )
        ll1 = loglike(init1,y=y_spl,X=X_spl,varlist=list(kinM,chM,identM),gxe_int=gxe_int,intX=intX,famid=famid)    #by default there will be all three VCs
        ll2 = loglike(init2,y=y_spl,X=X_spl,varlist=list(kinM,chM,identM),gxe_int=gxe_int,intX=intX,famid=famid)
        if(ll1>ll2) init = init1 else init = init2
    }else{
        time1 <- system.time( init <- em(c(beta,rep(totvar/n_var,n_var),bint),y=y_spl,X_spl=X_spl,varlist=eval(parse(text=paste0('list(',mats,')'))),
                                         gxe_int = gxe_int,inv_xtx_xt=inv_xtx_xt,max_iter=max_em_iter,verbose = verbose_em) )
    }
       

    if(gxe_int){
        opt_cmd = paste0('fs(theta=init,y=y_spl,X=X_spl,varlist=list(',mats,'),gxe_int=gxe_int,intX=intX,famid=famid)')
    }else{
        opt_cmd = paste0('fs(theta=init,y=y_spl,X=X_spl,varlist=list(',mats,'))')
    }

    time_fs <- system.time(result <- eval(parse(text=opt_cmd)))
    
    if(gxe_int){
        result$time_init_pos <- time1        #time for finding initial for positive bint
        result$time_init_neg <- time2        #time for finding initial for negative bint
    }else{
        result$time_init <- time1        #time for finding initial
    }
    
    result$time_opt <- time_fs        #time for optimization
    result$init <- init
    
    result
}
