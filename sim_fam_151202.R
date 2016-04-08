# created: 01/30/2014
# include a single function to simulate family data, which can accept parameters
# simulate G E effects in N extended families (n_snp G 1E + mGxE multiple interactions)
# modified for sim_multi_GxE_extFam.R
# updated: 10/13/2014
# add back sig_a, whenever the exact additive var is fixed and no need to simulate each SNP
# updated: 11/19/2014
# SNPs and sig_a can be specified and added together
# updated: 07/21/2015
# add another argument to control whether simulate the case where there is regular environment factor that has interaction with G, rather than CHE
# this environment is sampled iid from some Normal distribution
# updated: 12/02/2015
# add a second CHE to see the robustness of the estimation model


sim_fam <- function(n_fam=100,n_chd=1,n_snp=10,n_isnp=1,sig_a=NULL,sig_E=20,sig_e=10,p_a=FALSE,b_G=FALSE,
                    rEnv = FALSE,b_GxE=FALSE,p_a_rg=c(0.1,0.3),b_G_rg=c(0.5,2),b_GxE_rg=c(0.1,0.5),sig_E2=NULL,seed=FALSE){
    #n_fam: number of nuclear families
    #n_chd: number of children per family
    #n_snp: number of SNPs
    #n_isnp: number of SNPs having interaction with E
    #sig_a: sigma of additive genetic effect, overides SNPs if specified
    #sig_E: sigma of environment
    #sig_e: sigma of residual
    #p_a: MAF explicitly specified if not FALSE; should have same length with n_snp
    #b_G: beta for G explicitly specified if not FALSE; should have same length with n_snp
    #rEnv: if the environment factor is a regular one or not (CHE)
    #b_GxE: beta for interaction explicitly specified if not FALSE; should have same length with n_isnp
    #p_a_rg: MAF range using a unif; or we could sample using a beta dist; can provide fixed number
    #b_G_rg: range of beta for G; Can try negative values, eg -2:2; can provide fixed number
    #b_GxE_rg: range of beta for interaction; can provide fixed number
    #seed: random number generator seed, must be a number. if not provided, a different dataset would be generated each time
    
    #== check ==#

    args = lapply(as.list(match.call())[-1],eval)        #eval to find arguments' values

    if(identical(p_a,F)){
        args$p_a=NULL
    }else if(length(args$p_a)!=n_snp){
        stop("The lengths of n_snp and p_a are not consistent!") 
    }
    
    if(identical(b_G,F)){
        args$b_G=NULL
    }else if(length(args$b_G)!=n_snp){
        stop("The lengths of n_snp and b_G are not consistent!") 
    }
    
    if(identical(b_GxE,F)){
        args$b_GxE=NULL
    }else if(length(args$b_GxE)!=n_isnp){
        stop("The lengths of n_isnp and b_GxE are not consistent!") 
    }
    
    if(!all(sapply(args[!names(args)=='rEnv'],is.numeric))){
        stop("Check the arguments!")
    }

    #== load package ==#

    require("MASS")
    require('kinship2')
    
    #== parameters ==#

    n_indPerFam = 2+n_chd        #number of individuals per family (+2 parents)


    p_trans = matrix(
        c(1,0.5,0,0.25,0,0,
            0,0.5,1,0.5,0.5,0,
            0,0,0,0.25,0.5,1),
        nrow=3,byrow=T)        #transmission prob
    #         aaxaa    aaxAa    aaxAA    AaxAa    AaxAA    AAxAA    (parents)
    # aa
    # Aa
    # AA
    # (child)
    
    p_trans_colidx = matrix(
        c(6,5,3,
        5,4,2,
        3,2,1),
        nrow=3)        #column index for p_trans according to parents geno
    #         AA    Aa    aa    (dad)
    # AA
    # Aa
    # aa
    # (mom)
    

    #== Simulation ==#
    #pedigree
    famid <- rep(1:n_fam,each=n_indPerFam)
    id <- paste0(famid,0,rep(1:n_indPerFam,times=n_fam))
    fa <- rep(NA,times=n_indPerFam*n_fam)
    chld_idx = rep(seq(2,by=n_indPerFam,length.out=n_fam),times=n_chd)+rep(seq(1,n_chd),each=n_fam)        #indeces of children
    fa[chld_idx] = id[seq(1,by=n_indPerFam,length.out=n_fam)]
    mo <- rep(NA,times=n_indPerFam*n_fam)
    mo[chld_idx] = id[seq(2,by=n_indPerFam,length.out=n_fam)]
    sex <- rep(c(1,2,rep(1,n_chd)),times=n_fam)

    gped <- kinship2:::pedigree(id, fa, mo, sex=sex, famid=famid)

    #G
    add_var = 0
    if(n_snp != 0){
        #G1: sample genotypes
        geno_sam <- function(p,n_fam){    #p:MAF, n_fam:number of families
            #sample genotypes for a single SNP
            geno_vec = rep(NA,3*n_fam)        #init
                    
            p_G = c(p^2,2*p*(1-p),(1-p)^2)        #G freq for aa,Aa,AA, assuming HWE
            geno_pa = matrix(sample(2:0,2*n_fam,replace=T,prob=p_G),nrow=2)        #genotype of parents
                    
                    
            geno_vec[c(1,2)+rep(seq(0,n_indPerFam*(n_fam-1),by=n_indPerFam),each=2)] = geno_pa
                        
            for(i in 1:n_chd){
                geno_child = apply(geno_pa,2,function(v) {col = p_trans_colidx[v[1]+1,v[2]+1];
                                    sample(2:0,1,prob=p_trans[,col])})
                geno_vec[seq(2+i,by=n_indPerFam,length.out=n_fam)] = geno_child
            }
                    
            geno_vec
        }
            
        if(!identical(seed,F)){
            set.seed(seed)
        }
            
        if(identical(p_a,F)){
            p_a = runif(n_snp,min=min(p_a_rg),max=max(p_a_rg))
        }
        geno = sapply(p_a,geno_sam,n_fam=n_fam)        #genotypes of SNPs
            
        #G2: sample marginal effect sizes for all SNPs
        if(identical(b_G,F)){
            b_G = runif(n_snp,min=min(b_G_rg),max=max(b_G_rg))
        }
                
        #G3: combine the marginal effect from all SNPs
        geno_marg = geno%*%b_G
            
        #G4: calculate additive var
        add_var = sum(2*p_a*(1-p_a)*b_G^2)
    }else{
        geno_marg = 0
        geno = NULL
    }
    
    if(!is.null(sig_a)){
        kmat <- kinship(gped[1])
        kmat <- 2*kmat
        geno_add = as.vector(t(mvrnorm(n_fam,rep(0,n_indPerFam),kmat*sig_a)))
        add_var = add_var+sig_a
    }else{
        geno_add = 0
    }
    
    #E
    if(rEnv){
        env = rep(rnorm(n_fam*n_indPerFam,sd=sqrt(sig_E)))
    }else{
        env = rep(rnorm(n_fam,sd=sqrt(sig_E)),each=n_indPerFam)
    }
    
    if(!is.null(sig_E2)){
        env2 = rep(rnorm(n_fam,sd=sqrt(sig_E2)),each=n_indPerFam)
    }else{
        env2 = 0
    }
    
    #GxE: the first n_isnp SNPs having interaction with E
    #1. sample effect sizes for these SNPs
    if(identical(b_GxE,F)){
        b_GxE = runif(n_isnp,min=min(b_GxE_rg),max=max(b_GxE_rg))
    }
            
    #2. combine the interation effect from all these SNPs
    gxe = as.matrix(geno[,1:n_isnp])%*%b_GxE*env

    #e
    err <- rnorm(n_indPerFam*n_fam,sd=sqrt(sig_e))
    
    #combine to be a dataset
    gdata <- data.frame(famid,id,fa,mo,sex,geno,geno_marg,geno_add,env,gxe,err,env2,stringsAsFactors = F)
    if(n_snp==1){names(gdata)[6] <- "X1"}   #if only one SNP
    
    gdata$ph_gei = with(gdata,geno_marg+geno_add+env+gxe+err+env2)        #pheno: G+E+interaction
    gdata$ph_g = with(gdata,geno_marg+geno_add+err+env2)        #only add effect
    gdata$ph_ge = with(gdata,geno_marg+geno_add+env+err+env2)        #add+env
    
    attr(gdata,"add_var") = add_var

    gdata
}



index_parents <- function(n_fam=100,n_chd=1){
    c(1,2)+rep(0:(n_fam-1)*(2+n_chd),each=2)
}
