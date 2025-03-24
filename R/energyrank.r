rank_energy_stat <- function(x,sampletimes,dis.beta=1){
    n <- nrow(x);d = ncol(x)
    if(n<30){
        min_size = 2; #default min_size = 30
    }else{
        min_size = 10
    }
    #############
    #rank-permute#
    #############
    grid <- matrix(halton(n,dim=d),ncol=d);
    distmat <- Rfast::dista(xnew=x,x=grid)^2;
    solution<-transport::transport(rep(1/n, n), rep(1/n, n), p = 2,
                                method = 'networkflow',
                                costm = distmat); ###O(n^3)
    rm(distmat)
    gc()
    co_rank <- matrix(grid[solution[,2],],ncol = d);
    D_rank = (Rfast::Dist(co_rank))^{dis.beta}
    #normalization for Y(t)#
    stats1_rank_eig <- split_re(s_ = 1,e_ = n,D_=D_rank,min_size_ = min_size)
    #
    stats_rank_permut <- split_re(s_ = 1,e_ = n,D_= D_rank,min_size_ = min_size)
    stat_rank_permut <- stats_rank_permut[2]#max(stats)
    loc_rank_permut <- stats_rank_permut[1]#which.max(stats)
    stat0_rank_permut <- sapply(1:sampletimes, function(x){
        shuffle <- sample(1:nrow(D_rank));
        stat <- split_re(s_ = 1,e_ = n,D_=D_rank[shuffle,shuffle],min_size_ = min_size)
        return(stat[2])
    })
    eig_Hn_rank <- eigencompute(D_rank,n)

    p.value_rank_permut <- ( sum( stat0_rank_permut >= stat_rank_permut ) + 1 ) / (sampletimes + 1)
    ####
    #test by asymptotic approximation#
    if(n>=50){
        eig_num = 50;#50 or bigger, can compare#
    }else{
        eig_num=n
    }
    #
    time_ind <- seq(0,1,length=1000)
    ttt = time_ind*(1-time_ind)
    eigs_rank <- eig_Hn_rank$values[order(abs(eig_Hn_rank$values),decreasing = T)[1:eig_num]]
    ##
    AA<-lapply(1:sampletimes, function(NN){
        rr <- sapply(1:eig_num,function(x){
            return(sde::BBridge(x=0, y=0, t0=0, T=1, N=999)^2)
        })
        tmp <- t(ttt - rr)
        Yt_rank <- colsums(abs(eigs_rank)*tmp)
        return(max(abs(Yt_rank)))#,max(abs(Yt1))))
    })
    AA <- do.call(c,AA)
    ###rank_eig###
    stat0_rank_eig <- AA;#cutoffs;
    stat_rank_eig <- stats1_rank_eig[2];
    loc_rank_eig <- stats1_rank_eig[1];#which.max(stats1)
    p.value_rank_eig <- ( sum( stat0_rank_eig >= stat_rank_eig ) + 1 ) / (sampletimes + 1)
    return(
        c(p.value_rank_permut,p.value_rank_eig,
            loc_rank_permut,loc_rank_eig)
    )
}

energy_stat <- function(x,sampletimes,dis.beta=1){
    n <- nrow(x);d = ncol(x)
    D = (Rfast::Dist(x))^{dis.beta}
    if(n<30){
        min_size = 2; #default min_size = 30
    }else{
        min_size = 10
    }
    #normalization for Y(t)#
    stats1_ori_eig <- split_re(s_ = 1,e_ = n,D_=D,min_size_ = min_size)
    ##
    stats_ori_permut <- split_re1(s_ = 1,e_ = n,D_= D,min_size_ = min_size)
    stat_ori_permut <- stats_ori_permut[2]#max(stats)
    loc_ori_permut<- stats_ori_permut[1]#which.max(stats)
    #
    stat0_ori_permut <- sapply(1:sampletimes, function(x){
        shuffle <- sample(1:n);
        stat <- split_re1(s_ = 1,e_ = n,D_=D[shuffle,shuffle],min_size_ = min_size)
        return(stat[2])
    })
    eig_Hn <- eigencompute(D,n)
    p.value_ori_permut <- ( sum( stat0_ori_permut >= stat_ori_permut ) + 1 ) / (sampletimes + 1)

    ####
    #test by asymptotic approximation#
    if(n>=50){
        eig_num = 50;#50 or bigger#
    }else{
        eig_num=n
    }
                         #
    time_ind <- seq(0,1,length=1000)
    ttt = time_ind*(1-time_ind)
    eigs <- eig_Hn$values[order(abs(eig_Hn$values),decreasing = T)[1:eig_num]]
    ##
    AA<-lapply(1:sampletimes, function(NN){
        rr <- sapply(1:eig_num,function(x){
            return(sde::BBridge(x=0, y=0, t0=0, T=1, N=999)^2)
        })
        tmp <- t(ttt - rr)
        Yt <- colsums(abs(eigs)*tmp)
        return(max(abs(Yt)))#,max(abs(Yt1))))
    })
    AA <- do.call(c,AA)
    ###ori_eig###
    stat0_ori_eig <- AA;
    stat_ori_eig <- stats1_ori_eig[2];
    loc_ori_eig <- stats1_ori_eig[1];#which.max(stats1)
    p.value_ori_eig <- ( sum( stat0_ori_eig >= stat_ori_eig ) + 1 ) / (sampletimes + 1)
    return(
        c(p.value_ori_permut,p.value_ori_eig,
            loc_ori_permut,loc_ori_eig)
    )
}

eigencompute <- function(D,n){
    cc <- diag(n)-1/n*matrix(rep(1,n),ncol=1)%*%t(matrix(rep(1,n),ncol=1))
    #centering the kernel matrix D#
    Hn <- cc%*%D%*%t(cc)/n
    eig_Hn <- eigen(Hn,only.values = T);
    return(eig_Hn)
}

