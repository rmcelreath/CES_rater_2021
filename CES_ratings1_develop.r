# CES ratings model
# N talks
# M judges
# K judges per talk
# each talk receives a single category score 1 through 5 from K=3 random judges
# talks do not have "true" scores, but rather evoke a range of scores
# judges have unique mappings of underlying scores onto the 1-5 reporting scale, such that some judges assign only a very narrow range of the scale, and some have higher/lower means than others

# model framework: ordered logit with random cut-points (by judge) and different means for each talk (like an IRT)

# SIMULATION

sim_talks <- function( N=50 , M=10 , K=4 , L=4 , Q=2 , RHO=NULL , verbose=FALSE , flip=c() ) {

    require(rethinking)

    #N <- 50 # talks
    #M <- 10 # judges
    #K <- 4 # judges per talk
    #L <- 4 # number of cut-points (scale 1-5)
    #Q <- 2 # number of features each talk is rated for

    # simulate features of talk
    if ( Q==1 ) {
        talks <- rnorm(N,0,1.5)
        talks <- matrix( talks , nrow=N , ncol=1 )
    } else {
        # more than one feature, so talks needs to be a matrix
        if ( is.null(RHO) ) RHO <- rlkjcorr(1,Q,eta=2)
        talks <- rmvnorm2( N , rep(0,Q) , rep(1.5,Q) , Rho=RHO )
    }

    # assign judges to talks
    # each judge in M needs to rate at least K*N/M talks, so we can sample without replacement across talks, if we populate a sampling list with each judge ID K*N/M+1 times

    judge_counts <- rep( K*ceiling(N/M) , times=M )
    judge_by_talk <- matrix( NA , nrow=N , ncol=K )

        for ( i in 1:N ) {
            j2 <- 1:M
            j2 <- j2[judge_counts>0]
            if (length(j2)>1) {
                judge_by_talk[i,] <- sample( j2 , size=K , prob=exp(judge_counts[judge_counts>0]*2) )
            } else {
                judge_by_talk[i,] <- j2
            }
            # now deplete judges list by who was sampled
            judge_counts[judge_by_talk[i,]] <- judge_counts[judge_by_talk[i,]] - 1
            if ( verbose==TRUE ) print(judge_counts)
        }

    # simulate features of judges
    # each judge has a unique set of cut-points for each feature
    # but cut-points are correlated within judges
    if ( Q==0 ) {
        judges <- matrix(NA,nrow=M,ncol=L)
        for ( i in 1:M ) judges[i,] <- sort(rnorm(L,0,2))
    } else {
        # multiple features, so need an array
        # [judge,feature,cut-points]
        # want correlation between matching cut-points
        # with L=2 Q=2 we have e.g. 8 cut-points per judge [ a b A B ]
        # we want pairs with same letter to be correlated R
        # and we want each set [ab] [AB] to have own correlation structure
        #     a    b    A    B
        # a [ 1    (ab) R    0    ]
        # b [ (ab) 1    0    R    ]
        # A [ R    0    1    (AB) ]
        # B [ 0    R    (AB) 1    ]
        judges <- array( NA , dim=c(M,Q,L) )

        # algorithm 1: free
        if ( FALSE ) {
            for ( i in 1:M )
                for ( q in 1:Q )
                    judges[i,q,] <- sort(rnorm(L,0,2))
        } else {
        # algorithm 2: constrained
            for ( i in 1:M ) {
                for ( q in 1:Q )
                    judges[i,q,] <- sort(rnorm(L,0,2))
                # now add correlations to matched positions
                vL <- rnorm(L,0,1) # L match positions
                for ( l in 1:L )
                    judges[i,,l] <- judges[i,,l] + vL[1]
                # make sure still sorted
                for ( q in 1:Q ) judges[i,q,] <- sort(judges[i,q,])
            }#i
        }
    }

    if ( FALSE ) {
        # code to plot judgment distributions for judges
        blank2()
        plot( NULL , xlim=c(1,4) , ylim=c(0,1) , xlab="cut point" , ylab="prob" )
        for ( i in 1:M ) {
            p <- inv_logit( judges[i,] )
            lines( 1:4 , p )
        }
    }

    # now sim rating of each talk
    ratings_long <- matrix(NA,nrow=N*K,ncol=Q)
    jid <- rep(NA,N*K)
    tid <- rep(NA,N*K)
    r <- 1
    for ( i in 1:N ) { # judges
        for ( j in 1:K ) { # talks
            for ( q in 1:Q ) # features
                ratings_long[r,q] <- rordlogit( 1 , phi=talks[i,q] , a=judges[ judge_by_talk[i,j] , q , ] )
            # flipped scale?
            if ( i %in% flip ) {
                # flip the ratings for this judge
                ratings_long[r,] <- L + 2 - ratings_long[r,]
            }
            jid[r] <- judge_by_talk[i,j]
            tid[r] <- i
            r <- r + 1
        }#j
    }#i

    dat <- list(
        N=length(jid),
        n_talks=N,
        M=M,
        K=K,
        L=L,
        Q=Q,
        y = ratings_long,
        jid = jid,
        tid = tid )

    return( list( dat=dat , 
        truth=list(
            talks=talks,
            judges=judges,
            RHO=RHO) ) )

}#sim_talks

X <- sim_talks( verbose=TRUE )
dat <- X$dat
dat$weights <- rep(1,dat$Q)

# Stan model

m1 <- cstan( file="model1.stan" , data=dat , chains=3 , cores=3 , control=list(max_treedepth=15) )

precis(m1,3,pars=c("score","RHO_talks"))

# summarize RHO_cuts
post <- extract.samples(m1)
rc <- apply( post$RHO_cuts , 2:3 , mean )
round(rc,2)

##############
# m1 check

post <- extract.samples(m1)
score <- apply(post$score,2:3,mean)
rbPal <- colorRampPalette(c('black',2))
talks_col <- rbPal(10)[as.numeric(cut(talks,breaks = 10))]
talks_col <- matrix( talks_col , ncol=2 )

blank(ex=2)

par(mfrow=c(2,2))
for ( q in 1:2 ) {
    plot( talks[,q] , score[,q] , xlab="true score" , ylab="estimated score" , lwd=1.5 , col=talks_col[,q] )
    ci <- apply(post$score[,,q],2,PCI)
    for ( i in 1:ncol(ci) ) lines( c(talks[i,q],talks[i,q]) , ci[,i] , lwd=0.5 , col=talks_col[i,q] )
    mtext( concat("feature ",q) , 3 )
}#q
# scores for first two features
plot( talks[,1:2] , xlab="true feature 1" , ylab="true feature 2" , lwd=1.5 , col="white" )
text( talks[,1] , talks[,2] , labels=1:N , cex=0.8 )
plot( score[,1:2] , xlab="post mean 1" , ylab="post mean 2" , lwd=1.5 , col="white" )
text( score[,1] , score[,2] , labels=1:N , cex=0.8 )

#########
# example csv data for committee

X <- sim_talks( 
    N=54 ,
    M=18 ,
    K=4 ,
    L=4 ,
    Q=2 ,
    RHO=matrix(c(1,0.5,0.5,1),2,2) ,
    verbose=FALSE )

dat <- X$dat

out <- data.frame( talk=dat$tid , judge=dat$jid , feature1=dat$y[,1] , feature2=dat$y[,2] )

write.csv( out , file="example_data.csv" , row.names=FALSE )
