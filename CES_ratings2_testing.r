# We have 54 candidate presentations and 18 raters. Each rater watches
# 12 videos (=216 ratings) so that each presentation is rated by 4
# raters.

# Raters evaluate each presentation with two criteria: scientific
# quality and presentation quality. Raters are required to describe the
# reasons for evaluations for each presentation for letting raters think
# deeply.

# The exec plans to give prizes to four categories:
# A. students from LMIc / Global South = 4 out of 54.
# B. students from locals (Japan) = 9 out of 54.
# C. students giving poster presentations = 5 out of 54 (1 also belongs
# to category A)
# D. overall students  = 54.

# First, rank 54 students and pick up the top 3 and give them prizes.
# Second, pick up the best from LMIc/Global South, the best from locals,
# and the best from posters. Give them prizes.

# If necessary, we can provide the demographic information of raters:
# career stage (e.g., professor, assistant professor), sex, country of
# affiliations, ethnicity. We can also give you the demographic
# information of presenters (i.e., countries and gender).

###################################################################
# NEED sim_talks() function definition from CES_ratings1_develop.r
###################################################################

set.seed(1)

X <- sim_talks( 
    N=61 ,
    M=18 ,
    K=4 ,
    L=4 ,
    Q=2 ,
    RHO=matrix(c(1,0.5,0.5,1),2,2) ,
    verbose=FALSE )

dat <- X$dat

# count number of times each judge appears
table(dat$jid)

# run model
dat$weights <- rep(1,dat$Q)
m1 <- cstan( file="model1.stan" , data=dat , chains=3 , cores=3 )
# m0 <- cstan( file="model_null.stan" , data=dat , chains=3 , cores=3 )

precis(m1,3,pars=c("score","RHO_talks"))
precis(m1,3,pars=c("total_score"))

post <- extract.samples(m1)
score <- apply(post$score,2:3,mean)
talks <- X$truth$talks
N <- dat$n_talks

# total score
blank()
total_score_true <- sapply( 1:N , function(i) sum( dat$weights * talks[i,] ) )
total_score_est <- apply( post$total_score , 2,  mean )
rbPal <- colorRampPalette(c('black',2))
talks_col <- rbPal(10)[as.numeric(cut(total_score_true,breaks = 10))]
talks_col <- matrix( talks_col , ncol=2 )
plot( total_score_true , total_score_est , ylim=range(c(total_score_true,total_score_est)) , col=talks_col )
ci <- apply(post$total_score,2,PCI)
for ( i in 1:ncol(ci) ) lines( rep(total_score_true[i],2) , ci[,i] , lwd=0.8 , col=talks_col[i] )
abline( a=0 , b=1 , lty=2 , lwd=0.5 )

# compare ranks
rank_true <- rank( -total_score_true )
rank_est <- rank( -total_score_est )
plot( rank_true , rank_est )
cor( rank_true , rank_est )

# show difference between top-two ranked talks
ps1 <- post$total_score[, which(rank_est==3) ]
ps2 <- post$total_score[, which(rank_est==4) ]
dens( ps1 - ps2 , xlab="difference rank 3 and 4" )
abline( v=0 , lty=2 )

plot( precis( data.frame(ps1,ps2) ) , xlab="posterior total score" )
sum(ps1>ps2)/length(ps1)

# distinguishable talks by posterior diff threshold
delta <- 0.1
diffs <- matrix(NA,nrow=N,ncol=N)
dt <- diffs
for ( i in 1:N ) {
    for ( j in 1:N ) {
        dij <- post$total_score[ , which(rank_est==i) ] - post$total_score[ , which(rank_est==j) ]
        diffs[i,j] <- sum(dij > 0)/length(dij)
        dt[i,j] <- ifelse( diffs[i,j] > delta , 1 , 0 )
    }
}

round(diffs[1:10,1:10],2)

# features
blank(ex=2)

rbPal <- colorRampPalette(c('black',2))
talks_col <- rbPal(10)[as.numeric(cut(X$truth$talks,breaks = 10))]
talks_col <- matrix( talks_col , ncol=2 )

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

############################
# batch simulation

f <- function(N=54,M=18,K=4,L=4,Q=2) {
    X <- sim_talks( 
    N=N ,
    M=M ,
    K=K ,
    L=L ,
    Q=Q ,
    verbose=FALSE )

    dat <- X$dat

    dat$weights <- rep(1,dat$Q)
    m1 <- cstan( file="model1.stan" , data=dat , chains=3 , cores=3 , refresh=500 )

    post <- extract.samples(m1)
    talks <- X$truth$talks
    total_score_true <- sapply( 1:N , function(i) sum( dat$weights * talks[i,] ) )
    total_score_est <- apply( post$total_score , 2,  mean )
    rank_true <- rank( -total_score_true )
    rank_est <- rank( -total_score_est )
    
    return( cor( rank_true , rank_est ) )
}

n_sims <- 20
batch <- list()

for ( b in c(1,2,4,8) ) batch[[b]] <- replicate( n_sims , try( f(K=b) ) )


blank(w=1.3)
plot( NULL , xlim=c(0,1) , ylim=c(0.5,8.5) , xlab="rank correlation" , ylab="number of judges per item" , yaxt="n" )
axis( 2 , at=c(1,2,4,8) )
for ( b in c(1,2,4,8) ) points( batch[[b]] , rep(b,n_sims) , cex=1.6 , lwd=1.5 , col=2 )
