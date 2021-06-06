# Produce posterior scores and ranks for CES 2021

library(rethinking)
library(cmdstanr)

# data need to be provided - placeholder here
dat_raw <- read.csv( "example_data.csv" )

dat_input <- list(
        N=nrow(dat_raw),
        n_talks=length(unique(dat_raw$talk)),
        M=length(unique(dat_raw$judge)),
        L=max(dat_raw[,3:4])-1,
        Q=2,
        y = dat_raw[,3:4],
        jid = dat_raw$judge,
        tid = dat_raw$talk,
        weights=c(1,1) )

# count number of times each judge appears
table(dat_input$jid)

# count number of times each talk appears (number of judges rating each talk)
table(dat_input$tid)

# run model
m_ces <- cstan( file="model1.stan" , data=dat_input , chains=4 , cores=4 , iter=4000 )

# check chains

x <- precis(m_ces,3)
plot( x$n_eff , x$Rhat4 , xlab="n_eff" , ylab="Rhat (v4)" , lwd=2 , col=col.alpha(2,0.8) )
# optional line to manually label outliers
# identify( x$n_eff , x$Rhat4 , labels=rownames(x) )

precis(m_ces,3,pars=c("RHO_talks"))
precis(m_ces,3,pars=c("total_score"))

# extract scores

post <- extract.samples(m_ces)
score <- apply(post$score,2:3,mean)
N <- dat_input$n_talks

# total score, sorted by rank
blank(w=2)

total_score_est <- apply( post$total_score , 2,  mean )
rank_est <- rank( -total_score_est )
rbPal <- colorRampPalette(c('black',2))
talks_col <- rbPal(10)[as.numeric(cut(total_score_est,breaks = 10))]
talks_col <- matrix( talks_col , ncol=2 )
ci <- apply(post$total_score,2,PCI)

plot( rank_est , total_score_est , ylim=range(ci) , col=talks_col )
for ( i in 1:ncol(ci) ) lines( rep(rank_est[i],2) , ci[,i] , lwd=0.8 , col=talks_col[i] )

# posterior rank distribution
blank(w=2)

rank_post <- sapply( 1:length(post$lp__) , function(i) rank( -post$total_score[i,] ) )
rank_mean <- apply( rank_post , 1 , mean )
o <- rank( rank_mean )
ci <- apply(rank_post,1,PCI)
plot( o , rank_mean , ylim=range(ci) , col=2 , xlab="mean rank" , ylab="posterior rank" )
for ( i in 1:ncol(ci) ) lines( rep(o[i],2) , ci[,i] , lwd=0.8 , col=2 )


# features
blank(ex=2)
post <- extract.samples(m_ces)
score <- apply(post$score,2:3,mean)
plot( score[,1:2] , xlab="post mean 1" , ylab="post mean 2" , lwd=1.5 , col="white" )
text( score[,1] , score[,2] , labels=1:N , cex=0.8 )
mtext( "features of each talk" )
