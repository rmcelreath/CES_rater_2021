# Produce posterior scores and ranks for CES 2021

library(rethinking)
library(cmdstanr)

# raw anonymized ratings
dat_raw <- read.csv( "CES_report/raw/data_v2.csv" )

# remove talk 61 (was not uploaded, cannot be rated)
w61 <- which( dat_raw$talk==61 )
dat_raw <- dat_raw[ -w61 , ]

# prep data input
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

# extract scores
post <- extract.samples(m_ces)
score <- apply(post$score,2:3,mean)
total_score_est <- apply( post$total_score , 2,  mean )
N <- dat_input$n_talks

# compute ranks
rank_post <- sapply( 1:length(post$lp__) , function(i) rank( -post$total_score[i,] ) )
rank_mean <- apply( rank_post , 1 , mean )

# posterior means of each feature
features <- apply(post$score,2:3,mean)

# report
report <- data.frame( talk=1:dat_input$n_talks , rank=NA , score=NA , feature1=NA , feature2=NA )
report$rank <- rank_mean
report$score <- total_score_est
report$feature1 <- features[,1]
report$feature2 <- features[,2]
# sort report by rank
report_o <- report[ order( rank_mean ) , ]
write.csv( report_o , file="report.csv" , row.names=FALSE )

# save RData - reload with load("CES2021_talks.RData")
save.image(file = "CES2021_talks.RData")


#################################################################
##### PLOTS #####

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
blank()
score_naive <- sapply( 1:60 , function(i) mean( unlist(dat_input$y[dat_input$tid==i,]) ) )
ci <- apply(rank_post,1,PCI)
plot( score_naive , rank_mean , ylim=range(ci) , col=2 , xlab="naive average score" , ylab="posterior rank" )
for ( i in 1:ncol(ci) ) lines( rep(score_naive[i],2) , ci[,i] , lwd=0.8 , col=2 )

# features
blank()
post <- extract.samples(m_ces)
score <- apply(post$score,2:3,mean)
plot( score[,1:2] , xlab="posterior mean feature 1" , ylab="posterior mean feature 2" , lwd=2 , col=2 , cex=1.5 )
# text( score[,1] , score[,2] , labels=1:N , cex=0.8 )
mtext( "features of each talk" )

######################
# supplemental plots

# naive average score approach
score_naive <- sapply( 1:60 , function(i) mean( unlist(dat_input$y[dat_input$tid==i,]) ) )
plot( score_naive , total_score_est , xlab="naive score" , ylab="posterior mean score" , lwd=2 , col=2 )

plot( rank( -score_naive ) , rank_mean , xlab="naive rank" , ylab="posterior mean rank" , lwd=2 , col=2 )

# distribution of raw scores
simplehist(dat_input$y , xlab="score" , col=c(1,2) )
text( 2 , 50 , "feature 1" )
text( 2 , 46 , "feature 2" , col=2 )

# show judge variation
# distribution of raw scores for specific judge
simplehist( dat_input$y[ dat_input$jid==1 , ] , xlab="score" , col=c(1,2) )
text( 2 , 50 , "feature 1" )
text( 2 , 46 , "feature 2" , col=2 )

# average judge cut-points
dcuts <- apply( post$dcuts , 2:3 , mean )
# covert to logit scale cut-points
lcuts <- matrix( NA , nrow=2 , ncol=6 )
lcuts[,1] <- dcuts[,1] # these are already logit scale
for ( j in 2:6 ) {
        lcuts[,j] <- lcuts[,j-1] + exp( dcuts[,j] )
}
mu_cuts <- lcuts
#plot
blank()
plot( NULL , xlim=c(1,7) , ylim=c(0,1) , xlab="value" , ylab="cumulative probability" )
lines( 1:7 , c( inv_logit(lcuts[1,]) , 1 ) , lwd=2 )
lines( 1:7 , c( inv_logit(lcuts[2,]) , 1 ) , lwd=2 , col=2 )

# sample judge cut-points
blank()

plot( NULL , xlim=c(1,7) , ylim=c(0,1) , xlab="value" , ylab="cumulative probability" )
for ( i in 1:16 ) {
        lcuts <- apply( post$cuts_jid[,i,,] , 2:3 , mean )
        lines( 1:7 , c( inv_logit(lcuts[1,]) , 1 ) , lwd=1 )
        #lines( 1:7 , c( inv_logit(lcuts[2,]) , 1 ) , lwd=1 , col=2 )
}

# again by individual judge
blank(ex=3.2)
par(mfrow=c(4,4))
for ( i in 1:16 ) {
plot( NULL , xlim=c(1,7) , ylim=c(0,1) , xlab="" , ylab="" )
lcuts <- apply( post$cuts_jid[,i,,] , 2:3 , mean )
lines( 1:7 , c( inv_logit(lcuts[1,]) , 1 ) , lwd=2 )
lines( 1:7 , c( inv_logit(lcuts[2,]) , 1 ) , lwd=2 , col=2 )
mtext( concat("Rater ",i) )
lines( 1:7 , c( inv_logit(mu_cuts[1,]) , 1 ) , lwd=0.5 , lty=2 )
#lines( 1:7 , c( inv_logit(mu_cuts[2,]) , 1 ) , lwd=0.5 , lty=2 , col=2 )
}

# illustrate posterior rank distribution using top two talks (38 and 21)
simplehist( rank_post[38,] , xlim=c(1,60) , xlab="rank" , col=2 )
simplehist( rank_post[21,] , ylim=c(0,3500) , xlim=c(1,60) , xlab="rank" )

dr <- rank_post[38,]-rank_post[21,]
simplehist( dr , xlab="rank difference" )
sum(dr<0)/8000
