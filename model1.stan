// model with multiple rating features
data{
    int N;       // number of talk/judge combinations
    int n_talks; // number of unique talks in the sample
    int M;       // number of unique judges in sample
    int L;       // number of cut-points (scale levels minus one)
    int Q;       // number of features each talk is rated on
    int jid[N];  // judge ids
    int y[N,Q];  // ratings with features in columns
    int tid[N];  // talk ids
    vector[Q] weights; // weight of each feature in final score
}
parameters{
    // model represents cut-points using positive steps between them, so the ordering is done right up front
    // judge parameters - means
    vector[L] dcuts[Q]; // logit scale offsets from previous cuts - [1,] are absolute logit, not offset
    cholesky_factor_corr[L*Q] L_cuts;
    vector<lower=0>[L*Q] sigma_cuts;
    matrix[ L*Q , M ] z_cuts; // log scale offsets from dcuts vector for each judge

    // talk parameters
    cholesky_factor_corr[Q] L_talks;
    matrix[ Q , n_talks ] z;
    vector<lower=0>[Q] tau; // scale for talk features
}
transformed parameters{
    vector[L] cuts_jid[M,Q];
    matrix[n_talks,Q] score;

    // compute judge cut-points on proper scale
    {
      matrix[M,L*Q] V;
      V = (diag_pre_multiply( sigma_cuts , L_cuts ) * z_cuts)';
      for ( i in 1:M ) {
          for ( q in 1:Q ) {
            int offset = (q-1)*L;
            // first cut-point is just logit scale with offset for judge - no transformations
            cuts_jid[i,q,1] = dcuts[q,1] + V[ i , 1 + offset ];
            // other cut-points are offsets from first one, with exp transforms to enforce positive offsets
            for ( j in 2:L )
                cuts_jid[i,q,j] = cuts_jid[i,q,j-1] + exp( dcuts[q,j] + V[i,j+offset] );
          }//q
      }//i
    }

    // compute talk scores
    score = (diag_pre_multiply( tau , L_talks ) * z)';
}
model{
    // priors
    tau ~ exponential( 1 );
    sigma_cuts ~ exponential(1);
    to_vector(z) ~ normal( 0 , 1 );
    to_vector(z_cuts) ~ normal(0,1);
    L_cuts ~ lkj_corr_cholesky(2);
    L_talks ~ lkj_corr_cholesky(2);

    // priors for cut-point means
    // first of each set absolute logit scale, rest are differences on log scale
    to_vector(dcuts[ 1:Q , 1 ]) ~ normal(-1,1);
    for ( j in 2:L )
        to_vector(dcuts[ 1:Q , j ]) ~ normal(0,0.5);

    // prob of ratings
    for ( i in 1:N ) {
        real phi = 0;
        for ( q in 1:Q ) {
            phi = score[tid[i],q];
            y[i,q] ~ ordered_logistic( phi , cuts_jid[ jid[i] , q ] );
        }//q
    }//i
}
generated quantities {
    matrix[Q,Q] RHO_talks;
    matrix[L*Q,L*Q] RHO_cuts;
    vector[n_talks] total_score;
    for ( i in 1:n_talks ) total_score[i] = dot_product( weights , score[i,] );
    RHO_talks = multiply_lower_tri_self_transpose(L_talks);
    RHO_cuts = multiply_lower_tri_self_transpose(L_cuts);
}
