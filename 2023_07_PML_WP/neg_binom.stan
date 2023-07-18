data {
int N;
int r;
int n[N];
}

parameters {
real<lower=0, upper = 1> beta;
}

model {
beta ~ beta(1,1);
n ~ neg_binomial(r, beta);
}

generated quantities {
 int yrep[N];
 yrep = neg_binomial_rng(rep_vector(r, N), rep_vector(beta, N));
}
