data {

     int y_spec_1;
     int n_spec_1;
     int y_sens_1;
     int n_sens_1;
     int N;
     int y_A;
     int y_B;
     int y_AB;
     real sigma;
}
parameters {
	   real<lower=0, upper=1> s1;
   	 real<lower=0, upper=1> s2;
   	 real<lower=0, upper=1> s3;
   	 
	   real<lower=0, upper=1> c1;
   	 real<lower=0, upper=1> c2;
   	 real<lower=0, upper=1> c3;
     real<lower=0, upper=1> p1  ;
     real<lower=0, upper=1> p2  ;
     real<lower=0, upper=1> p_12_a;
} transformed parameters{
  real p_12;
   real ps[2];
   real pmins[2];
   ps[1] = p1;
   ps[2] = p2; 
   pmins[1] = 0;
   pmins[2] = p1 + p2 - 1;
   p_12 = max(pmins) + p_12_a*(min(ps) - max(pmins));

}
model {

      y_A ~ binomial(N,p1*s1 + (1-p1)*(1-c1));
      y_B ~ binomial(N,p2*s2 + (1-p2)*(1-c2));
      y_AB ~ binomial(N,p_12*(s3) + (1-p_12)*((1-c3)));

      y_spec_1 ~ binomial(n_spec_1, c1);
      y_sens_1 ~ binomial(n_sens_1, s1);
      s2 ~ normal(s1, sigma);
      s3 ~ normal(s1, sigma);
      c2 ~ normal(c1, sigma);
      c3 ~ normal(c1, sigma);
      }
generated quantities{
  real p = p1 + p2 - p_12;
}