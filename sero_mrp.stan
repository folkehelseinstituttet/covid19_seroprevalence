data {
     int<lower = 0> N; // number of tests in the sample (3330 for Santa Clara)
     int<lower = 0, upper = 1> y[N]; // 1 if positive, 0 if negative
     vector<lower = 0, upper = 1>[N] male; // 0 if female, 1 if male
     int<lower = 1, upper = 5> age[N]; 
     int<lower = 0> N_fylke; // number of fylke codes (58 in this case)
     int<lower = 1, upper = N_fylke> fylke[N]; // fylke codes 1 through 58
     int<lower = 0> y_spec;
     int<lower = 0> n_spec;
     int<lower = 0> y_sens;
     int<lower = 0> n_sens;
     int<lower = 0> J; // number of population cells, J = 2*4*4*58
     vector<lower = 0>[J] N_pop; // population sizes for poststratification
     real intercept_prior_mean;
     real<lower = 0> intercept_prior_scale;
     real<lower = 0> coef_prior_scale;
}
parameters {
	   real<upper=1> spec;
	   real sens;
	   vector[2] b; // intercept, coef for male
	   real<lower = 0> sigma_age;
	   real<lower = 0> sigma_fylke;
	   vector<multiplier = sigma_age>[5] a_age; // varying intercepts for age category
	   vector<multiplier = sigma_fylke>[N_fylke] a_fylke; // varying intercepts for fylke code
}
model {
      vector[N] p = inv_logit(b[1]
      		+ b[2] * male
		+ a_age[age]
		+ a_fylke[fylke]);
       vector[N] p_sample = p * sens + (1 - p) * (1 - spec);
       y ~ bernoulli(p_sample);
       y_spec ~ binomial(n_spec, spec);
       y_sens ~ binomial(n_sens, sens);
       a_age ~ normal(0, sigma_age);
       a_fylke ~ normal(0, sigma_fylke);
       // prior on centered intercept
       b[1] + b[2] * mean(male) 
       	    ~ normal(intercept_prior_mean, intercept_prior_scale);
       b[2] ~ normal(0, coef_prior_scale);
       sigma_age ~ normal(0, coef_prior_scale);
       sigma_fylke ~ normal(0, coef_prior_scale);
}
generated quantities {
	  real p_avg;
	  vector[J] p_pop; // population prevalence in the J poststratification cells
	  int count;
	  count = 1;
	  for (i_fylke in 1:N_fylke) {
	      for (i_age in 1:5) {
	           for (i_male in 0:1) {
		      	  p_pop[count] = inv_logit(b[1]
			  	       + b[2] * i_male
				       + a_age[i_age]
				       + a_fylke[i_fylke]);
				       count += 1;
	           }

 		}
	}
	p_avg = sum(N_pop .* p_pop) / sum(N_pop);
}
