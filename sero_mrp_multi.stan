data {
     int<lower = 0> N; // number of tests in the sample 
     int<lower = 0, upper = 1> y_1[N]; // 1 if positive, 0 if negative for test 1
     int<lower = 0, upper = 1> y_2[N]; // 1 if positive, 0 if negative for test 2
     int<lower = 0, upper = 1> y_12[N]; // 1 if positive, 0 if negative for test 1 and 2
     vector<lower = 0, upper = 1>[N] male; // 0 if female, 1 if male
     int<lower = 1, upper = 5> age[N]; 
     int<lower = 0> N_fylke; // number of fylke codes 
     int<lower = 1, upper = N_fylke> fylke[N]; // fylke 
     int<lower = 0> y_spec;
     int<lower = 0> n_spec;
     int<lower = 0> y_sens;
     int<lower = 0> n_sens;
     int<lower = 0> J; // number of population cells, 
     vector<lower = 0>[J] N_pop; // population sizes for poststratification
     real intercept_prior_mean;
     real<lower = 0> intercept_prior_scale;
     real<lower = 0> coef_prior_scale;
     real sigma;
}
parameters {
	   real<upper=1> spec_1;
	   real<upper=1> sens_1;
	   real<upper=1> spec_2;
	   real<upper=1> sens_2;
	   real<upper=1> spec_3;
	   real<upper=1> sens_3;
	   vector[2] b_1;
      vector[2] b_2;
      vector[2] b_3;
       real<lower = 0> sigma_age;
	   real<lower = 0> sigma_fylke;
	   vector<multiplier = sigma_age>[5] a_age_1; // varying intercepts for age category
	   vector<multiplier = sigma_fylke>[N_fylke] a_fylke_1; // varying intercepts for fylke code
       vector<multiplier = sigma_age>[5] a_age_2; // varying intercepts for age category
	   vector<multiplier = sigma_fylke>[N_fylke] a_fylke_2; // varying intercepts for fylke code
       vector<multiplier = sigma_age>[5] a_age_3; // varying intercepts for age category
	   vector<multiplier = sigma_fylke>[N_fylke] a_fylke_3; // varying intercepts for fylke code
       
}
model {
     vector[N] p_12;
       real ps[2];
       real pmins[2];
      
      vector[N] p_1 = inv_logit(b_1[1]
      		+ b_1[2] * male
		+ a_age_1[age]
		+ a_fylke_1[fylke]);
       vector[N] p_sample_1 = p_1 * sens_1 + (1 - p_1) * (1 - spec_1);
     
       vector[N] p_2 = inv_logit(b_2[1]
      		+ b_2[2] * male
		+ a_age_2[age]
		+ a_fylke_2[fylke]);
       vector[N] p_sample_2 = p_2 * sens_2 + (1 - p_2) * (1 - spec_2);
     
        vector[N] p_3 = inv_logit(b_3[1]
      		+ b_3[2] * male
		+ a_age_3[age]
		+ a_fylke_3[fylke]);
       vector[N] p_sample_3 = p_3 * sens_3 + (1 - p_3) * (1 - spec_3);
      for(i in 1:N){
         ps[1] = p_1[i];
         ps[2] = p_2[i]; 
         pmins[1] = 0;
         pmins[2] = p_1[i] + p_2[i] - 1;
         p_12[i] = max(pmins) + p_3[i]*(min(ps) - max(pmins));
      }
      y_2 ~ bernoulli(p_sample_2);
      y_1 ~ bernoulli(p_sample_1);
      y_12 ~ bernoulli(p_12);
      y_spec ~ binomial(n_spec, spec_1);
      y_sens ~ binomial(n_sens, sens_1);
      sens_2 ~ normal(sens_1, sigma);
      sens_3 ~ normal(sens_1, sigma);
      spec_2 ~ normal(spec_1, sigma);
      spec_3 ~ normal(spec_1, sigma);
      a_age_1 ~ normal(0, sigma_age);
      a_fylke_1 ~ normal(0, sigma_fylke);
       // prior on centered intercept
      b_1[1] + b_1[2] * mean(male) 
       	    ~ normal(intercept_prior_mean, intercept_prior_scale);
      b_1[2] ~ normal(0, coef_prior_scale);
      a_age_1 ~ normal(0, sigma_age);
      a_fylke_1 ~ normal(0, sigma_fylke);
       // prior on centered intercept
      b_2[1] + b_2[2] * mean(male) 
       	    ~ normal(intercept_prior_mean, intercept_prior_scale);
      b_2[2] ~ normal(0, coef_prior_scale);
      a_age_2 ~ normal(0, sigma_age);
      a_fylke_2 ~ normal(0, sigma_fylke);
       // prior on centered intercept
      b_3[1] + b_3[2] * mean(male) 
       	    ~ normal(intercept_prior_mean, intercept_prior_scale);
      b_3[2] ~ normal(0, coef_prior_scale);
      sigma_age ~ normal(0, coef_prior_scale);
      sigma_fylke ~ normal(0, coef_prior_scale);
}
generated quantities {
    real ps[2];
       real pmins[2];
	  real p_avg;
	  vector[J] p_pop; 
      vector[J] p_pop_1;
      vector[J] p_pop_2;
      vector[J] p_pop_3;
      vector[J] p_pop_12;
	  int count;
	  count = 1;
	  for (i_fylke in 1:N_fylke) {
	      for (i_age in 1:5) {
	           for (i_male in 0:1) {
		      	  p_pop_1[count] = inv_logit(b_1[1]
			  	       + b_1[2] * i_male
				       + a_age_1[i_age]
				       + a_fylke_1[i_fylke]);
                  p_pop_2[count] = inv_logit(b_1[1]
			  	       + b_1[2] * i_male
				       + a_age_2[i_age]
				       + a_fylke_2[i_fylke]);
                  p_pop_3[count] = inv_logit(b_3[1]
			  	       + b_3[2] * i_male
				       + a_age_3[i_age]
				       + a_fylke_3[i_fylke]);
                    ps[1] = p_pop_1[count];
                    ps[2] = p_pop_2[count]; 
                    pmins[1] = 0;
                    pmins[2] = p_pop_1[count] + p_pop_2[count] - 1;
                       p_pop_12[count] = max(pmins) + p_pop_3[count]*(min(ps) - max(pmins));
                        p_pop[count] = p_pop_1[count] + p_pop_2[count] - p_pop_12[count];
				       count += 1;
	           }

 		}
	}
	p_avg = sum(N_pop .* p_pop) / sum(N_pop);
}
