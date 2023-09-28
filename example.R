library(rstan)
library(data.table)


data <- data.frame(
  p=c(140), # Positive antigen 1
  p2=c(100), # Positive antigen 2
  p_12=c(80),# Positive antigen 1 and 2
  N=c(200) #number of samples
)

model_multi <- stan_model("sero_prev_multi.stan")
  data_stan <- list(
      y_A=data$p,
      y_B = data$p2,
      y_AB=data$p_12,
      N=data$N,
      y_sens_1=190,
      n_sens_1=200,
      y_spec_1=98,
      n_spec_1=100,
      sigma=0.05
    )
    iter = 3000
    chains = 3
    options(mc.cores = parallel::detectCores())
    fit <- sampling(model_multi, data = data_stan, iter=iter, chains=chains, verbose=F, control=list(max_treedepth=15))
    values <- extract(fit, "p")$p
    glue::glue("{format(round(mean(values)*100, 1), nsmall=1)} ({format(round(quantile(values, prob=0.025)*100, 1),nsmall=1)} - {format(round(quantile(values, prob=0.975)*100, 1), nsmall=1)})")
    values <- extract(fit, "p1")$p1
    glue::glue("{format(round(mean(values)*100, 1), nsmall=1)} ({format(round(quantile(values, prob=0.025)*100, 1),nsmall=1)} - {format(round(quantile(values, prob=0.975)*100, 1), nsmall=1)})")



y_spec=180
n_spec=200
y_sens=98
n_sens=100

n <- 1000
d_mrp <- data.frame(
    p1=sample(c(TRUE, FALSE), n, replace=TRUE),
    p2=sample(c(TRUE, TRUE, FALSE), n, replace=TRUE),
    kjonn=sample(c("m", "k"), n, replace=TRUE),
    age_group=sample(c(1,1,1,2,3,4,5), n, replace=TRUE),
    fylke=sample(c(1,2), n, replace=TRUE)

)
pop_group <- rep(1000, 20)

mrp_stan_model <- stan_model("sero_mrp.stan")

mrp_stan_multi_model <- stan_model("sero_mrp_multi.stan")



setDT(d_mrp)
mrp_data <- list(
    N=nrow(d_mrp),
    y=as.integer(d_mrp[, p1]),
    male=as.integer(d_mrp[, kjonn=="m"]),
    age=as.integer(factor(d_mrp[,age_group])),
    N_fylke=length(unique(d_mrp[, fylke])),
    fylke=as.integer(factor(d_mrp[, fylke])),
    y_spec=y_spec,
    n_spec=n_spec,
    y_sens=y_sens,
    n_sens=n_sens,
    J=2*5*length(unique(d_mrp[, fylke])),
    N_pop=pop_group,
    intercept_prior_mean=0,
    intercept_prior_scale=2.5,
    coef_prior_scale=0.5
)

mrp_data_multi <- list(
    N=nrow(d_mrp),
    y_1=as.integer(d_mrp[, p1]),
    y_2=as.integer(d_mrp[, p2]),
    y_12=as.integer(d_mrp[, p1==TRUE & p2==TRUE]),
    male=as.integer(d_mrp[, kjonn=="m"]),
    age=as.integer(factor(d_mrp[,age_group])),
    N_fylke=length(unique(d_mrp[, fylke])),
    fylke=as.integer(factor(d_mrp[, fylke])),
    y_spec=y_spec,
    n_spec=n_spec,
    y_sens=y_sens,
    n_sens=n_sens,
    J=2*5*length(unique(d_mrp[, fylke])),
    N_pop=pop_group,
    intercept_prior_mean=0,
    intercept_prior_scale=2.5,
    coef_prior_scale=0.5,
    sigma=0.05
)



iter = 2000
chains = 4
options(mc.cores = parallel::detectCores())
fit <- sampling(mrp_stan_model, data = mrp_data, iter=iter, chains=chains, verbose=F)
values <- extract(fit, "p_avg")$p_avg

fit_comb <- sampling(mrp_stan_multi_model, data = mrp_data_multi, iter=iter, chains=chains, verbose=F)

values_comb <- extract(fit_comb, "p_avg")$p_avg

glue::glue("One: {round(median(values)*100, 2)} ({round(quantile(values, prob=0.025)*100, 2)} - {round(quantile(values, prob=0.975)*100, 2)})")
glue::glue("COMB: {round(median(values_comb)*100, 2)} ({round(quantile(values_comb, prob=0.025)*100, 2)} - {round(quantile(values_comb, prob=0.975)*100, 2)})")
