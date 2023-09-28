# covid19_seroprevalence
Method for calculating seroprevalence using antibodies against two variants of covid-19. The method is described in article under review. The code for estimating seroprevalence against one variant, with or without the MRP, is directly from "Gelman A, Carpenter B. Bayesian analysis of tests with unknown specificity and sensitivity. Journal of the Royal Statistical Society: Series C (Applied Statistics). 2020;69(5):1269-83" This is the code in the files sero_prev.stan and sero_prev_mrp.stan. The remaing two stan files is the code for the new method taking two variants into account. 

example.R includes some simple example of how to run the stan code


