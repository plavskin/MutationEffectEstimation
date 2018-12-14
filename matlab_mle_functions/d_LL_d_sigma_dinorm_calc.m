function d_LL_d_sigma_1 = d_LL_d_sigma_dinorm_calc(x,lambda,mu_1,sigma_1,mu_2,sigma_2)
	% finds gradient of log likelihood of data drawn from a mixture of
		% two normal distributions (x) relative to sigma_1
		% (the sd of the first distribution)
	% x can be a scalar or a vector of data, the other parameters must
		% be either scalars or vectors of same size as x
	
	full_density = lambda*pdf('normal',x,mu_1,sigma_1)+...
        (1-lambda)*pdf('normal',x,mu_2,sigma_2);
    single_dist_density = pdf('normal',x,mu_1,sigma_1);
    d_LL_single_dist_d_sigma = d_LL_d_sigma_norm_calc(x, mu_1, sigma_1);
    d_LL_d_sigma_1 = lambda*d_LL_single_dist_d_sigma.*single_dist_density./full_density;
