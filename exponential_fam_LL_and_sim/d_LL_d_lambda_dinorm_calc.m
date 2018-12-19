function d_LL_d_lambda = d_LL_d_lambda_dinorm_calc(x,lambda,mu_1,sigma_1,mu_2,sigma_2)
	% finds gradient of log likelihood of data drawn from a mixture of
		% two normal distributions (x) relative to lambda
		% (the proportion drawn from the first distribution)
	% x can be a scalar or a vector of data, the other parameters must
		% be either scalars or vectors of same size as x

	full_density = lambda*pdf('normal',x,mu_1,sigma_1)+...
        (1-lambda)*pdf('normal',x,mu_2,sigma_2);
    d_LL_d_lambda = (pdf('normal',x,mu_1,sigma_1)-pdf('normal',x,mu_2,sigma_2))./...
    	full_density;