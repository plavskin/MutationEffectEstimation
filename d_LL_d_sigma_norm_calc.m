function d_LL_d_sigma = d_LL_d_sigma_norm_calc(x,mu,sigma)
	% finds gradient of log likelihood of normally-distributed data (x)
		% relative to sigma (s.d.)
	% x can be a scalar or a vector of data, mu and sigma must be
		% either scalars or vectors of same size as x
	mean_diff = x - mu;
	d_LL_d_sigma = (-1./sigma).*(1-(1./(sigma.^2)).*mean_diff.^2);