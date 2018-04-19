function d_LL_d_mu = d_LL_d_mu_norm_calc(x,mu,sigma)
	% finds gradient of log likelihood of normally-distributed data (x)
		% relative to mu (mean)
	% x can be a scalar or a vector of data, mu and sigma must be
		% either scalars or vectors of same size as x
	mean_diff = x - mu;
	d_LL_d_mu = (1./(sigma.^2)).*mean_diff;