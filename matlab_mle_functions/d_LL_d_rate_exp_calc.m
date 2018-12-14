function d_LL_d_rate = d_LL_d_rate_exp_calc(x, rate)
	% finds gradient of log likelihood of exponentially-distributed data (x)
		% relative to rate
	% x can be a scalar or a vector of data, mu must be either a scalars or
		% a vectors of same size as x
	d_LL_d_rate = (1./rate) - x;