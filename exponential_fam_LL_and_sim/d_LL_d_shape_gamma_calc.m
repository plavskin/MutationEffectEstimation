function d_LL_d_shape = d_LL_d_shape_gamma_calc(x, shape, rate)
	% finds gradient of log likelihood of gamma-distributed data (x)
		% relative to shape
	% x can be a scalar or a vector of data, shape and rate must be either a scalars or
		% vectors of same size as x
	d_LL_d_shape = log(x) + log(rate) - psi(shape);