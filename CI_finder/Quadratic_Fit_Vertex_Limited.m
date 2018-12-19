function y_diff = Quadratic_Fit_Vertex_Limited(coeff_vector,x_vals,y_vals)

	% EP 18-11-27
	% Calculates difference between y_vals and quadratic(x_vals), where
		% the parameters in coeff_vector are functions of a, b, and c,
		% so that constraints can be placed on the x- and y- position of
		% the vertex, and the direction of the parabola, by placing
		% constraints on coeff_vector

	a = coeff_vector(1);
	b = -2*coeff_vector(2)*a;
	c = coeff_vector(3)+(b^2)/(4*a);

	y_fit = a * x_vals.^2 + b * x_vals + c;
	y_diff = y_vals - y_fit;