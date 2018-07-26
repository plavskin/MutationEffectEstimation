function Quadratic_Bound_Finder(key_list, value_list)

	% EP 07-02-18

	% Takes list of parameter values, as well as cdf values based on
		% LRT corresponding to each of these parameter values, and fits
		% the parameter value corresponding to cdf_bound
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get parameter values
    parameter_dict = containers.Map(key_list,value_list);

    cdf_bound = parameter_dict('cdf_bound');
    	% cdf_bound = 1-p_value
    mle_param_val = parameter_dict('mle_param_val');
    parameter_vals = parameter_dict('parameter_values');
    cdf_vals = parameter_dict('cdf_vals');
    output_file = parameter_dict('output_file');
    quad_fit_file = parameter_dict('fit_file');
    profile_side = parameter_dict('profile_side');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % quadratic_fit_fun reframes quadratic function to allow bounds to
    	% be placed on coefficients
	quadratic_fit_fun = @(coeff_vector,x_vals,y_vals) y_vals-...
		(coeff_vector(1)*x_vals.^2-...
		2*coeff_vector(1)*coeff_vector(2)*x_vals+coeff_vector(3));
		% coeff_vector = [a, 2*a*b, c]

	% initialize starting coefficients and set constrains
	% a must be negative
	% x_max > -b/2a or x_max < -b/2a, depending on direction of CI bound
		% this is the axis of the parabola; this requirement means that
			% all parameter_vals are fitted on a single,
			% 'ascending' side of the parabola
	direct_fit_coefficients = polyfit(parameter_vals,cdf_vals,2);
	starting_coefficients = NaN([1,3]);
	starting_coefficients(1) = min(direct_fit_coefficients(1),(mle_param_val-10^-50));
		% make sure a is negative, i.e. parabola 'faces' down
%	if min(parameter_vals)<mle_param_val
	if strcmp(profile_side,'lower')
		current_lb = [-Inf,-Inf,-Inf];
		current_ub = [0,min(parameter_vals),Inf];
		starting_coefficients(2) = max(-direct_fit_coefficients(2)/(2*starting_coefficients(1)),min(parameter_vals));
			% second coefficient is -b/2a
%	else
	elseif strcmp(profile_side,'upper')
		current_lb = [-Inf,max(parameter_vals),-Inf];
		current_ub = [0,Inf,Inf];
		starting_coefficients(2) = min(-direct_fit_coefficients(2)/(2*starting_coefficients(1)),max(parameter_vals));
	else
		disp('Profile side selected:')
		disp(profile_side)
		error('Invalid profile_side')
	end
	starting_coefficients(3) = direct_fit_coefficients(3);

	% Fit quadratic with constraints above
	optimal_coefficients = lsqnonlin(@(coeffs) quadratic_fit_fun(coeffs,parameter_vals,cdf_vals),...
		starting_coefficients,current_lb,current_ub);

	a = optimal_coefficients(1);
	b = -optimal_coefficients(2)*2*a;
	c = optimal_coefficients(3);
	quadratic_fit = [a,b,c];
    
%    x_plot_vals = linspace(min(parameter_vals),max(parameter_vals),50);
%    figure; plot(parameter_vals,cdf_vals,'ob'); hold on;
%    plot(x_plot_vals,polyval(quadratic_fit,x_plot_vals),'-r'); hold off;
    
	possible_bounds = roots(quadratic_fit+[0,0,-cdf_bound]);
		% here, cdf_bound is subtracted from the value for c to
			% find the x-intercepts of the parabola at height cdf_bound

	% find the point out of possible_bounds that is closest to mle_param_val:
		% this is the intercept on the correct side of the parabola

	[~,bound_index] = min(abs(possible_bounds-mle_param_val));
	parameter_bound = possible_bounds(bound_index);

	dlmwrite(output_file,parameter_bound,'delimiter',',','precision',9);
	dlmwrite(quad_fit_file,quadratic_fit,'delimiter',',','precision',9);

end