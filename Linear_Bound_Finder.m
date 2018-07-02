function Linear_Bound_Finder(key_list, value_list)

	% EP 07-02-18

	% Takes list of parameter values, as well as cdf values based on
		% LRT corresponding to each of these parameter values, and fits
		% the parameter value corresponding to 1-cutoff_pval
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get parameter values
    parameter_dict = containers.Map(key_list,value_list);

    cutoff_pval = str2num(parameter_dict('cutoff_pval'));
    mle_param_val = str2num(parameter_dict('mle_param_val'));
    parameter_vals = parameter_dict('parameter_values');
    cdf_vals = parameter_dict('cdf_vals');
    output_file = parameter_dict('combined_max_array');
    output_file = parameter_dict('output_file');
    linear_fit_file = parameter_dict('linear_fit_file');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	optimal_coefficients = polyfit(parameter_vals,cdf_vals,1);
	
	m = optimal_coefficients(1);
	b = optimal_coefficients(2);
	linear_fit = [m,b];
	cdf_bounds = roots(linear_fit+[0,-(1-cutoff_pval)]);
		% here, (1-cutoff_pval) is subtracted from the value for b to
			% find the x-intercepts of the line at height
			% 1-cutoff_pval

	dlmwrite(output_file,[cdf_bound],'delimiter',',','precision',9);
	dlmwrite(linear_fit_file,linear_fit,'delimiter',',','precision',9);

end