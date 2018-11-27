function Linear_Bound_Finder(key_list, value_list)

	% EP 07-02-18

	% Takes list of parameter values, as well as cdf values based on
		% LRT corresponding to each of these parameter values, and fits
		% the parameter value corresponding to cdf_bound
	tic;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get parameter values
    parameter_dict = containers.Map(key_list,value_list);

    cdf_bound = parameter_dict('cdf_bound');
    	% cdf_bound = 1-p_value
    mle_param_val = parameter_dict('mle_param_val');
    parameter_vals = parameter_dict('parameter_values');
    cdf_vals = parameter_dict('cdf_vals');
    output_file = parameter_dict('output_file');
    linear_fit_file = parameter_dict('fit_file');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	optimal_coefficients = polyfit(parameter_vals,cdf_vals,1);
	
	m = optimal_coefficients(1);
	b = optimal_coefficients(2);
	linear_fit = [m,b];
	parameter_bounds = roots(linear_fit+[0,-cdf_bound]);
		% here, cdf_bound is subtracted from the value for b to
			% find the x-intercepts of the line at height cdf_bound

    
%    figure; plot(parameter_vals,cdf_vals,'ob'); hold on;
%    plot(parameter_vals,polyval(linear_fit,parameter_vals),'-r'); hold off;
            
            
	dlmwrite(output_file,parameter_bounds,'delimiter',',','precision',9);
	dlmwrite(linear_fit_file,linear_fit,'delimiter',',','precision',9);

	runtime = toc;

    pause_at_end = true;

    if pause_at_end & runtime < 120
        pausetime=120-runtime;
        pause(pausetime)
    end

end