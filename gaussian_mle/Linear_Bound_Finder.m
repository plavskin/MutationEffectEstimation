function Linear_Bound_Finder(key_list, value_list)

	% EP 07-02-18

	% Takes list of parameter values, as well as cdf values based on
		% LRT corresponding to each of these parameter values, and fits
		% the parameter value corresponding to cdf_bound
    min_scaled_log_p_val = log(10^-50);
	tic;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get parameter values
    parameter_dict = containers.Map(key_list,value_list);

    cdf_bound = parameter_dict('cdf_bound');
    	% cdf_bound = 1-p_value
%    mle_param_val = parameter_dict('mle_param_val');
    parameter_vals = parameter_dict('parameter_values');
    cdf_vals = parameter_dict('cdf_vals');
    output_file = parameter_dict('output_file');
    linear_fit_file = parameter_dict('fit_file');
    pause_at_end = parameter_dict('pause_at_end');
    
    log_p_bound = log(1-cdf_bound);
    	% cdf_bound = 1-p_value
    log_p_vals = log(1-cdf_vals);
    % prevent errors due to parameter_dict('cdf_vals') == 1
    log_p_vals(log_p_vals < min_scaled_log_p_val) = min_scaled_log_p_val;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	optimal_coefficients = polyfit(parameter_vals,log_p_vals,1);
	
	m = optimal_coefficients(1);
	b = optimal_coefficients(2);
	linear_fit = [m,b];
	parameter_bounds = roots(linear_fit+[0,-log_p_bound]);
		% here, log_p_bound is subtracted from the value for b to
			% find the x-intercepts of the line at height log_p_bound

    
%    figure; plot(parameter_vals,log_p_vals,'ob'); hold on;
%    plot(parameter_vals,polyval(linear_fit,parameter_vals),'-r'); hold off;
            
            
	dlmwrite(output_file,parameter_bounds,'delimiter',',','precision',9);
	dlmwrite(linear_fit_file,linear_fit,'delimiter',',','precision',9);

	runtime = toc;

    if pause_at_end & runtime < 120
        pausetime=120-runtime;
        pause(pausetime)
    end

end