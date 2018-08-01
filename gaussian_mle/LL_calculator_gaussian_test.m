function [neg_combined_LL,global_gradient_vector_partial] = LL_calculator_gaussian_test(param_vals_partial,...
    global_fixed_parameter_indices,global_fixed_parameter_values,...
    test_data,max_neg_LL_val)
    % EP 17-11-07

    % Calculates likelihood and gradient of test_data from a gaussian, given mu and sigma parameters

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % identify parameter values

    % some parameters may be 'fixed' (i.e. not fitted by current iteration of
        % MLE); these are provided to the function in a separate list
    % compile a list of all parameters, fixed or fitted, and identify values
        % belonging to each individual parameter using that list
    param_vals = NaN(size(global_fixed_parameter_indices));
    param_vals(global_fixed_parameter_indices) = global_fixed_parameter_values(~isnan(global_fixed_parameter_values));
    param_vals(~global_fixed_parameter_indices) = param_vals_partial;
    
    mu = param_vals(1);
        % mean of distribution
    sigma = param_vals(2);
        % s.d. of distribution

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate likelihood of observing test_data given current global parameters
    tic;
    data_likelihoods = pdf('normal',test_data,mu,sigma);
    neg_combined_LL = -sum(log(data_likelihoods));

    d_LL_d_mu = d_LL_d_mu_norm_calc(test_data,mu,sigma);
    d_LL_d_sigma = d_LL_d_sigma_norm_calc(test_data,mu,sigma);

    global_gradient_vector = [-sum(d_LL_d_mu),-sum(d_LL_d_sigma)];
    
    global_gradient_vector_partial = global_gradient_vector(~global_fixed_parameter_indices);
    global_gradient_vector_partial(global_gradient_vector_partial>max_neg_LL_val) = max_neg_LL_val;
    global_gradient_vector_partial(global_gradient_vector_partial<-max_neg_LL_val) = -max_neg_LL_val;

    if neg_combined_LL > max_neg_LL_val
        neg_combined_LL = max_neg_LL_val;
    end
    
    runtime = toc;

end
    
