function [neg_combined_LL,global_gradient_vector_partial] = LL_calculator_gaussian_test(param_vals_partial,...
    global_fixed_parameter_indices,global_fixed_parameter_values,...
    global_logspace_array, global_scaling_array, max_neg_LL_val, parameter_dict, pre_MLE_output_dict)
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
    param_vals = reverse_value_scaler(param_vals,global_logspace_array,global_scaling_array);
    
    lambda = param_vals(1);
        % proportion of total distribution coming from distribution 1
    mu_1 = param_vals(2);
        % mean of distribution 1
    sigma_1 = param_vals(3);
        % s.d. of distribution 1
    mu_2 = param_vals(4);
        % mean of distribution 2
    sigma_2 = param_vals(5);
        % s.d. of distribution 2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract test_data
    test_data = pre_MLE_output_dict('test_data');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate likelihood of observing test_data given current global parameters
    tic;
    data_likelihoods = lambda*pdf('normal',test_data,mu_1,sigma_1)+...
        (1-lambda)*pdf('normal',test_data,mu_2,sigma_2);
    neg_combined_LL = -sum(log(data_likelihoods));

    d_LL_d_mu_1 = d_LL_d_mu_dinorm_calc(test_data,lambda,mu_1,sigma_1,mu_2,sigma_2);
    d_LL_d_sigma_1 = d_LL_d_sigma_dinorm_calc(test_data,lambda,mu_1,sigma_1,mu_2,sigma_2);
    d_LL_d_mu_2 = d_LL_d_mu_dinorm_calc(test_data,(1-lambda),mu_2,sigma_2,mu_1,sigma_1);
    d_LL_d_sigma_2 = d_LL_d_sigma_dinorm_calc(test_data,(1-lambda),mu_2,sigma_2,mu_1,sigma_1);
    d_LL_d_lambda = d_LL_d_lambda_dinorm_calc(test_data,lambda,mu_1,sigma_1,mu_2,sigma_2);

    unscaled_global_gradient_vector = [-sum(d_LL_d_lambda),-sum(d_LL_d_mu_1),...
        -sum(d_LL_d_sigma_1),-sum(d_LL_d_mu_2),-sum(d_LL_d_sigma_2)];
    global_gradient_vector = gradient_value_rescaler(unscaled_global_gradient_vector,...
        param_vals,global_logspace_array,global_scaling_array);
    
    global_gradient_vector_partial = global_gradient_vector(~global_fixed_parameter_indices);
    global_gradient_vector_partial(global_gradient_vector_partial>max_neg_LL_val) = max_neg_LL_val;
    global_gradient_vector_partial(global_gradient_vector_partial<-max_neg_LL_val) = -max_neg_LL_val;

    if neg_combined_LL > max_neg_LL_val
        neg_combined_LL = max_neg_LL_val;
    end
    
    runtime = toc;

end
    
