function [combined_LL, unscaled_gradient_vector, grad_parameter_names] = LL_calculator_gaussian_test(param_vals,...
    input_value_dict, pre_MLE_output_dict)
    % EP 17-11-07

    % Calculates likelihood and gradient of test_data from a gaussian, given mu and sigma parameters

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % identify parameter values
    parameter_list = input_value_dict('parameter_list');
    parameter_dict = containers.Map(parameter_list,param_vals);
    
    mu = parameter_dict('mu');
        % mean of distribution
    sigma = parameter_dict('sigma');
        % s.d. of distribution

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract test_data
    test_data = pre_MLE_output_dict('test_data');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate likelihood of observing test_data given mu and sigma
    data_likelihoods = pdf('normal',test_data,mu,sigma);
    combined_LL = sum(log(data_likelihoods));

    d_LL_d_mu = d_LL_d_mu_norm_calc(test_data,mu,sigma);
    d_LL_d_sigma = d_LL_d_sigma_norm_calc(test_data,mu,sigma);

    unscaled_gradient_vector = [sum(d_LL_d_mu),sum(d_LL_d_sigma)];
    grad_parameter_names = {'mu','sigma'};

end
    
