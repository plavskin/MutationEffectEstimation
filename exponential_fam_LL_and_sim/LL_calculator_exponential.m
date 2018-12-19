function [combined_LL, unscaled_gradient_vector, grad_parameter_names] = LL_calculator_exponential(param_vals,...
    input_value_dict, pre_MLE_output_dict)
    % EP 17-11-07

    % Calculates likelihood and gradient of test_data from a gaussian, given mu and sigma parameters

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % identify parameter values
    parameter_list = input_value_dict('parameter_list');
    parameter_dict = containers.Map(parameter_list,param_vals);
    
    rate = parameter_dict('rate');
    mu = 1/rate;
        % mean of distribution

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract test_data
    test_data = pre_MLE_output_dict('test_data');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate likelihood of observing test_data given mu and sigma
    data_likelihoods = pdf('exponential',test_data,mu);
    combined_LL = sum(log(data_likelihoods));

    d_LL_d_rate = d_LL_d_rate_exp_calc(test_data, rate);

    unscaled_gradient_vector = [sum(d_LL_d_rate)];
    grad_parameter_names = {'rate'};

end
    
