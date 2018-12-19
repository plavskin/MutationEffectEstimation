function [combined_LL, unscaled_gradient_vector, grad_parameter_names] = LL_calculator_multi_exponential(param_vals,...
    input_value_dict, pre_MLE_output_dict)
    % EP 17-11-07

    % Calculates likelihood and gradient of test_data from a gaussian, given mu and sigma parameters

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % identify parameter values
    parameter_list = input_value_dict('parameter_list');
    parameter_dict = containers.Map(parameter_list,param_vals);
    
    rate = parameter_dict('rate');
    shape = parameter_dict('shape');
    % translate into gamma distribution parameters
    scale = 1/rate;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract test_data
    test_data = pre_MLE_output_dict('test_data');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate likelihood of observing test_data given mu and sigma
    data_likelihoods = pdf('gamma', test_data, shape, scale);
    combined_LL = sum(log(data_likelihoods));

    d_LL_d_shape = d_LL_d_shape_gamma_calc(test_data, shape, rate);
    d_LL_d_rate = d_LL_d_rate_gamma_calc(test_data, shape, rate);

    unscaled_gradient_vector = [sum(d_LL_d_shape),sum(d_LL_d_rate)];
    grad_parameter_names = {'shape','rate'};

end
    
