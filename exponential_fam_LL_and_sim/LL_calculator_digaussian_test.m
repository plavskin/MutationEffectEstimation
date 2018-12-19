function [combined_LL, unscaled_gradient_vector, grad_parameter_names] = LL_calculator_gaussian_test(param_vals,...
    input_value_dict, pre_MLE_output_dict)
    % EP 17-11-07

    % Calculates likelihood and gradient of test_data from a gaussian, given mu and sigma parameters

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % identify parameter values
    parameter_list = input_value_dict('parameter_list');
    parameter_dict = containers.Map(parameter_list,param_vals);
    
    lambda = parameter_dict('lambda');
        % proportion of total distribution coming from distribution 1
    mu_1 = parameter_dict('mu');
        % mean of distribution 1
    sigma_1 = parameter_dict('sigma');
        % s.d. of distribution 1
    rel_mu2 = parameter_dict('rel_mu2');
    mu_2 = mu_1 + rel_mu2;
        % mean of distribution 2
    sigma_2 = parameter_dict('sigma2');
        % s.d. of distribution 2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % extract test_data
    test_data = pre_MLE_output_dict('test_data');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate likelihood of observing test_data given current parameters

    data_likelihoods = lambda*pdf('normal',test_data,mu_1,sigma_1)+...
        (1-lambda)*pdf('normal',test_data,mu_2,sigma_2);
    combined_LL = sum(log(data_likelihoods));

    d_LL_d_sigma_1 = d_LL_d_sigma_dinorm_calc(test_data,lambda,mu_1,sigma_1,mu_2,sigma_2);
    d_LL_d_mu_2 = d_LL_d_mu_dinorm_calc(test_data,(1-lambda),mu_2,sigma_2,mu_1,sigma_1);
    d_LL_d_sigma_2 = d_LL_d_sigma_dinorm_calc(test_data,(1-lambda),mu_2,sigma_2,mu_1,sigma_1);
    d_LL_d_lambda = d_LL_d_lambda_dinorm_calc(test_data,lambda,mu_1,sigma_1,mu_2,sigma_2);

    d_LL_d_rel_mu2 = d_LL_d_mu_2;
    d_LL_d_mu_1 = d_LL_d_mu_dinorm_calc(test_data,lambda,mu_1,sigma_1,mu_2,sigma_2) + d_LL_d_mu_2;

    unscaled_gradient_vector = [sum(d_LL_d_lambda),sum(d_LL_d_mu_1),...
        sum(d_LL_d_sigma_1),sum(d_LL_d_rel_mu2),sum(d_LL_d_sigma_2)];
    grad_parameter_names = {'lambda','mu','sigma','rel_mu2','sigma2'};

end
    
