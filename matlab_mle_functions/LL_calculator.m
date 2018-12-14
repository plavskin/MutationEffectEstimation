function [neg_combined_LL,global_gradient_vector_partial] = LL_calculator(param_vals_partial,...
    global_fixed_parameter_indices,global_fixed_parameter_values,...
    global_logspace_array, global_scaling_array, max_neg_LL_val, input_value_dict, pre_MLE_output_dict)
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
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get name of function to calculate LL and convert it to function handle
    LL_calculator_name = input_value_dict('LL_calculator');
    current_LL_calculator = str2func(LL_calculator_name);
    gradient_specification = input_value_dict('gradient_specification');

    % Calculate likelihood of observing test_data given current global parameters
    tic;
    
    if gradient_specification
        [combined_LL, unscaled_global_gradient_vector, unscaled_global_gradient_key] = ...
            current_LL_calculator(param_vals, input_value_dict, pre_MLE_output_dict);
        % reorder unscaled_global_gradient_vector to match order of parameters in parameter_list
        parameter_list = input_value_dict('parameter_list');
        [~, parameter_order] = ismember(unscaled_global_gradient_key, parameter_list);
        unscaled_global_gradient_vector_ordered = unscaled_global_gradient_vector(parameter_order);
        % convert gradient vector to negative and rescale
        neg_unscaled_global_gradient_vector = -unscaled_global_gradient_vector_ordered;
        global_gradient_vector = gradient_value_rescaler(neg_unscaled_global_gradient_vector,...
            param_vals,global_logspace_array,global_scaling_array);
        % return only those values in gradient vector that correspond to fitted parameters
        global_gradient_vector_partial = global_gradient_vector(~global_fixed_parameter_indices);
        % deal with overflow
        global_gradient_vector_partial(global_gradient_vector_partial>max_neg_LL_val) = max_neg_LL_val;
        global_gradient_vector_partial(global_gradient_vector_partial<-max_neg_LL_val) = -max_neg_LL_val;
    else
        [combined_LL] = ...
            current_LL_calculator(param_vals, input_value_dict, pre_MLE_output_dict)
        global_gradient_vector_partial = []
    end

    % convert LL to negative LL
    neg_combined_LL = -combined_LL;

    if neg_combined_LL > max_neg_LL_val
        neg_combined_LL = max_neg_LL_val;
    end
    
    runtime = toc;

end
    
