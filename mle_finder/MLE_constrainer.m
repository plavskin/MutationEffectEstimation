function [c, ceq, c_gradient, ceq_gradient] = MLE_constrainer(param_vals_partial,...
    global_mle_parameter_names, global_fixed_parameter_indices,global_fixed_parameter_values,...
    global_logspace_array, global_scaling_array, max_neg_LL_val, input_value_dict, pre_MLE_output_dict)

    % Rescales parameters for feeding into nonlinear constraint function,
        % and inverts calculated likelihood and gradients to allow
        % minimization in maximum likelihood search

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
    
    % pass a list of fitted parameters to LL calculation code in case they
        % are used to determine which gradients to calculate
    fitted_parameters = ...
        global_mle_parameter_names(~global_fixed_parameter_indices);
    pre_MLE_output_dict('fitted_parameters') = fitted_parameters;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get name of function to calculate LL and convert it to function handle
    nonlinear_constraint_function = input_value_dict('nonlinear_constraint_function');
    current_MLE_constrainer = str2func(nonlinear_constraint_function);
    gradient_specification = input_value_dict('gradient_specification');

    % Calculate likelihood of observing test_data given current global parameters
    
    if gradient_specification
        [c, ceq, unscaled_c_gradient_vals, unscaled_ceq_gradient_vals, ...
            unscaled_gradient_key] = ...
            current_MLE_constrainer(param_vals, input_value_dict, ...
                pre_MLE_output_dict);
        % reorder unscaled_c_gradient_vals and unscaled_ceq_gradient_vals
            % to match order of parameters in global_mle_parameter_names
        [~, parameter_order] = ...
            ismember(unscaled_gradient_key, global_mle_parameter_names);
%            ismember(global_mle_parameter_names, unscaled_gradient_key);
        
        parameter_order = parameter_order(parameter_order > 0);
        unscaled_grad_array = ...
            {unscaled_c_gradient_vals, unscaled_ceq_gradient_vals};
        scaled_partial_grad_array = {[], []};
        for counter = 1:2
            current_unscaled_grad_vals = unscaled_grad_array{counter};
            if ~isempty(current_unscaled_grad_vals)
                current_unscaled_grad_ordered = ...
                    zeros(size(global_mle_parameter_names));
                current_unscaled_grad_ordered(parameter_order) = ...
                    current_unscaled_grad_vals;
                % rescale gradients
                current_scaled_grad = ...
                    gradient_value_rescaler(current_unscaled_grad_ordered,...
                        param_vals, global_logspace_array, ...
                        global_scaling_array);
                % return only those values in gradient vector that
                    % correspond to fitted parameters
                current_grad_partial = ...
                    current_scaled_grad(~global_fixed_parameter_indices);
                % deal with overflow
                current_grad_partial(current_grad_partial > ...
                    max_neg_LL_val) = max_neg_LL_val;
                current_grad_partial(current_grad_partial < ...
                    -max_neg_LL_val) = -max_neg_LL_val;
                scaled_partial_grad_array{counter} = current_grad_partial;
            end
        end
        c_gradient = scaled_partial_grad_array{1}';
        ceq_gradient = scaled_partial_grad_array{2}';
    else
        [c, ceq] = ...
            current_MLE_constrainer(param_vals, input_value_dict, ...
                pre_MLE_output_dict);
        c_gradient = [];
        ceq_gradient = [];
    end
    
end
    
