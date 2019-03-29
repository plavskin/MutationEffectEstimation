function [neg_combined_LL,global_gradient_vector_partial] = ...
    LL_calculator(param_vals_partial, global_mle_parameter_names, ...
        global_fixed_parameter_indices,global_fixed_parameter_values, ...
        global_logspace_array, global_scaling_array, max_neg_LL_val, ...
        input_value_dict, pre_MLE_output_dict)
    % EP 17-11-07

    % Rescales parameters for feeding into likelihood function, and inverts
        % calculated likelihood and gradients to allow minimization in
        % maximum likelihood search

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
    LL_calculator_name = input_value_dict('LL_calculator');
    current_LL_calculator = str2func(LL_calculator_name);
    gradient_specification = input_value_dict('gradient_specification');

    % Calculate likelihood of observing test_data given current global parameters
    
    if gradient_specification
        [combined_LL, unscaled_global_gradient_vector, unscaled_global_gradient_key] = ...
            current_LL_calculator(param_vals, input_value_dict, pre_MLE_output_dict);
        % reorder unscaled_global_gradient_vector to match order of parameters in global_mle_parameter_names
        [~, parameter_order] = ismember(global_mle_parameter_names, unscaled_global_gradient_key);
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
            current_LL_calculator(param_vals, input_value_dict, pre_MLE_output_dict);
        global_gradient_vector_partial = [];
    end

    % convert LL to negative LL
    neg_combined_LL = -combined_LL;

    if neg_combined_LL > max_neg_LL_val
        neg_combined_LL = max_neg_LL_val;
    end

    % if checkpointing, save checkpoint file with current unscaled
        % parameter vals, assuming combined_LL is higher than previous
        % saved LL
    write_checkpoint = input_value_dict('write_checkpoint');
    if write_checkpoint
        checkpoint_file = input_value_dict('checkpoint_file');
        % check what previous checkpoint_LL was, if it existed
        if exist(checkpoint_file, 'file') == 2
            checkpoint_table = readtable(checkpoint_file);
            if combined_LL < checkpoint_table.LL(1)
                write_checkpoint = false;
            end
        end
        if write_checkpoint
            start_time = input_value_dict('global_start_time');
            current_runtime = toc(start_time);
            previous_checkpoint_time = input_value_dict('checkpoint_time');
            total_runtime = current_runtime + previous_checkpoint_time;
            checkpoint_data = ...
                num2cell([combined_LL, total_runtime, param_vals]');
            checkpoint_table = table(checkpoint_data{:}, ...
                'VariableNames', ...
                ['LL', 'runtime_in_secs', global_mle_parameter_names]);
            writetable(checkpoint_table, checkpoint_file);
        end
    end
    
end
    
