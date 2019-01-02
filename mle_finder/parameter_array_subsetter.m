function [param_number, mle_parameter_names, logspace_array, ...
    scaling_array, fixed_parameter_values, fixed_parameter_indices, ...
    lower_bounds_fitted, upper_bounds_fitted, start_vals_fitted] = ...
    parameter_array_subsetter(mle_parameters, input_value_dict)

    combined_fixed_parameter_array = input_value_dict('tempfixed_parameter_bool');
    combined_min_array_unscaled = input_value_dict('min_parameter_vals');
    combined_max_array_unscaled = input_value_dict('max_parameter_vals');
    combined_length_array = input_value_dict('profile_point_num_list');
    combined_position_array = cellfun(@str2num,input_value_dict('combined_position_array'));
    combined_start_values_array_unscaled = input_value_dict('starting_parameter_vals');
    combined_scaling_array = input_value_dict('scaling_array');
    combined_profile_ub_array_unscaled = input_value_dict('profile_upper_limits');
    combined_profile_lb_array_unscaled = input_value_dict('profile_lower_limits');
    combined_logspace_parameters = input_value_dict('logspace_profile_parameters');
    parameter_list = input_value_dict('parameter_list');

    % process parameter name arrays into bool arrays
    combined_logspace_array = parameter_identifier(parameter_list,combined_logspace_parameters);

    % rescale parameters and convert to logspace as needed
    combined_min_array = value_rescaler(combined_min_array_unscaled,combined_logspace_array,combined_scaling_array);
    combined_max_array = value_rescaler(combined_max_array_unscaled,combined_logspace_array,combined_scaling_array);
    combined_start_values_array = value_rescaler(combined_start_values_array_unscaled,combined_logspace_array,combined_scaling_array);
    combined_profile_ub_array = value_rescaler(combined_profile_ub_array_unscaled,combined_logspace_array,combined_scaling_array);
    combined_profile_lb_array = value_rescaler(combined_profile_lb_array_unscaled,combined_logspace_array,combined_scaling_array);

    if length(combined_position_array)==1
        combined_position_array = repmat(combined_position_array(1),size(parameter_list));
    end

    param_bool_list = parameter_identifier(parameter_list, mle_parameters);
    param_number = sum(param_bool_list);
    if param_number == 0
        mle_parameter_names = parameter_list;
        param_bool_list = true(size(param_bool_list));
        param_number = sum(param_bool_list);
    else
        mle_parameter_names = parameter_list(param_bool_list);
    end

    fixed_parameter_array = ...
        combined_fixed_parameter_array(param_bool_list);
    min_array = ...
        combined_min_array(param_bool_list);
    max_array = ...
        combined_max_array(param_bool_list);
    length_array = ...
        combined_length_array(param_bool_list);
    position_array = ...
        combined_position_array(param_bool_list);
    start_values = ...
        combined_start_values_array(param_bool_list);

    profile_lb_array = ...
        combined_profile_lb_array(param_bool_list);
    profile_ub_array = ...
        combined_profile_ub_array(param_bool_list);

    logspace_array = ...
        combined_logspace_array(param_bool_list);
    scaling_array = ...
        combined_scaling_array(param_bool_list);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up a list of fixed effect values for 'global' parameters

    % Create an array indicating which fixed parameters need to be
        % created on a log scale, rather than a linear scale
    fixed_parameter_values = ...
        fixed_parameter_processor(fixed_parameter_array,profile_lb_array,...
            profile_ub_array,length_array,position_array);

    fixed_parameter_indices = logical(fixed_parameter_array);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up lower and upper bounds for parameters

    lower_bounds_fitted = min_array(~fixed_parameter_indices);
    upper_bounds_fitted = max_array(~fixed_parameter_indices);
    start_vals_fitted = start_values(~fixed_parameter_indices);
end
