function output_vals = reverse_value_scaler(input_vals, logspace_array, scaling_array)

    % EP 18-08-08
    % Descales input_vals and then converts those that were in logspace back to linear space

    descaled_outputs = input_vals./scaling_array;
    output_vals = descaled_outputs;
    output_vals(logspace_array) = exp(descaled_outputs(logspace_array));

end