function output_vals = value_rescaler(input_vals, logspace_array, scaling_array)

    % EP 18-08-07
    % Converts input_vals for which logspace_array is true to the ln() of those values
    % Multiplies each value (after logspace conversion where necessay) by scaling array

    logspace_converted_vals = input_vals;
    logspace_converted_vals(logspace_array) = log(input_vals(logspace_array));
    output_vals = logspace_converted_vals.*scaling_array;

end