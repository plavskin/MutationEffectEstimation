function rescaled_gradient = gradient_value_rescaler(gradient_vals, parameter_vals, logspace_array, scaling_array)

    % EP 18-08-08
    % Gradients are computed in the space of x, but MLE runs in the space of scaled (and potentially log-transformed) x
    % This function rescales the gradient vector to the space relevant to the MLE
    rescaling_vector = 1./scaling_array;
    delog_vector = ones(size(logspace_array));
    delog_vector(logspace_array) = parameter_vals(logspace_array);
    rescaling_delog_vector = rescaling_vector.*delog_vector;
    rescaled_gradient = gradient_vals.*rescaling_delog_vector;

end