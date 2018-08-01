function [fixed_parameter_values] = fixed_parameter_processor(fixed_parameter_array,min_array,max_array,...
	length_array,position_array,logspace_array)


	fixed_parameter_indices = false(size(fixed_parameter_array));
    fixed_parameter_values = NaN(size(fixed_parameter_array));

    for parameter_counter = 1:length(fixed_parameter_array)

%        current_fixed_parameter = fixed_parameter_array(parameter_counter);

        if ~fixed_parameter_array(parameter_counter)
            fixed_parameter_values(parameter_counter) = NaN;
            %current_param_val = NaN;
        else
            current_parameter_min = min_array(parameter_counter);
            current_parameter_max = max_array(parameter_counter);
            current_parameter_array_length = length_array(parameter_counter);
            current_parameter_position = position_array(parameter_counter);
        
%            current_parameter_index = strcmp(parameter_names,current_fixed_parameter);
%            fixed_parameter_indices = fixed_parameter_indices | current_parameter_index;
                % updates fixed_parameter_indices to include true at
                    % position of current param
            if logspace_array(parameter_counter)
            	current_param_val_array = logspace(log10(current_parameter_min),log10(current_parameter_max),current_parameter_array_length);
            else
            	current_param_val_array = linspace(current_parameter_min,current_parameter_max,current_parameter_array_length);
            end

%            % use linear space if fitted parameter is q or p0; otherwise, use log space
%            if (strcmp('q_1',current_fixed_parameter) || strcmp('p0_1',current_fixed_parameter))
%                current_param_val_array = linspace(current_parameter_min*linear_multiplier,current_parameter_max*linear_multiplier,current_parameter_array_length);
%            elseif (strcmp('combined_mean_2',current_fixed_parameter) || strcmp('combined_sd_2',current_fixed_parameter))
%                current_param_val_array = linspace(current_parameter_min,current_parameter_max,current_parameter_array_length);
%            else
%                current_param_val_array = log(logspace(log10(current_parameter_min),log10(current_parameter_max),current_parameter_array_length));
%            end

            current_param_val = current_param_val_array(current_parameter_position);
            fixed_parameter_values(parameter_counter) = current_param_val;
        end

    end

end