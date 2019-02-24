function output_dict = pre_MLE_nested(input_value_dict)
	% Generates pre_MLE output dictionary for MLE_finder function being
		% run within an LL_calculator function inside another
		% MLE_finder_function, which had itself already calculated a
		% pre_MLE_output_dict; that pre_MLE_output_dict can then be
		% added to the input_value_dict, and returned through this
		% function

    output_dict = input_value_dict('pre_MLE_output_dict');

end