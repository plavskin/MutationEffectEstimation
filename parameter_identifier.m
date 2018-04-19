function parameter_bool_list =  parameter_identifier(parameter_list,parameter_sublist)
	% 18-04-04
	% Identifies indices in parameter_list that correspond to members of parameter_sublist

	[~,sublist_indices] = intersect(parameter_list,parameter_sublist);
    parameter_bool_list = false(size(parameter_list));
    parameter_bool_list(sublist_indices) = true;
