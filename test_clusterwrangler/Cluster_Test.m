function Cluster_Test(key_list, value_list)

    % Writes text to output_file

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get input values
    input_value_dict = containers.Map(key_list,value_list);

    external_counter = str2num(input_value_dict('external_counter'));
    output_file = input_value_dict('output_file');
    code_options = input_value_dict('code_options');

    external_counter = repmat(external_counter, size(code_options'));
    code_options = code_options';
    
    data_table = table(external_counter, code_options);
    writetable(data_table, output_file);