function output_dict = pre_MLE_gaussian_test(key_list, value_list)
    % EP 18-12-04

    % Reads in data and outputs test_data
    parameter_dict = containers.Map(key_list,value_list);
    phenotype_file = parameter_dict('phenotype_file');

    data_table = readtable(phenotype_file);
    test_data = data_table.data;

    output_dict = containers.Map(['test_data'],[test_data]);