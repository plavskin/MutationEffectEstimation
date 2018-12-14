function output_dict = pre_MLE_gaussian_test(input_value_dict)
    % EP 18-12-04

    % Reads in data and outputs test_data
    phenotype_file = input_value_dict('phenotype_file');

    data_table = readtable(phenotype_file);
    test_data = data_table.data;

    output_dict = containers.Map(['test_data'],[test_data]);