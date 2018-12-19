function sim_gaussian_test(key_list, value_list)

    % takes input from a data file, and simulates the same number of datapoints

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get parameter values
    input_value_dict = containers.Map(key_list,value_list);
    
    external_counter = str2num(input_value_dict('external_counter'));
    combined_start_values_array_unscaled = input_value_dict('combined_start_values_array');
    parameter_list = input_value_dict('parameter_list');
    original_phenotype_file = input_value_dict('original_phenotype_file');
    phenotype_file = input_value_dict('phenotype_file');
    pause_at_end = input_value_dict('pause_at_end');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Use a random seed that's the sum of the current time in seconds and
        % external_counter, so that mutliple jobs starting at the same time have
        % different seeds
    rng('shuffle')
    rng_shuffle = rng;
    random_seed = rng_shuffle.Seed + external_counter;
        % note that if random_seed exceeds 2^32, it maxes out
    if random_seed > 2^32
        random_seed = rng_shuffle.Seed - external_counter;
    end
    rng(random_seed);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get data
    tic;
    data_table = readtable(original_phenotype_file);
    test_data = data_table.data;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parameter_dict = containers.Map(parameter_list,...
        combined_start_values_array_unscaled);
    
    rate = parameter_dict('rate');
    mu = 1/rate;
    
    sim_data = exprnd(mu, size(test_data));

    data_table.data = sim_data;
    writetable(data_table,phenotype_file);

    runtime = toc;

    if pause_at_end & runtime < 120
        pausetime=120-runtime;
        pause(pausetime)
    end

end

