function sim_gaussian_test(key_list, value_list)

    % takes input from a data file, and simulates the same number of datapoints

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get parameter values
    parameter_dict = containers.Map(key_list,value_list);
    
    external_counter = str2num(parameter_dict('external_counter'));
    combined_start_values_array_unscaled = parameter_dict('combined_start_values_array');
    parameter_list = parameter_dict('parameter_list');
    original_phenotype_file = parameter_dict('original_phenotype_file');
    phenotype_file = parameter_dict('phenotype_file');


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Use a random seed that's the sum of the current time in seconds and
        % external_counter, so that mutliple jobs starting at the same time have
        % different seeds
    rng_shuffle = rng('shuffle');
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
    mu = combined_start_values_array_unscaled(1);
    sigma = combined_start_values_array_unscaled(2);
    
    sim_data = normrnd(mu, sigma, size(test_data));

    data_table.data = sim_data;
    writetable(data_table,phenotype_file);

    runtime = toc;

    pause_at_end = true;

    if pause_at_end & runtime < 120
        pausetime=120-runtime;
        pause(pausetime)
    end

end

