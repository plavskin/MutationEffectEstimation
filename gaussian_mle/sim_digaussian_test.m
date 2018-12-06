function sim_digaussian_test(key_list, value_list)

    % takes input from a data file, and simulates the same number of datapoints

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get parameter values
    parameter_dict = containers.Map(key_list,value_list);
    
    external_counter = str2num(parameter_dict('external_counter'));
    combined_start_values_array_unscaled = parameter_dict('combined_start_values_array');
    parameter_list = parameter_dict('parameter_list');
    original_phenotype_file = parameter_dict('original_phenotype_file');
    phenotype_file = parameter_dict('phenotype_file');
    pause_at_end = parameter_dict('pause_at_end');

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
    lambda = combined_start_values_array_unscaled(1);
    mu_1 = combined_start_values_array_unscaled(2);
    sigma_1 = combined_start_values_array_unscaled(3);
    mu_2 = combined_start_values_array_unscaled(4);
    sigma_2 = combined_start_values_array_unscaled(5);

    n = max(size(test_data));
    n_1 = binornd(n, lambda);
    n_2 = n-n_1;

    sim_data_1 = normrnd(mu_1, sigma_1, [n_1 1]);
    sim_data_2 = normrnd(mu_2, sigma_2, [n_2 1]);
    sim_data = [sim_data_1; sim_data_2];

    data_table.data = sim_data;
    writetable(data_table,phenotype_file);

    runtime = toc;

    if pause_at_end & runtime < 120
        pausetime=120-runtime;
        pause(pausetime)
    end

end

