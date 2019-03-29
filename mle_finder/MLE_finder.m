function MLE_finder(key_list, value_list)

    % EP 18-12-4

    % takes input from a data file, estimates best-fit parameters
    input_value_dict = containers.Map(key_list,value_list);

    post_MLE_function_name = input_value_dict('post_MLE_function');
    external_counter = str2num(input_value_dict('external_counter'));
    output_file = input_value_dict('output_file');
    pause_at_end = input_value_dict('pause_at_end');
    parallel_processors = input_value_dict('parallel_processors');
    checkpoint_file = input_value_dict('checkpoint_file');
    write_checkpoint = input_value_dict('write_checkpoint');
    ms_positions = input_value_dict('multistart_positions');
    global_mle_parameters = input_value_dict('top_level_parameters');
    parameter_list = input_value_dict('parameter_list');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % To avoid overflow issues, set a max value that log likelihood can
        % be equal to throughout the MLE
    max_neg_LL_val = 10^300;

    % add max_neg_LL_val to input_value_dict
    input_value_dict('max_neg_LL_val') = max_neg_LL_val;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if doing MLE_search with one starting point, check whether there's a
        % checkpoint file specified; if so, read starting parameter values
        % from that file
    if write_checkpoint & ms_positions == 1 & ...
        exist(checkpoint_file, 'file') == 2
        checkpoint_table = readtable(checkpoint_file);
        checkpoint_starting_param_vals = ...
            checkpoint_table{1, global_mle_parameters};
        original_starting_param_vals = ...
            input_value_dict('starting_parameter_vals');
        starting_param_dict = ...
            containers.Map(parameter_list, original_starting_param_vals);
        for param_counter = 1:length(global_mle_parameters)
            current_param = global_mle_parameters{param_counter};
            current_val = checkpoint_starting_param_vals(param_counter);
            starting_param_dict(current_param) = current_val;
        end
        input_value_dict('starting_parameter_vals') = ...
            cell2mat(values(starting_param_dict, parameter_list));
        % see how much time was elapsed up to previous checkpoint
        input_value_dict('checkpoint_time') = ...
            checkpoint_table.runtime_in_secs(1);
    else
        input_value_dict('checkpoint_time') = 0;
    end
    % Run MLE search
    global_start_time = tic;

    % save global_start_time in input_value_dict to use for checkpointing
    input_value_dict('global_start_time') = global_start_time;
    
    [vout_current_corrected, fval_current, ~, ~, solutions_current, ...
        global_mle_parameter_names, ~, pre_MLE_output_dict] = ...
        MLE_search_executer(input_value_dict);

    runtime = toc(global_start_time);
    
    export_data = [-fval_current,runtime,vout_current_corrected];

    % if doing mutlistart or global search and one solution didn't converge, replace export_data with NA
    if exist('solution_counter','var') == 1
        for solution_counter = 1:length(solutions_current)
            temp_exit = solutions_current(solution_counter).Exitflag;
            if temp_exit < 1
                export_data = NaN(size(export_data));
            end
        end
    end
    
    %dlmwrite(output_file,export_data,'delimiter',',','precision',9);

    % modify runtime to include runtime before checkpoint
    runtime_incl_precheckpoint = ...
        runtime + input_value_dict('checkpoint_time');
    table_data = num2cell([-fval_current,runtime_incl_precheckpoint,vout_current_corrected,external_counter]);
    table_headers = {'LL','runtime_in_secs',global_mle_parameter_names{:},'point_num'};
    T = table(table_data{:},'VariableNames',table_headers);
    writetable(T,output_file);

    if ~isnan(post_MLE_function_name)
        post_MLE_function = str2func(post_MLE_function_name);
        post_MLE_function(input_value_dict, pre_MLE_output_dict, T);
    end

    % close current parallel pool
    if parallel_processors > 1
        delete(gcp('nocreate'));
    end

    
    if pause_at_end & runtime < 120
        pausetime=120-runtime;
        pause(pausetime)
    end

end

