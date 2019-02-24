function MLE_finder(key_list, value_list)

    % EP 18-12-4

    % takes input from a data file, estimates best-fit parameters
    input_value_dict = containers.Map(key_list,value_list);

    post_MLE_function_name = input_value_dict('post_MLE_function');
    external_counter = str2num(input_value_dict('external_counter'));
    output_file = input_value_dict('output_file');
    pause_at_end = input_value_dict('pause_at_end');
    parallel_processors = input_value_dict('parallel_processors');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % To avoid overflow issues, set a max value that log likelihood can
        % be equal to throughout the MLE
    max_neg_LL_val = 10^50;

    % add max_neg_LL_val to input_value_dict
    input_value_dict('max_neg_LL_val') = max_neg_LL_val;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run MLE search
    global_start_time = tic;

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

    table_data = num2cell([-fval_current,runtime,vout_current_corrected,external_counter]);
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

