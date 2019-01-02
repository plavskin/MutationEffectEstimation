function MLE_finder(key_list, value_list)

    % EP 18-12-4

    % takes input from a data file, estimates best-fit parameters

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % To avoid overflow issues, set a max value that log likelihood can
        % be equal to throughout the MLE
    max_neg_LL_val = 10^50;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get input values
    input_value_dict = containers.Map(key_list,value_list);
    external_counter = str2num(input_value_dict('external_counter'));
    parameter_list = input_value_dict('parameter_list');
    output_file = input_value_dict('output_file');
    parallel_processors = input_value_dict('parallel_processors');
    ms_positions = input_value_dict('multistart_positions');
    ms_grid_parameter_array = input_value_dict('multistart_grid_parameters');
    global_mle_parameters = input_value_dict('top_level_parameters');
    tolx_val = input_value_dict('x_tolerance');
    tolfun_val = input_value_dict('fun_tolerance');
    pause_at_end = input_value_dict('pause_at_end');
    pre_MLE_function_name = input_value_dict('pre_MLE_function');
    post_MLE_function_name = input_value_dict('post_MLE_function');
    gradient_specification = input_value_dict('gradient_specification');
    model_code_location = input_value_dict('model_code_location');
    addpath(genpath(model_code_location));

    % add max_neg_LL_val to input_value_dict
    input_value_dict('max_neg_LL_val') = max_neg_LL_val;
    
    % process parameter name arrays into bool arrays
    indices_to_multistart = parameter_identifier(parameter_list,ms_grid_parameter_array);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    [global_param_number, global_mle_parameter_names, ...
        global_logspace_array, global_scaling_array, global_fixed_parameter_values, ...
        global_fixed_parameter_indices, global_lower_bounds_fitted, ...
        global_upper_bounds_fitted, global_start_vals_fitted] = ...
        parameter_array_subsetter(global_mle_parameters, input_value_dict);

    input_value_dict('global_mle_parameter_names') = global_mle_parameter_names;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % run pre_MLE_function, if it's supplied
    if ~isnan(pre_MLE_function_name)
        pre_MLE_function = str2func(pre_MLE_function_name);
        pre_MLE_output_dict = pre_MLE_function(input_value_dict);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % if any parameters are fitted, run MLE; otherwise, just run
        % LL_calculator directly
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global_start_time = tic;
    if sum(~global_fixed_parameter_indices) > 0
        % Set up global search structure

        gs = GlobalSearch;
        gs.TolFun = tolfun_val;
        gs.TolX = tolx_val;
        gs.NumStageOnePoints = 2^global_param_number; % default = 200; sufficient is 50
        gs.NumTrialPoints = 3^global_param_number; % default = 1000; sufficient is 200
        gs.StartPointsToRun='bounds';
        gs.Display='final';

        ms = MultiStart(gs);
            % assign same Tol, StartPointsToRun, and Display properties to ms as to gs
        if parallel_processors > 1
            ms.UseParallel = 1;
        else
            ms.UseParallel = 0;
        end

        fmincon_opts = optimoptions('fmincon','TolX',tolx_val,'TolFun',tolfun_val,...
            'Algorithm','interior-point','MaxIter',5000,'MaxFunEvals',12000,...
            'SpecifyObjectiveGradient',gradient_specification,'CheckGradients',...
            false,'Display','off');

        min_problem_fixed_params = createOptimProblem('fmincon','objective',...
            @(v) LL_calculator(v,...
                global_mle_parameter_names, global_fixed_parameter_indices,global_fixed_parameter_values,...
                global_logspace_array,global_scaling_array, max_neg_LL_val, input_value_dict, pre_MLE_output_dict),...
            'x0',global_start_vals_fitted,'lb',global_lower_bounds_fitted,'ub',global_upper_bounds_fitted,...
            'options',fmincon_opts);

        if ms_positions == 0
            % use globalsearch algorithm
            [vout_current,fval_current,exitflag_current,output_current,solutions_current]=...
                run(gs,min_problem_fixed_params);
        elseif ms_positions == 1
            % only use single startpoint, at global_start_values
            [vout_current,fval_current,exitflag_current]=...
                fmincon(min_problem_fixed_params);
        else

            if parallel_processors > 1
                parpool('local',parallel_processors);
            end

            %ms_positions_across_dimensions = ms_positions^global_param_number+1;
            indices_to_multistart_fitted = indices_to_multistart(~global_fixed_parameter_indices);

            ms_starting_point_mat = NaN;
            startpoints = MS_Startposition_Generator_v2(indices_to_multistart_fitted,ms_positions,...
                global_start_vals_fitted,global_lower_bounds_fitted,global_upper_bounds_fitted,...
                ms_starting_point_mat);
            [vout_current,fval_current,exitflag_current,output_current,solutions_current]=...
                run(ms,min_problem_fixed_params,startpoints);
        end

        vout_current_corrected = zeros(size(global_fixed_parameter_indices));
        vout_current_corrected(global_fixed_parameter_indices) = global_fixed_parameter_values(~isnan(global_fixed_parameter_values));
        vout_current_corrected(~global_fixed_parameter_indices) = vout_current;
        vout_current_corrected = reverse_value_scaler(vout_current_corrected,global_logspace_array,global_scaling_array);
    else
        
        [fval_current] = LL_calculator([],...
                global_fixed_parameter_indices,...
                global_fixed_parameter_values, global_logspace_array,...
                global_scaling_array, max_neg_LL_val, input_value_dict,...
                pre_MLE_output_dict);
        vout_current_corrected = reverse_value_scaler(global_fixed_parameter_values,...
            global_logspace_array,global_scaling_array);
    end
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

