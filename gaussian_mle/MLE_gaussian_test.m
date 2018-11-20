function MLE_gaussian_test(key_list, value_list)

    % EP 18-4-4

    % takes input from a data file, estimates best-fit gaussian parameters

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % To avoid overflow issues, set a max value that log likelihood can
        % be equal to throughout the MLE
    %max_neg_LL_val = realmax/(strain_number*2);
    max_neg_LL_val = 10^50;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get parameter values
    parameter_dict = containers.Map(key_list,value_list);

    external_counter = str2num(parameter_dict('external_counter'));
    combined_fixed_parameter_array = parameter_dict('combined_fixed_parameter_array');
    combined_min_array_unscaled = parameter_dict('combined_min_array');
    combined_max_array_unscaled = parameter_dict('combined_max_array');
    combined_length_array = parameter_dict('combined_length_array');
    combined_position_array = cellfun(@str2num,parameter_dict('combined_position_array'));
    combined_start_values_array_unscaled = parameter_dict('combined_start_values_array');
    combined_scaling_array = parameter_dict('combined_scaling_array');
    parameter_list = parameter_dict('parameter_list');
    output_file = parameter_dict('output_file');
    parallel_processors = parameter_dict('parallel_processors');
    ms_positions = parameter_dict('ms_positions');
    combined_profile_ub_array_unscaled = parameter_dict('combined_profile_ub_array');
    combined_profile_lb_array_unscaled = parameter_dict('combined_profile_lb_array');
    ms_grid_parameter_array = parameter_dict('ms_grid_parameter_array');
    combined_logspace_parameters = parameter_dict('combined_logspace_parameters');
    tolx_val = parameter_dict('tolx_val');
    tolfun_val = parameter_dict('tolfun_val');
    % optional parameters
    phenotype_file = parameter_dict('phenotype_file');
    % process parameter name arrays into bool arrays
    combined_logspace_array = parameter_identifier(parameter_list,combined_logspace_parameters);
    indices_to_multistart = parameter_identifier(parameter_list,ms_grid_parameter_array);

    % rescale parameters and convert to logspace as needed
    combined_min_array = value_rescaler(combined_min_array_unscaled,combined_logspace_array,combined_scaling_array);
    combined_max_array = value_rescaler(combined_max_array_unscaled,combined_logspace_array,combined_scaling_array);
    combined_start_values_array = value_rescaler(combined_start_values_array_unscaled,combined_logspace_array,combined_scaling_array);
    combined_profile_ub_array = value_rescaler(combined_profile_ub_array_unscaled,combined_logspace_array,combined_scaling_array);
    combined_profile_lb_array = value_rescaler(combined_profile_lb_array_unscaled,combined_logspace_array,combined_scaling_array);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get data
    disp(phenotype_file)
    data_table = readtable(phenotype_file);

    test_data = data_table.data;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    global_param_number = 2;

    if length(combined_position_array)==1
        combined_position_array = repmat(combined_position_array(1),size(parameter_list));
    end

    global_fixed_parameter_array = ...
        combined_fixed_parameter_array(1:global_param_number);
    global_min_array = ...
        combined_min_array(1:global_param_number);
    global_max_array = ...
        combined_max_array(1:global_param_number);
    global_length_array = ...
        combined_length_array(1:global_param_number);
    global_position_array = ...
        combined_position_array(1:global_param_number);
    global_start_values = ...
        combined_start_values_array(1:global_param_number);

    global_profile_lb_array = ...
        combined_profile_lb_array(1:global_param_number);
    global_profile_ub_array = ...
        combined_profile_ub_array(1:global_param_number);

    global_logspace_array = ...
        combined_logspace_array(1:global_param_number);
    global_scaling_array = ...
        combined_scaling_array(1:global_param_number);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up a list of fixed effect values for 'global' parameters

    % Create an array indicating which fixed parameters need to be
        % created on a log scale, rather than a linear scale
    global_fixed_parameter_values = ...
        fixed_parameter_processor(global_fixed_parameter_array,global_profile_lb_array,...
            global_profile_ub_array,global_length_array,global_position_array);

    global_fixed_parameter_indices = logical(global_fixed_parameter_array);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set up lower and upper bounds for parameters

    global_lower_bounds_fitted = global_min_array(~global_fixed_parameter_indices);
    global_upper_bounds_fitted = global_max_array(~global_fixed_parameter_indices);
    global_start_vals_fitted = global_start_values(~global_fixed_parameter_indices);

    global_param_number_fitted = length(global_lower_bounds_fitted);

    global_profile_lb_fitted = global_profile_lb_array(~global_fixed_parameter_indices);
    global_profile_ub_fitted = global_profile_ub_array(~global_fixed_parameter_indices);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
        'SpecifyObjectiveGradient',true,'CheckGradients',false,'Display','off');

    % To find parameters corresponding to each strain, use simple
        % space search with one start point or global search with
        % multiple starts?
    %strainwise_search_type = 'local';   
    strainwise_search_type = 'global';    

    global_start_time = tic;

    min_problem_fixed_params = createOptimProblem('fmincon','objective',...
        @(v) LL_calculator_gaussian_test(v,...
            global_fixed_parameter_indices,global_fixed_parameter_values,...
            test_data,max_neg_LL_val,global_logspace_array,global_scaling_array),...
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

        startpoints = MS_Startposition_Generator_v2(indices_to_multistart_fitted,ms_positions,...
            global_start_vals_fitted,global_lower_bounds_fitted,global_upper_bounds_fitted);
        [vout_current,fval_current,exitflag_current,output_current,solutions_current]=...
            run(ms,min_problem_fixed_params,startpoints);
    end
    
    vout_current_corrected = zeros(size(global_fixed_parameter_indices));
    vout_current_corrected(global_fixed_parameter_indices) = global_fixed_parameter_values(~isnan(global_fixed_parameter_values));
    vout_current_corrected(~global_fixed_parameter_indices) = vout_current;
    vout_current_corrected = reverse_value_scaler(vout_current_corrected,global_logspace_array,global_scaling_array);

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
    table_headers = {'LL','runtime_in_secs',parameter_list{:},'point_num'};
    T = table(table_data{:},'VariableNames',table_headers);
    writetable(T,output_file);

%    if runtime < 120
%        pausetime=120-runtime;
%        pause(pausetime)
%    end

end

