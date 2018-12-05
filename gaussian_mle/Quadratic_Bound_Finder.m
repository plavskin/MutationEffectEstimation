function Quadratic_Bound_Finder(key_list, value_list)

	% EP 07-02-18

	% Takes list of parameter values, as well as cdf values based on
		% LRT corresponding to each of these parameter values, and fits
		% the parameter value corresponding to cdf_bound
	tic;
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get parameter values
    parameter_dict = containers.Map(key_list,value_list);

    cdf_bound = log(1-parameter_dict('cdf_bound'));
    	% cdf_bound = 1-p_value
    mle_param_val = parameter_dict('mle_param_val');
    parameter_vals = parameter_dict('parameter_values');
    cdf_vals = log(1-parameter_dict('cdf_vals'));
    output_file = parameter_dict('output_file');
    quad_fit_file = parameter_dict('fit_file');
    profile_side = parameter_dict('profile_side');
    pause_at_end = parameter_dict('pause_at_end');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% initialize starting coefficients and set constrains
	% coeffs(1) [a]: must be negative
	% coeffs(2) [-b/2a]: x_max > -b/2a or x_max < -b/2a, depending on
		% direction of CI bound
		% this is the axis of the parabola; this requirement means that
			% all parameter_vals are fitted on a single,
			% 'ascending' side of the parabola
	% coeffs(3) [(-b^2 + 4ac)/(4a)]: cdf_bound <= (-b^2 + 4ac)/(4a)
		% this is the y-value of the vertex of the parabola; the
			% lower bound of this requirement ensures that if cdf_bound
			% needs to be extrapolated from the current data, it will
			% return a real value
	direct_fit_coefficients = polyfit(parameter_vals,cdf_vals,2);
	starting_coefficients = NaN([1,3]);
	starting_coefficients(1) = min(direct_fit_coefficients(1),(mle_param_val-10^-50));
		% make sure a is negative, i.e. parabola 'faces' down
%	if min(parameter_vals)<mle_param_val
	if strcmp(profile_side,'upper')
		current_lb = [-Inf,-Inf,cdf_bound];
		current_ub = [0,min(parameter_vals),Inf];
		starting_coefficients(2) = max(-direct_fit_coefficients(2)/(2*starting_coefficients(1)),min(parameter_vals));
			% second coefficient is -b/2a
%	else
	elseif strcmp(profile_side,'lower')
		current_lb = [-Inf,max(parameter_vals),cdf_bound];
		current_ub = [0,Inf,Inf];
		starting_coefficients(2) = min(-direct_fit_coefficients(2)/(2*starting_coefficients(1)),max(parameter_vals));
	else
		disp('Profile side selected:')
		disp(profile_side)
		error('Invalid profile_side')
	end
	starting_coefficients(3) = max(cdf_bound, ...
		(-starting_coefficients(2)^2 + 4 * starting_coefficients(1) * direct_fit_coefficients(3)) / ...
		(4 * starting_coefficients(1)));

	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	% Perform global search for optimal parameters
	gs = GlobalSearch;
    gs.NumStageOnePoints = 20; % default = 200; sufficient is 50
    gs.NumTrialPoints = 100; % default = 1000; sufficient is 200
    gs.StartPointsToRun='bounds';
    gs.Display='final';

    fmincon_opts = optimoptions('fmincon', 'Algorithm','interior-point', 'Display','off');

    global_start_time = tic;

    quad_diff_squared = @(coefficients, x_vals, y_vals) ...
    	sum((Quadratic_Fit_Vertex_Limited(coefficients,x_vals,y_vals)).^2);

    min_problem_fixed_params = createOptimProblem('fmincon','objective',...
        @(coeffs) quad_diff_squared(coeffs, parameter_vals, cdf_vals),...
        'x0',starting_coefficients,'lb',current_lb,'ub',current_ub,...
        'options',fmincon_opts);

    optimal_coefficients = run(gs,min_problem_fixed_params);

	a = optimal_coefficients(1);
	b = -optimal_coefficients(2)*2*a;
	c = optimal_coefficients(3)+(b^2)/(4*a);
	quadratic_fit = [a,b,c];
    
%    x_plot_vals = linspace(min(parameter_vals),max(parameter_vals),50);
%    figure; plot(parameter_vals,cdf_vals,'ob'); hold on;
%    plot(x_plot_vals,polyval(quadratic_fit,x_plot_vals),'-r'); hold off;
    
	possible_bounds = roots(quadratic_fit+[0,0,-cdf_bound]);
		% here, cdf_bound is subtracted from the value for c to
			% find the x-intercepts of the parabola at height cdf_bound

	% find the point out of possible_bounds that is closest to mle_param_val:
		% this is the intercept on the correct side of the parabola

%	[~,bound_index] = max(abs(possible_bounds-mle_param_val));
%   parameter_bound = possible_bounds(bound_index);
	% look for parameter_bound on the correct side of the
        % parabola's axis
    if strcmp(profile_side,'upper')
        parameter_bound = ...
            possible_bounds(possible_bounds > optimal_coefficients(2));
    elseif strcmp(profile_side,'lower')
        parameter_bound = ...
            possible_bounds(possible_bounds < optimal_coefficients(2));
    end

	dlmwrite(output_file,parameter_bound,'delimiter',',','precision',9);
	dlmwrite(quad_fit_file,quadratic_fit,'delimiter',',','precision',9);

	runtime = toc;

    if pause_at_end & runtime < 120
        pausetime=120-runtime;
        pause(pausetime)
    end

end