function results = run_for_parameters(f, params, varargin)
%RUN_FOR_PARAMETERS Evaluate the given function for different paramters.
%Syntax:
%   results = RUN_FOR_PARAMETERS(f, params[, ...]);
%Inputs:
%   f - Function handle whose must return a single value of the same type
%       each time it is called.
%   params - An n x k matrix or cell array where each row represents a
%            parameter vector. The given function handle should accept
%            k parameters.
%   ... - Options:
%           'Parallel' - Whether to evaluate the function in parallel.
%                        Default value is false.
%           'ParallelLevel' - Set to 0 to run each task for each parameter
%                             vector. Set to 1 to run each task for each
%                             repetition. Default value is 0.
%           'Repeat' - Specify how many times the function will be
%                      evaluated for each parameter vector. This is useful
%                      in running Monte Carlo simulations. Default value is
%                      one.
%           'ShowProgress' - Shows a progress bar. Default value is false.
%Output:
%   results - Depending on the return value of the function. Let r be the
%             number of repeats.
%             1. If f returns a scalar, the returned value will be an n x r
%                matrix;
%             2. If f returns a struct, the returned value will be an n x r
%                struct;
%             3. Otherwise, the returned value will be an n x r cell array.
use_parallel = false;
parallel_level = 0;
show_progress = false;
repeat = 1;
for ii = 1:2:nargin-2
    option_name = varargin{ii};
    option_value = varargin{ii+1};
    switch lower(option_name)
        case 'parallel'
            use_parallel = option_value;
        case 'parallellevel'
            parallel_level = option_value;
        case 'repeat'
            repeat = option_value;
        case 'showprogress'
            show_progress = option_value;
        otherwise
            error('Unknown option "%s".', option_name);
    end
end

results = [];
n = size(params, 1);
is_params_cell = iscell(params);
first_run = true;
if use_parallel
    pool = gcp();
    if show_progress
        if parallel_level == 0
            progressbar('reset', n);
        else
            progressbar('reset', n*repeat);
        end
    end
    first_run = true;
    if parallel_level == 0
        % one task per parameter
        % TODO: split into batches when n is large
        for nn = 1:n
            if is_params_cell
                cur_args = params(nn, :);
            else
                cur_args = num2cell(params(nn, :));
            end
            ff(nn) = parfeval(pool, @run_repetitions , 1, f, cur_args, repeat);
        end
        % fetch the results
        for nn = 1:n
            [idx, cur_results] = fetchNext(ff);
            % we need to store a batch of results here
            if first_run
                results = allocate_results(cur_results(1), n, repeat);
                first_run = false;
            end
            results(idx,:) = cur_results;
            if show_progress
                progressbar('advance');
            end
        end
    else
        % one task per repetition
        for nn = 1:n
            if is_params_cell
                cur_args = params(nn, :);
            else
                cur_args = num2cell(params(nn, :));
            end
            for rr = 1:repeat
                ff(rr) = parfeval(pool, f, 1, cur_args{:});
            end
            % we need to fetch the results for each repetition here
            for rr = 1:repeat
                [idx, cur_result] = fetchNext(ff);
                if first_run
                    results = allocate_results(cur_result, n, repeat);
                    first_run = false;
                end
                if isnumeric(results) || isstruct(results)
                    results(nn,idx) = cur_result;
                else
                    results{nn,idx} = cur_result;
                end
                if show_progress
                    progressbar('advance');
                end
            end
        end
    end
    if show_progress
        progressbar('end');
    end
else
    if show_progress
        progressbar('reset', n * repeat);
    end
    for nn = 1:n
        if is_params_cell
            cur_args = params(nn, :);
        else
            cur_args = num2cell(params(nn, :));
        end
        for rr = 1:repeat
            cur_result = f(cur_args{:});
            if first_run
                results = allocate_results(cur_result, n, repeat);
                first_run = false;
            end
            if isnumeric(results) || isstruct(results)
                results(nn,rr) = cur_result;
            else
                results{nn,rr} = cur_result;
            end
            if show_progress
                progressbar('advance');
            end
        end
    end
    if show_progress
        progressbar('end');
    end
end

end

function results = run_repetitions(f, args, repeat)
% A wrapper function to evaluate f multiple times.
first_run = true;
for rr = 1:repeat
    cur_result = f(args{:});
    if first_run
        results = allocate_results(cur_result, 1, repeat);
        first_run = false;
    end
    if isnumeric(results) || isstruct(results)
        results(rr) = cur_result;
    else
        results{rr} = cur_result;
    end
end
end

function results = allocate_results(sample_result, n, r)
% determine return type
if isnumeric(sample_result) && isscalar(sample_result)
    results = zeros(n, r);
elseif isstruct(sample_result)
    field_names = fieldnames(sample_result);
    args = reshape([field_names'; cell(1, length(field_names))], 1, []);
    results = repmat(struct(args{:}), n, r);
else
    results = cell(n, r);
end
end

