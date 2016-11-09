function progressbar(option, varargin)
%PROGRESSBAR Prints a progress bar.
%Syntax:
%   PROGRESSBAR('reset', total_progress)
%   PROGRESSBAR('resettimer')
%   PROGRESSBAR('advance'[, step])
%   PROGRESSBAR('barwidth', bar_width)
%   PROGRESSBAR('displaymode', mode)
%   PROGRESSBAR('minimalupdateinterval', min_interval);

persistent total_progress;
persistent current_progress;
persistent last_text_width;
persistent start_time;
persistent total_time;
persistent bar_width;
persistent display_mode;
persistent min_interval;
persistent last_print_time;

if isempty(total_progress)
    total_progress = 0;
    current_progress = 0;
    last_text_width = 0;
    bar_width = 20;
    display_mode = 0;
    min_interval = 0.3;
end

switch lower(option)
    case 'minimalupdateinterval'
        min_interval = max(double(varargin{1}), 0);
    case 'displaymode'
        if nargin == 2
            switch lower(varargin{1})
                case 'replace'
                    display_mode = 0;
                case 'append'
                    display_mode = 1;
                otherwise
                    error('Display mode should be either "replace" (default) or "append".');
            end
        else
            error('Display mode not specified.');
        end
    case 'barwidth'
        if nargin == 2
            if ischar(varargin{1}) && strcmpi(varargin{1}, 'default')
                bar_width = 20;
            elseif isreal(varargin{1}) && varargin{1} > 0 && floor(varargin{1}) == varargin{1}
                bar_width = varargin{1};
            else
                error('Bar width must be a positive integer or "default".');
            end
        else
            error('Bar width is not specified.');
        end
    case 'resettimer'
        start_time = datetime();
        total_time = [];
    case 'reset'
        start_time = datetime();
        total_time = [];
        current_progress = 0;
        last_text_width = 0;
        if nargin == 2
            total_progress = varargin{1};
            if total_progress <= 0 || floor(total_progress) ~= total_progress
                error('Positive integer expected for the total progress.');
            end
        end
    case 'advance'
        if nargin >= 2
            step = varargin{1};
            if isscalar(step) && isreal(step) && floor(step) == step
                current_progress = current_progress + step;
            else
                error('Step size should be an integer.');
            end
        else
            current_progress = current_progress + 1;
        end
        if current_progress > total_progress
            current_progress = total_progress;
        end
        now_time = datetime();
        % check if this function is called too frequently
        if current_progress < total_progress && ...
                ~isempty(last_print_time) && ...
                seconds(now_time - last_print_time) < min_interval
            return;
        end
        last_print_time = now_time;
        % print progress
        percentage = current_progress / total_progress;
        if percentage == 1.0        
            % completed
            if isempty(total_time)
                total_time = now_time - start_time;
            end
            text2print = sprintf('%s (%d/%d) %.1f%%%% [%s]', ...
                get_progress_bar(percentage, bar_width), ...
                current_progress, total_progress, ...
                percentage * 100, seconds2str(seconds(total_time)));
        else
            % incomplete, generate ETA string
            if last_text_width == 0
                eta_str = '-';
            else
                seconds_elapsed = seconds(now_time - start_time);
                seconds_remaining = seconds_elapsed * ((1.0 - percentage) / percentage);
                eta_str = seconds2str(seconds_remaining);
            end
            text2print = sprintf('%s (%d/%d) %.1f%%%% ETA: %s', ...
                get_progress_bar(percentage, bar_width), ...
                current_progress, total_progress, ...
                percentage * 100, ...
                eta_str);
        end
        
        if last_text_width > 0 && display_mode == 0
            % not first print
            fprintf(repmat('\b', 1, last_text_width - 1));
        end
        fprintf(text2print);
        last_text_width = length(text2print);
    otherwise
        error('Unknown option.');
end

end

function bar = get_progress_bar(percent, width)
n_complete = floor(width * percent);
n_incomplete = width - n_complete;
if n_incomplete == 0
    bar = ['[' repmat('=', 1, n_complete) ']'];
else
    bar = ['[' repmat('=', 1, n_complete) '>' repmat(' ', 1, n_incomplete-1) ']'];
end
end

