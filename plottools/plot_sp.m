function hf = plot_sp(sp, varargin)
%PLOT_SP Plots spectrum.
options = opt2struct(varargin);
if isfield(options, 'reusefigure') && options.reusefigure
    hf = gcf;
else
    hf = figure;
end
if ~isfield(options, 'plottype')
    options.plottype = 'cartesian';
end
if isfield(sp, 'z')
    % 2d DOA
    error('Not implemented.');
else
    plot_1d_impl(sp, options);
end
end

function [hsp, htp] = plot_1d_impl(sp, options)
    % 1d DOA
    % check sp
    if length(sp.x) ~= size(sp.y, 1)
        error('Length mismatch.');
    end
    % check unit
    switch lower(sp.x_unit)
        case {'radian', 'rad'}
            x_range = [-pi/2 pi/2];
            x_label = '\theta/rad';
        case {'degree', 'deg'}
            x_range = [-90 90];
            x_label = '\theta/deg';
        case {'sin'}
            x_range = [-1 1];
            x_label = 'sin\theta';
        otherwise
            error('Unknown x unit "%s".', sp.x_unit);
    end
    if ~isfield(options, 'usedbscale')
        use_db = false;
    else
        use_db = options.usedbscale;
    end
    if use_db
        y_label = 'Normalized Spectrum (dB)';
    else
        y_label = 'Normalized Spectrum';
    end
    % normalize
    sp.y = sp.y / max(max(sp.y), eps);
    if use_db
        sp.y(sp.y < eps) = eps;
        sp.y = 10*log10(sp.y);
        y_range = [-inf 0];
    else
        y_range = [0 1];
    end
    % plot
    switch lower(options.plottype)
        case 'cartesian'
            if isfield(sp, 'true_positions')
                htp = stem(sp.true_positions, ones(size(sp.true_positions)) * y_range(2), '--r');
                hold on;
            end
            if sp.discrete
                hsp = stem(sp.x, sp.y);
            else
                hsp = plot(sp.x, sp.y);
            end
            hold off;
            if isfield(options, 'xlabel')
                xlabel(options.xlabel);
            else
                xlabel(x_label);
            end
            if isfield(options, 'ylabel')
                ylabel(options.ylabel);
            else
                ylabel(y_label);
            end
            axis([x_range y_range]);
        case 'polar'
            % draw the unit circle
            theta = linspace(0, pi*2, 360);
            plot3(cos(theta), sin(theta), ones(1,360) * y_range(1), '-k');
            hold on;
            if isfield(sp, 'true_positions')
                theta = sp.true_positions / x_range(2) * pi/2;
                htp = stem3(cos(theta), sin(theta), ones(size(sp.true_positions)) * y_range(2), '--r');
                hold on;
            end
            theta = sp.x / x_range(2) * pi/2;
            if sp.discrete
                hsp = stem3(cos(theta), sin(theta), sp.y, 'x');
            else
                hsp = plot3(cos(theta), sin(theta), sp.y);
            end
            hold off;
            axis('equal');
            axis([0 1 -1 1 y_range]);
        otherwise
            error('Unknown plot type "%s".', options.plottype);
    end
    % decorate
    
    if isfield(options, 'title')
        title(options.title);
    else
        title('Spectrum Plot');
    end
end