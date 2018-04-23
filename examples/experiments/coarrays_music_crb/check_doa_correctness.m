function correct = check_doa_correctness(actual, est, tolerance)
%CHECK_DOA_CORRECTNESS Checks if estimated doa is reasonably close to the
%actual ones. Input must be in ascending order within range of
%(-pi/2, pi/2).
%More precisely, suppose the true DOAs are x_1, x_2, ..., x_k, the estimated
%DOAs must fall within the following regions:
% (-pi/2, x_1 + (x_2 - x_1) / 2),
% (x_1 + (x_2 - x_1) / 2, x_2 + (x_3 - x_2) / 2), ...
% (x_k - (x_k - x_{k-1}) / 2, pi/2)
%If tolerance is not inf, the estimated DOAs must also fall within the following
% regions:
% (x_1 - tolerance, x_1 + tolerance), ...
% (x_k - tolerance, x_k + tolerance)

if nargin == 2
    tolerance = inf;
end
n = length(actual);
if n ~= length(est)
    correct = false;
    return;
end
if n == 1
    correct = abs(actual(1) - est(1)) < tolerance;
elseif n > 1
    spaces = diff(actual);
    max_deviation = min(spaces/2, tolerance);
    correct = true;
    % first
    correct = correct && ...
        (est(1) >= max(-pi/2, actual(1)-tolerance) && ...
        est(1) < actual(1) + max_deviation(1)); 
    % last
    correct = correct && ...
        (est(end) <= min(pi/2, actual(end)+tolerance) && ...
        est(end) > actual(end) - max_deviation(end));
    % middle
    for ii = 2:(n-1)
        correct = correct && ...
        (est(ii) > actual(ii) - max_deviation(ii-1) && ... 
        est(ii) < actual(ii) + max_deviation(ii));
    end
end

end

