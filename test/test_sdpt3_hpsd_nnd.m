clear();
rng(1337);

fprintf('Detecting CVX... ');
test_against_cvx = true;
try
    cvx_clear;
    cvx_begin;
    cvx_end;
    fprintf('Success!\n');
catch ex
    test_against_cvx = false;
    fprintf('CVX not available. Skipping tests using CVX.\n');
end

%% Test 1: identity matrices
% they have closed form solutions
lambda = 0.5;
fprintf('Testing using identity matrices...\n');

for n = 8:4:24
    R0 = eye(n);
    
    tolerance = 1e-8;
    R = sdpt3_hpsd_nnd(R0, lambda);
    R_truth = eye(n) * (1 - lambda/sqrt(n));
    discrepency = norm(R(:) - R_truth(:), 1);
    fprintf('Dim %d Constrainted: %.3e\n', n, discrepency);
    assert(discrepency < tolerance, sprintf('Discrepency %.3e (> %.3e) is too large.', discrepency, tolerance));
    
    tolerance = 1e-3;
    R = sdpt3_hpsd_nnd(R0, lambda, 'Formulation', 'Penalized');
    R_truth = eye(n) * (1 - lambda/2);
    discrepency = norm(R(:) - R_truth(:), 1);
    fprintf('Dim %d Penalized: %.3e\n', n, discrepency);
    assert(discrepency < tolerance, sprintf('Discrepency %.3e (> %.3e) is too large.', discrepency, tolerance));
end

%% Test 2: comparison with CVX, unweighted
if test_against_cvx
    fprintf('Testing against results obtained from CVX...\n');
    % test 5 cases
    m = 20;
    l = 5;
    noise_var = 0.01;
    lambda = 4;
    for ii = 1:5
        L = randn(m, l);
        R0 = L*L';
        R = R0 + randn(m, m) * sqrt(noise_var);
        
        cvx_clear();
        cvx_begin sdp quiet;
            variable R1(m,m) complex semidefinite;
            minimize norm_nuc(R1);
                subject to
                norm(R(:) - R1(:)) <= lambda;
        cvx_end;
        [R2, info] = sdpt3_hpsd_nnd(R, lambda);
        
        discrepency = norm(R1(:) - R2(:), 1);
        obj_rel_diff = (cvx_optval - info.obj(1)) / cvx_optval * 100;
        fprintf('Case %d Constrained: || vec(R) - vec(R_cvx) ||_1 = %.3e, obj_rel_diff = %.2f%%\n', ...
            ii, discrepency, obj_rel_diff);
        fprintf('        CPU time: %fs (CVX = %fs)\n', info.cputime, cvx_cputime);
        assert(obj_rel_diff < 1e-5, 'Objective function differs too much.');
        
        cvx_clear();
        cvx_begin sdp quiet;
            variable R1(m,m) complex semidefinite;
            minimize lambda * norm_nuc(R1) + sum_square_abs(R(:) - R1(:));
        cvx_end;
        [R2, info] = sdpt3_hpsd_nnd(R, lambda, 'Formulation', 'Penalized');
        
        discrepency = norm(R1(:) - R2(:), 1);
        obj_rel_diff = (cvx_optval - info.obj(1)) / cvx_optval * 100;
        fprintf('Case %d Penalized: || vec(R) - vec(R_cvx) ||_1 = %.3e, obj_rel_diff = %.2f%%\n', ...
            ii, discrepency, obj_rel_diff);
        fprintf('        CPU time: %fs (CVX = %fs)\n', info.cputime, cvx_cputime);
        assert(obj_rel_diff < 1e-5, 'Objective function differs too much.');
    end
end

%% Test 3: comparison with CVX, weighted
if test_against_cvx
    fprintf('Testing scenarios with a weighting matrix...\n');

    lambda = 2;
    r = exp(1j*((1:20)'*[0.1 0.2 0.3])*pi);
    R0 = r * r'; % R0 is low rank
    R0 = 0.5*(R0 + R0');
    m = size(R0, 1);
    W = kron(R0.', R0) / 100 + eye(m*m) / 100;
    W = W / norm(W, 'fro');
    W = 0.5*(W+W');
    W_sqrt = sqrtm(W);
    noise = W_sqrt * (randn(m*m,1) + 1j*randn(m*m,1)) / sqrt(2);
    R = R0 + reshape(noise, m, m);
    W_est = W;

    W_chol = chol(W_est);
    cvx_clear();
    cvx_begin sdp quiet;
        variable R1(m,m) complex semidefinite;
        minimize norm_nuc(R1);
        subject to
            norm(W_chol*(R1(:) - R(:))) <= lambda;
    cvx_end;
    [R2, info] = sdpt3_hpsd_nnd(R, lambda, 'W', W_est);

    discrepency = norm(R2(:) - R1(:), 1);
    obj_rel_diff = (cvx_optval - info.obj(1)) / cvx_optval * 100;
    fprintf('Constrained: || vec(R) - vec(R_cvx) ||_1 = %.3e, obj_rel_diff = %.2f%%\n', ...
                discrepency, obj_rel_diff);
    fprintf('             CPU time: %fs (CVX = %fs)\n', info.cputime, cvx_cputime);
    assert(obj_rel_diff < 1e-5, 'Objective function differs too much.');

    cvx_clear();
    cvx_begin sdp quiet;
        variable R1(m,m) complex semidefinite;
        minimize lambda * norm_nuc(R1) + sum_square_abs(W_chol*(R1(:) - R(:)));
    cvx_end;
    [R2, info] = sdpt3_hpsd_nnd(R, lambda, 'W', W_est, 'Verbose', false, ...
        'Formulation', 'Penalized');

    discrepency = norm(R2(:) - R1(:), 1);
    obj_rel_diff = (cvx_optval - info.obj(1)) / cvx_optval * 100;
    fprintf('Penalized: || vec(R) - vec(R_cvx) ||_1 = %.3e, obj_rel_diff = %.2f%%\n', ...
                discrepency, obj_rel_diff);
    fprintf('           CPU time: %fs (CVX = %fs)\n', info.cputime, cvx_cputime);
    assert(obj_rel_diff < 1e-5, 'Objective function differs too much.');
end
    