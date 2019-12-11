function [results] = my_richardson(phi_not,A,f,err_red,opt)
% [results] = my_richardson(phi_not,A,f,err_red)
%
% perform the 3 step richardson method 
%
% input:
%    phi_not : initial condition vector [Nx1]
%    A       : preconditioned matrix [NxN]
%    f       : preconditioned right hand side vector [Nx1]
%    err_red : error reduction thresholds (err1, err2, ..., err_m),
%              must be specified in going smaller order
%              i.e. err1 > err2 > ... > err_m
%    opts    : additional options, stored as struct
%       max_t   : (optional) max iterations, default: 10000
%       verbose : (optional) return residual history of every iteration, default: 0
%
% output:
%   results.phi_n  : solution at each error tolerance [Nxm]
%   results.n      : number of iterations to reach each err (n1,n2, ... , n_m)

if (~issorted(err_red,'strictdescend'))
    error('error reduction tolerances not in strictly descending order');
end

% check for additional options
if isfield(opt,'verbose')
    verbose = opt.verbose;
else
    verbose = 0;
end
if isfield(opt,'max_t')
    max_t = opt.max_t;
else
    max_t = 10000;
end

% H is the same as point jacobi method
H = -diag(diag(A));

% compute iteration matrix G
N = length(phi_not);
G = H\A;
K = H\f;

% pre-allocate results vecors
m = length(err_red);
n = nan(1,m);
phi_n = nan(N,m);
m_counter = 1;
if verbose
    r = nan(max_t,1);
else
    r = [];
end
% compute initial residual
r_not = norm(A*phi_not - f);
if verbose
    % store residual
    r(1) = r_not;
end
% step sizes
h_s = [4/(6-sqrt(3)) 2/3 4/(6+sqrt(3))];
h_counter = 1;
% iterate
phi_prev = phi_not;
n_counter = 2;
while (n_counter < max_t)
    % perform a forward euler step
    phi_next = phi_prev + h_s(h_counter)*G*phi_prev - h_s(h_counter)*K;
    r_next = A*phi_next - f;
    residual = norm(r_next);
    if verbose
        % store residual
        r(n_counter) = residual;
    end
    if (residual/r_not)<=err_red(m_counter)
        phi_n(:,m_counter) = phi_next;
        n(m_counter) = n_counter-1;
        m_counter = m_counter + 1;
        
        if (m_counter > m)
            % reached end of err_red, break out of loop
            break;
        end
    end

    phi_prev = phi_next;
    n_counter = n_counter + 1;
    h_counter = h_counter + 1;
    if h_counter == 4
        h_counter = 1;
    end
end
if verbose
    r = r(1:n_counter);
end
% store results into a struct
results = struct('phi_n',phi_n,'n',n,'r',r);

end

