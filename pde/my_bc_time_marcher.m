function [u_final,u] = my_bc_time_marcher(u_not,A,dt,T,bc,bc_f,method)
% [u_final,u] = my_time_marcher(u_not,A,T,dt,method)
%
% general time marching function to perform many time marching methods
% assuming system of the form du/dt = Au + bc.*f(t), where A is independent of u,t
% and bc is indepent of u,t
% input:
%    u_not : initial condition vector [nx1]
%    A     : "stiffness"? matrix [nxn]
%    dt    : time step
%    T     : end time assuming time marching from 0:dt:T
%    bc    : coeffiecients to bc.*f(t) [nx1]
%    bc_f  : handle to function f(t), assumes f is continuous from [0:T]
%    method: string indicating marching method
%
% output:
%   u_final: u(T) [nx1]
%   u      : time history of u [nxlength(0:dt:T)]

N = round(T/dt) + 1;
n = length(u_not);
u = nan(n,N);
u(:,1) = u_not;

switch method
    case 'fe' % forward euler
        for i = 2 : N
            u(:,i) = u(:,i-1) + dt*A*u(:,i-1) + dt*bc_f((i-2)*dt)*bc;
        end
    case 'be' % backward euler
        for i = 2 : N
            u(:,i) = (eye(n)-dt*A)\(u(:,i-1) + dt*bc_f((i-1)*dt)*bc);
        end
    case 'ab2' % AB2

    case 'trap' % trapezoidal

    case 'rk4' % rk4
        
        for i = 2 : N
            % confusing but matlab uses 1 as first index
            % in time marching we assume t0 = 0
            bc_i = bc_f((i-2)*dt)*bc;
            bc_i_half = bc_f((i-1.5)*dt)*bc;
            bc_i_one = bc_f((i-1)*dt)*bc;
            
            v1 = u(:,i-1) + dt/2*(A*u(:,i-1) + bc_i);
            v2 = u(:,i-1) + dt/2*(A*v1 + bc_i_half);
            v3 = u(:,i-1) + dt*(A*v2 + bc_i_half);
            
            u(:,i) = u(:,i-1) + 1/6*dt*((A*u(:,i-1) + bc_i) + 2*(A*v1 + bc_i_half) + 2*(A*v2 + bc_i_half)+ (A*v3 + bc_i_one));
        end
    otherwise
      error('Unrecognized time marchining method')  
end

u_final = u(:,end);
end

