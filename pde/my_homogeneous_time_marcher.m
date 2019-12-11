function [u_final,u] = my_homogeneous_time_marcher(u_not,A,dt,T,method)
% [u_final,u] = my_time_marcher(u_not,A,T,dt,method)
%
% general time marching function to perform many time marching methods
% assuming system of the form du/dt = Au, where A is independent of u,t
%
% input:
%    u_not : initial condition vector [nx1]
%    A     : "stiffness"? matrix [nxn]
%    dt    : time step
%    T     : end time assuming time marching from 0:dt:T
%    method: string indicating marching method
%
% output:
%   u_final: u(T) [nx1]
%   u      : time history of u [nxlength(0:dt:T)]

N = length(0:dt:T);
n = length(u_not);
u = nan(n,N);
u(:,1) = u_not;

switch method
    case 'fe' % forward euler
        for i = 2 : N
            u(:,i) = u(:,i-1) + dt*A*u(:,i-1);
        end
    case 'be' % backward euler
        for i = 2 : N
            u(:,i) = (eye(n)-dt*A)\u(:,i-1);
        end
    case 'ab2' % AB2
        % as with insturctions from textbook, initialize u to be very small
        u_not2(:,1) = u_not/norm(u_not)/10;
        for i = 2 : N
            if i==2
                u(:,i) = u(:,i-1)+dt/2*(3*A*u(:,i-1)- A*u_not2);
            else
                u(:,i) = u(:,i-1)+dt/2*(3*A*u(:,i-1)- A*u(:,i-2));
            end
        end
    case 'trap' % trapezoidal
        for i = 2 : N
            u(:,i) = (eye(n)-dt/2*A)\(eye(n)+dt/2*A)*u(:,i-1);
        end
    case 'rk4' % rk4
        for i = 2 : N
            v1 = u(:,i-1) + 0.5*dt*A*u(:,i-1);
            v2 = u(:,i-1) + 0.5*dt*A*v1;
            v3 = u(:,i-1) + dt*A*v2;
            
            u(:,i) = u(:,i-1) + (1/6)*dt*(A*u(:,i-1) + 2*A*(v1 + v2) + A*v3);
        end
    otherwise
      error('Unrecognized time marchining method')  
end

u_final = u(:,end);
end

