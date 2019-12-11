function [A] = banded_matrix_periodic(N,list)
% [A] = banded_matrix_periodic(N,list)
%
% create a banded periodic matrix related to a finite difference
% approximation scheme
%
% input:
%    N     : desired size of matrix
%    list  : stencil of matrix, size must be odd 
%
% output:
%   A      : output banded matrix [NxN]

elements = length(list);
% check that length of list is odd
if ~bitand(elements,1)
    error('require odd number of elements in list');
end
if elements > N
    error('number of elements in list exceeeds N');
end

start_idx = ceil(elements/2);
temp2 = 1:1:elements;
j_s = nan(N,elements);
for i = 1 : N
    j_s(i,:) = mod(temp2-start_idx+i-1,N)+1;
end
i_s = repmat((1:1:N)',[1,elements]);
s_s = repmat(list,[N,1]);
A = sparse(i_s,j_s,s_s);

end

