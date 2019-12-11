function [A] = banded_matrix(N,list)
% [A] = banded_matrix(N,list)
%
% create a banded matrix related to a finite difference
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
bands = (1:1:elements) - start_idx;
e = ones(N,1);
A = spdiags(e*list,bands,N,N);

end

