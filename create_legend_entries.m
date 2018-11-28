function [entry] = create_legend_entries(pre_text,values,style)
%[entry] = create_legend_entries(pre_text,values,style)
%   create entries meant to be passed into legend
%
% example:
%   pre_text = 'Weight'
%   values = [1,2,3]
%   style  = '%d'
%
%  output: entry = {'Weight = 1','Weight = 2','Weight = 3'}
%
% currently no option for post_text, but might be more useful later
%
% Eugene Du
% Feb. 18, 2018



for i = 1: length(values{1})
    entry{i} = '';
    for k = 1 : length(pre_text)
        val = values{k};
        entry{i} = [entry{i} sprintf(['%s = ' style{k}],pre_text{k},val(i))];
    end
end

end

