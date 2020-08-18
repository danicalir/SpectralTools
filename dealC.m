function [varargout] = dealC(C)

% deals output to multiple variables (use in conjuction with cellfun)
%
% If X is a 1x4 cell array, the following returns each cell's first element:
%       [X1,X2,X3,X4]=dealC(cellfun(@(a) a(1),X,'UniformOutput',0));
%
% Note: this replaces anonymous function DC=@(C) deal(C{:});

v = [];
for i=1:length(C)
    v = [v ' varargout{',num2str(i),'}'];
end
eval(['[',v,'] = deal(C{:});']);