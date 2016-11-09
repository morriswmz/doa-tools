function s = opt2struct(opt_args)
%OPT2STRUCT Converts options specified by "varargin" to a struct.
s = struct();
if isempty(opt_args)
    return;
end
if ~iscell(opt_args)
    error('Cell array expected.');
elseif ~isvector(opt_args)
    error('1D cell array expected.');
elseif mod(length(opt_args),2) ~= 0
    error('Even number of elements expected');
end

for i = 1:2:length(opt_args)
    f = lower(opt_args{i});
    if isvarname(f)
        s.(f) = opt_args{i+1};
    end
end
end

