% add folders to search path
cur_dir = pwd;
folders = {'utils', 'array', 'estimator', 'plottool', 'performance', 'solvers'};
for ii = 1:length(folders)
    addpath(genpath(fullfile(cur_dir, folders{ii})));
end