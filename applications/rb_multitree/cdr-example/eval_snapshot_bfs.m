
files = dir('Runs/Stage6/Run4/offline_data_stage6/bf_*.txt');

Cumul_BF_0 = [];
sizes = zeros(size(files));
unique_sizes = zeros(size(files));
intersect_sizes = zeros(size(files));

row_str = 'rows';
for i=1:length(files)
    bf = importdata(strcat('Runs/Stage6/Run4/offline_data_stage6/', files(i).name), ' ' , 1);
    [pathstr, name, ext] = fileparts(files(i).name);
    v = genvarname(name);
    eval([v '= bf.data;']);
    sizes(i) = eval(['size(' v ', 1);']);
    
    oldcumulname = genvarname(sprintf('Cumul_BF_%u', i-1));
    cumulname = genvarname(sprintf('Cumul_BF_%u', i));
    eval([cumulname '= [' oldcumulname '; bf.data];']);
    
    unique_sizes(i) = eval(['size(unique(' cumulname '(:, 4:9), row_str ),1);']);
    
    oldisectname = genvarname(sprintf('Intersect_BF_%u', i-1));
    isectname = genvarname(sprintf('Intersect_BF_%u', i));
    
    if i==1
        eval([isectname '= [bf.data(:,4:9)];']);
    else
        eval([isectname '= intersect(bf.data(:, 4:9), ' oldisectname ', row_str );']);
    end
    intersect_sizes(i) = eval(['size(' isectname ', 1);']);
    
end

data = [(1:length(files))', sizes, unique_sizes, intersect_sizes]
