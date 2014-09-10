files = dir('Runs/Stage8/Run3/training_data_stage8/representors/*.txt');


for i=1:length(files)
    d = importdata(strcat('Runs/Stage8/Run3/training_data_stage8/representors/', files(i).name), ' ' , 1);
    [pathstr, name, ext] = fileparts(files(i).name);
    v = genvarname(name);
    eval([v '= d.data;']);
end

cumul_data = F_repr_before_it_5;
Nmax = (length(files)-1)/4;
unique_sizes = zeros(1,Nmax);
for N = 1:Nmax
    
    if N == 5
        cumul_data = [cumul_data; F_representor_1];
    end
    
    for i=1:4
        Arepr_name = genvarname(sprintf('A_representor_%u_%u', i, N));
        eval(['cumul_data =[cumul_data;' Arepr_name '];']);
    end
    
    size(cumul_data)
    
    unique_sizes(N) = size(unique(cumul_data(:, 4:9), 'rows' ),1);
    
end

data = [(1:Nmax)', unique_sizes']