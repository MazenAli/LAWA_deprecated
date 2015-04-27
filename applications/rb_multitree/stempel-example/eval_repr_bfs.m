clear all, close all, clc

%files = dir('Runs/Runs_NxD/Run10/training_data_stempel/representors/*.txt');

files = dir('Runs/Runs_NxD/Run1b/training_data_stempel/representors/*.txt');
for i=1:length(files)
    %d = importdata(strcat('Runs/Runs_NxD/Run10/training_data_stempel/representors/', files(i).name), ' ' , 1);
    d = importdata(strcat('Runs/Runs_NxD/Run1b/training_data_stempel/representors/', files(i).name), ' ' , 1);
    [pathstr, name, ext] = fileparts(files(i).name);
    v = genvarname(name);
    eval([v '= d.data;']);
end

Nmax = (length(files)-9)/4;


%%
files = dir('Runs/Runs_NxD/Run1b/training_data_stempel/representors/F_repr_before_it_12/*.txt');
for i=1:length(files)
    %d = importdata(strcat('Runs/Runs_NxD/Run10/training_data_stempel/representors/', files(i).name), ' ' , 1);
    d = importdata(strcat('Runs/Runs_NxD/Run1b/training_data_stempel/representors/F_repr_before_it_12/', files(i).name), ' ' , 1);
    [pathstr, name, ext] = fileparts(files(i).name);
    v = genvarname(strcat('Before_It_12_',name));
    eval([v '= d.data;']);
end

files = dir('Runs/Runs_NxD/Run1b/training_data_stempel/representors/F_repr_before_it_13/*.txt');
for i=1:length(files)
    %d = importdata(strcat('Runs/Runs_NxD/Run10/training_data_stempel/representors/', files(i).name), ' ' , 1);
    d = importdata(strcat('Runs/Runs_NxD/Run1b/training_data_stempel/representors/F_repr_before_it_13/', files(i).name), ' ' , 1);
    [pathstr, name, ext] = fileparts(files(i).name);
    v = genvarname(strcat('Before_It_13_',name));
    eval([v '= d.data;']);
end

files = dir('Runs/Runs_NxD/Run1b/training_data_stempel/representors/F_repr_before_it_14/*.txt');
for i=1:length(files)
    %d = importdata(strcat('Runs/Runs_NxD/Run10/training_data_stempel/representors/', files(i).name), ' ' , 1);
    d = importdata(strcat('Runs/Runs_NxD/Run1b/training_data_stempel/representors/F_repr_before_it_14/', files(i).name), ' ' , 1);
    [pathstr, name, ext] = fileparts(files(i).name);
    v = genvarname(strcat('Before_It_14_',name));
    eval([v '= d.data;']);
end


%%

cumul_data = [];
for i=1:9
    Frepr_name = genvarname(sprintf('Before_It_12_F_representor_%u', i));
    eval(['cumul_data =[cumul_data;' Frepr_name '];']);
end

size(cumul_data)

unique_sizes = zeros(1,Nmax);
for N = 1:Nmax
    
    if N == 12
        for i=1:9
            Frepr_name = genvarname(sprintf('Before_It_13_F_representor_%u', i));
            eval(['cumul_data =[cumul_data;' Frepr_name '];']);
        end
    end
    if N == 13
        for i=1:9
            Frepr_name = genvarname(sprintf('Before_It_14_F_representor_%u', i));
            eval(['cumul_data =[cumul_data;' Frepr_name '];']);
        end
    end
    if N >= 14
        for i=1:9
            Frepr_name = genvarname(sprintf('F_representor_%u', i));
            eval(['cumul_data =[cumul_data;' Frepr_name '];']);
        end
    end
    
    for i=1:4
        Arepr_name = genvarname(sprintf('A_representor_%u_%u', i, N));
        eval(['cumul_data =[cumul_data;' Arepr_name '];']);
    end
    
    size(cumul_data)
    
    unique_sizes(N) = size(unique(cumul_data(:, 4:9), 'rows' ),1);
    
end

data = [(1:Nmax)', unique_sizes']