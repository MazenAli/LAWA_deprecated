
Nmax = 10;
ntrain = 400;
format long

datastr = 'Runs/Runs_NxD/Run7c/stempel_greedy_info.txt_eps_';

Deltas = zeros(ntrain, Nmax);

for N = 0:Nmax
   eps_aff = load(strcat(datastr, 'aff_N_', num2str(N), '.txt')); 
   eps_bound = load(strcat(datastr, 'bound_N_', num2str(N), '.txt')); 
   
   Deltas(:,N+1) = eps_aff./eps_bound;
end

savestr = strcat(datastr, 'AllDeltas.txt');
save(savestr, 'Deltas', '-ascii')

MaxDeltas = max(Deltas)
MinDeltas = min(Deltas)
AvDeltas = mean(Deltas)
Data = [(1:Nmax+1)', MinDeltas', AvDeltas', MaxDeltas'];
savestr_max = strcat(datastr, 'AvDeltas.txt');
save(savestr_max, 'Data', '-ascii')
    