Errs = importdata('Acc1em4_Stage7_Run3_OnlineTest_Errors.txt', ' ', 1)
Errs = Errs.data
ErrEsts = importdata('Acc1em4_Stage7_Run3_OnlineTest_ErrorEstimators.txt', ' ', 1)
ErrEsts = ErrEsts.data

figure()

for N=1:7
    
    subplot(3,3, N)
    plot([1:20]', Errs(:,N+2), 'r-o');
    hold on
    plot([1:20]', ErrEsts(:,N+2), 'b-x');
    %legend('Errors', 'Error Bounds')
    title(sprintf('N= %u', N));
    hold off

end