
%% convergence analysis
for N = 1:10
    [freq, DOF] = compute_freq_FEM_phs2D(N,N, 1, 0 )
    errorvec(N) = freq-pi;
    DOFvec(N) = DOF; 
end

%%
for N = 1:10
    [freq, DOF] = compute_freq_FEM_phs2D(N,N, 1, 1 )
    errorvecP1P1(N) = freq-pi;
    DOFvecP1P1(N) = DOF; 
end

%%
figure(); loglog(DOFvec, abs(errorvec)); hold all;
loglog(DOFvecP1P1, abs(errorvecP1P1)); legend('P1P0','P1P1')

%%
for N = 1:10
    [freq, DOF] = compute_freq_FEM_phs2D(N,N, 1, 2 )
    errorvecP1P2(N) = freq-pi;
    DOFvecP1P2(N) = DOF; 
end
%%
figure('color','w'); loglog(DOFvec, abs(errorvec)/pi*100); hold all;
loglog(DOFvecP1P1, abs(errorvecP1P1)/pi*100);
loglog(DOFvecP1P2, abs(errorvecP1P2)/pi*100);
grid on;
legend('P1P0P0','P1P1P1', 'P1P2P2')
xlabel('Number of DOF');
ylabel('Error of the first natural frequency (%)')
axis([10 1000 1e-3 100]);
fig = gcf;
set(fig,'Units','Inches');
pos = get(fig,'Position');
fig.PaperPositionMode = 'auto';
set(fig,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)+0.2])
print('convergence2D','-dpdf')