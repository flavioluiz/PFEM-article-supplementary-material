using PortHamiltonian

Nelem_vec = 2:4:100;

#Norder1_vec = [1 1 1 2 2 2 2 3 3 3 3 4];
#Norder2_vec = [0 1 2 0 1 2 3 0 1 2 3 4];
Norder1_vec = [1 1 1 3];
Norder2_vec = [0 1 2 3];
f1 = zeros(length(Nelem_vec), length(Norder1_vec));
dof = zeros(length(Nelem_vec), length(Norder1_vec));
using PyPlot
for jj = 1:length(Norder1_vec)
	Norder1 = Norder1_vec[jj];
	Norder2 = Norder2_vec[jj];	
	i = 1;
	for Nelements = Nelem_vec	
		ph = PortHamiltonian.weak_phs_FEM(Nelements,Norder1+1, Norder2+1,0,1);
		freq = frequencies(ph);
		dof[i,jj] = size(ph.Q,1);
		f1[i,jj] = abs(freq[freq .> 0.1][2]-2*pi);
		i = i+1;
	end
	figure(1);
	loglog(dof[:,jj], f1[:,jj])
	figure(2);
	loglog(Nelem_vec, f1[:,jj])
end

ab =[];
for i  = 1:length(Norder1_vec)
	ab = [ab; string("P",Norder1_vec[i],"P", Norder2_vec[i])];
end
figure(1);
legend(ab);
xlabel("Number of DOF");
ylabel("Error of first natural frequency (%)");
ax1 = subplot(1,1,1)
ax1[:xaxis][:set_ticks]([collect(10:10:100); collect(100:100:1000)])
ax1[:set_xlim](xmin = 10, xmax = 1000)
grid("on")
filename = "convergence1D_nDOF";
savefig(string(filename,".pdf"));

figure(2);
legend(ab);
xlabel("Number of elements");
ylabel("Error of first natural frequency");
tight_layout();
filename = "convergence1D_Nelem_5thmode";
savefig(string(filename,".pdf"));
