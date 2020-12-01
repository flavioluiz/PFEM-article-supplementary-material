using PortHamiltonian
Nelements = 20;
Norder1 = 2;
Norder2 = 2;
ph = PortHamiltonian.weak_phs_FEM(Nelements,Norder1+1, Norder2+1,0,1);


function Hamiltonianintegrand(alpha1, alpha2,z)
	0.5*(alpha1^2 + alpha1*alpha2^2)
end
dz = 1/Nelements;
x1,~,~ = PortHamiltonian.lglnodes(Norder1,0,dz);
dz = 1/Nelements;
x2,~,~ = PortHamiltonian.lglnodes(Norder2,0,dz);
elements_domain = zeros(Nelements, 2);
for i = 1:Nelements
	elements_domain[i,:] = [dz*(i-1) dz*i];
end
z_int,w_int = PortHamiltonian.lgwt(Norder1+Norder2+2,0,dz);


function pol1(z)
	map(i->PortHamiltonian.leg_pol(z, x1, i), 1:(Norder1+1));
end
function pol2(z)
	map(i->PortHamiltonian.leg_pol(z, x2, i), 1:(Norder2+1));
end

M1_el = zeros(length(x1), length(x1));
for i = 1:length(x1)
	for j = 1:length(x1)
		for jj = 1:length(z_int)
		M1_el[i,j] = M1_el[i,j] + w_int[jj]*PortHamiltonian.leg_pol(z_int[jj], x1, i)*PortHamiltonian.leg_pol(z_int[jj], x1, j);
		end
	end
end
M1 = zeros((Norder1*Nelements+1),(Norder1*Nelements+1));
for i = 1:(Nelements)
	j = 1 + Norder1*(i-1);
	M1[j:(j+Norder1),j:(j+Norder1)] += M1_el;
end

M2_el = zeros(length(x2), length(x2));
for i = 1:length(x2)
	for j = 1:length(x2)
		for jj = 1:length(z_int)
		M2_el[i,j] = M2_el[i,j] + w_int[jj]*PortHamiltonian.leg_pol(z_int[jj], x2, i)*PortHamiltonian.leg_pol(z_int[jj], x2, j);
		end
	end
end


M2 = zeros((Norder2*Nelements+1),(Norder2*Nelements+1));
for i = 1:(Nelements)
	j = 1 + Norder2*(i-1);
	M2[j:(j+Norder2),j:(j+Norder2)] += M2_el;
end
el_mat2 = zeros(Norder2+1, Nelements)
for i = 1:(Nelements)
	j = 1 + Norder2*(i-1);
	el_mat2[:,i] = j:(j+Norder2)
end
el_mat1 = zeros(Norder1+1, Nelements)
for i = 1:(Nelements)
	j = 1 + Norder1*(i-1);
	el_mat1[:,i] = j:(j+Norder1)
end
function M122fun(k)
	
	M122 = zeros((Norder2*Nelements+1),(Norder2*Nelements+1));	
	ind_k = ind2sub(el_mat1, find(el_mat1.==k));
	#ind_k = el_mat1[el_mat1.==k];
	for ii =1:length(ind_k)
		
		kk = Integer(ind_k[ii][1]);
		M122_el = zeros(length(x2), length(x2));	
		#print(kk," ",ind_k[2][ii],"\n");
		for i = 1:length(x2)
			for j = 1:length(x2)
				for jj = 1:length(z_int)
				M122_el[i,j] = M122_el[i,j] + w_int[jj]*PortHamiltonian.leg_pol(z_int[jj], x1, kk)*PortHamiltonian.leg_pol(z_int[jj], x2, i)*PortHamiltonian.leg_pol(z_int[jj], x2, j);
				end
			end
		end
		j = 1 + Norder2*(ind_k[ii][2]-1);
		M122[j:(j+Norder2),j:(j+Norder2)] += M122_el;
	end
	M122
end

# the following function computes the nonlinear Hamiltonian gradient
# from the formulas presented in the paper
function gradHamiltonian(X)
	alpha1 = X[1:(Norder1*Nelements+1)];
	alpha2 = X[(Norder1*Nelements+2):end];
	gradNL = zeros(size(alpha1));
	gradNL2 = zeros(size(alpha2));
	for i = 1:length(gradNL)
		gradNL[i] = [0.5*alpha2'*M122fun(i)*alpha2][1];
	end
	for i = 1:length(gradNL2)
		gradNL2 = gradNL2 + alpha1[i]*M122fun(i)*alpha2;
	end
	[M1*alpha1+gradNL;
	 M2*alpha2*0 + gradNL2]
end


function DiscreteHamiltonian(X)
	alpha1 = X[1:(Norder1*Nelements+1)];
	alpha2 = X[(Norder1*Nelements+2):end];
	Hd = 0;	
	for i = 1:Nelements
		alpha1_el = alpha1[1+(i-1)*Norder1:1+i*Norder1];
		alpha2_el = alpha2[1+(i-1)*Norder2:1+i*Norder2];
		alpha1_elz = z-> (alpha1_el'*pol1(z))[1];
		alpha2_elz = z-> (alpha2_el'*pol2(z))[1];
		for j = 1:length(z_int)
			#Hd = Hd + w_int[j]'*Hamiltonianintegrand(alpha1_elz(z_int[j]), alpha2_elz(z_int[j]),z_int[j]+elements_domain[j,1]);
			Hd = Hd + w_int[j]'*Hamiltonianintegrand(alpha1_elz(z_int[j]), alpha2_elz(z_int[j]),z_int[j]*0);
		end
	end
	Hd
end
X = randn(size(ph.Q,1));

#grad = x-> ForwardDiff.gradient(DiscreteHamiltonian, ph.Q*x)
#hessian = x-> ForwardDiff.hessian(DiscreteHamiltonian, ph.Q*x)
#grad(X)

ph.Hamiltonian = DiscreteHamiltonian;
ph.GradHam = gradHamiltonian;

amplitude = 0.5;
freq = 2*pi/4*2;
function dynamics(t,x, xd)
	xd[:] = ph.J * ph.Q*ph.GradHam(ph.Q*x[:]) + amplitude*ph.B[:,1] *(t/(t+1))* sin(freq*t)+amplitude* ph.B[:,2] *(t/(t+1))* sin(freq*t)
end


function dynamics_linear(t,x, xd)
	xd[:] = ph.J * ph.Q*x[:] + amplitude*ph.B[:,1] *(t/(t+1))* sin(freq*t)+amplitude* ph.B[:,2] *(t/(t+1))* sin(freq*t)
end
using Sundials
xeq = [1+0.0*cos.(2*pi*linspace(0,1,(Norder1*Nelements+1)));zeros(length(X)-(Norder1*Nelements+1),1)][:];
xeq = inv(ph.Q)*xeq;
t = (linspace(0,1,500))
yout = Sundials.cvode(dynamics,xeq, collect(t))
yout_linear = Sundials.cvode(dynamics_linear,xeq*0, collect(t))
youtn = yout*(ph.Q);
youtn_linear = yout_linear*(ph.Q);
using PyPlot
figure(1)
surf(youtn[:,1:(Norder1*Nelements+1)]'); 
surf(youtn_linear[:,1:(Norder1*Nelements+1)]'+1); 
figure(2)
surf(youtn[:,(Norder1*Nelements+2):end])
surf(youtn_linear[:,(Norder1*Nelements+2):end])

figure(4);
subplot(2,2,1);
plot(youtn[100,1:(Norder1*Nelements+1)]);
plot(youtn_linear[100,1:(Norder1*Nelements+1)]+1);
subplot(2,2,2);
plot(youtn[200,1:(Norder1*Nelements+1)]);
plot(youtn_linear[200,1:(Norder1*Nelements+1)]+1);
subplot(2,2,3);
plot(youtn[300,1:(Norder1*Nelements+1)]);
plot(youtn_linear[300,1:(Norder1*Nelements+1)]+1);
subplot(2,2,4);
plot(youtn[400,1:(Norder1*Nelements+1)]);
plot(youtn_linear[400,1:(Norder1*Nelements+1)]+1);
legend("nonlinear", "linear")
# #yout, ypout = Sundials.idasol(dynamics2, [xeq;l0;0], [xeq;0;0]*0, [t])
# pfluid.GradHam = x-> Qlin*x
# youtlin = Sundials.cvode(dynamics,xeq, [t])
# #youtlin, ypoutlin = Sundials.idasol(dynamics2, [xeq;l0; 0], [xeq;0; 0]*0, [t])
# #pn = PortHamiltonian.constraint_elimination(pfluid)
 #abs(err[2]) < 1e-13

#ph_c_elim = constraint_elimination(ph_constrained)
#eigval, eigvec = eig(ph_c_elim);
#num_freq = frequencies(eigval)
#exact_freq = pi*collect(0:1:N-1);
#err = exact_freq - num_freq;

#@test abs(err[2]) < 1e-13