using PortHamiltonian
using ForwardDiff
Nelements = 20;
Norder1 = 2;
Norder2 = 2;
ph = PortHamiltonian.weak_phs_FEM(Nelements,Norder1+1, Norder2+1,0,1);


function Hamiltonianintegrand(alpha1, alpha2,z)
	0.5*(alpha1^2 + alpha1*alpha2^2)
end
dz = 1/Nelements;
x1,dmy1,dmy2 = PortHamiltonian.lglnodes(Norder1,0,dz);
dz = 1/Nelements;
x2,dmy1,dmy2 = PortHamiltonian.lglnodes(Norder2,0,dz);
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

println(DiscreteHamiltonian(ph.Q*X)-X'*ph.Q*X/2)
grad = x-> ForwardDiff.gradient(DiscreteHamiltonian, ph.Q*x)
hessian = x-> ForwardDiff.hessian(DiscreteHamiltonian, ph.Q*x)
#grad(X)

ph.Hamiltonian = DiscreteHamiltonian;
ph.GradHam = grad;


amplitude = 0.005;
freq = 2*pi/4*2;
using DifferentialEquations
xeq = [1+0.0*cos.(2*pi*linspace(0,1,(Norder1*Nelements+1)));zeros(length(X)-(Norder1*Nelements+1),1)][:];
xeq = inv(ph.Q)*xeq;

function dynamics(xd, x, p, t)
	xd[:] = ph.J * ph.Q*ph.GradHam(x[:]) + amplitude*ph.B[:,1] * (t/(t+1))*sin(freq*t)+amplitude* ph.B[:,2] *(t/(t+1))* sin(freq*t)
end
function dynamics_linear(xd, x, p, t)
	xd[:] = ph.J * ph.Q*x[:] + amplitude*ph.B[:,1] * (t/(t+1)) * sin(freq*t)+amplitude* ph.B[:,2] * (t/(t+1))*sin(freq*t)
end

prob = ODEProblem(dynamics,xeq,(0.,1))
problin = ODEProblem(dynamics_linear,xeq*0,(0.,1))
sol = solve(prob,Trapezoid(),reltol=1e-10,abstol=1e-10,saveat = 0.001);
sollin = solve(problin,Trapezoid(),reltol=1e-10,abstol=1e-10, saveat = 0.001);

youtNL = zeros(length(sol.u[1]),length(sol.u));
for i = 1:length(sol.u)
	youtNL[:,i] = ph.Q*sol.u[i];
end
youtL = zeros(length(sollin.u[1]),length(sollin.u));
for i = 1:length(sollin.u)
	youtL[:,i] = ph.Q*sollin.u[i];
end

using PyPlot;
tvec = sol.t;
xvec = collect(linspace(0,1,(Norder1*Nelements+1)));
TT = tvec' .+ xvec*0;
XX = xvec .+tvec'*0;
figure(11);
surf(XX,TT, youtNL[1:(Norder1*Nelements+1),:], rstride=2,edgecolors="k", cstride=4,
cmap=ColorMap("gray"), alpha=0.6, linewidth=0.1);
ylabel("time (s)"); xlabel("z position (m)"); zlabel("fluid height (m)");
savefig("fluid1DsimulationSmall.jpg");
figure(12);
surf(XX,TT,youtNL[(Norder1*Nelements+2):end,:], rstride=2,edgecolors="k", cstride=2,
cmap=ColorMap("gray"), alpha=0.8, linewidth=0.25)
youtNL = youtNL';
youtL = youtL';
figure(14);
subplot(2,2,1);
plot(xvec, youtNL[200,1:(Norder1*Nelements+1)]);
plot(xvec, youtL[200,1:(Norder1*Nelements+1)]+1,  linestyle = "dashed");
title("t = 0.2s")
xlabel("z position (m)"); ylabel("fluid height (m)");
subplot(2,2,2);
plot(xvec, youtNL[400,1:(Norder1*Nelements+1)]);
plot(xvec, youtL[400,1:(Norder1*Nelements+1)]+1, linestyle = "dashed");
title("t = 0.4s")
xlabel("z position (m)"); ylabel("fluid height (m)");
subplot(2,2,3);
plot(xvec, youtNL[600,1:(Norder1*Nelements+1)]);
plot(xvec, youtL[600,1:(Norder1*Nelements+1)]+1,
linestyle = "dashed");
title("t = 0.6s")
xlabel("z position (m)"); ylabel("fluid height (m)");
subplot(2,2,4);
plot(xvec, youtNL[800,1:(Norder1*Nelements+1)]);
plot(xvec, youtL[800,1:(Norder1*Nelements+1)]+1,linestyle = "dashed");
title("t = 0.8s")
xlabel("z position (m)"); ylabel("fluid height (m)");
legend(["nonlinear", "linear"])
tight_layout();
filename = "simulations_snapshot_small";
savefig(string(filename,".jpg"));
