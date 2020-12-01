using PortHamiltonian
using ForwardDiff
# domain dimensions
Lx = 1;
Ly = Lx; # for now, only Lx+Ly is allowed

# number of elements
Nex = 10;
Ney = Nex; # for now, only equal number of elements is allowed

# degree of the polynomial - elements
Pe1 = 1;
Pe23 = 0;
N1 = Pe1 + 1;
N2 = Pe23 + 1;
N3 = Pe23 + 1;

# element size
dx = Lx / Nex;
dy = Ly /Ney;

b = Lx;
a = 0;

x1,w1,P = PortHamiltonian.lglnodes(N1-1,(b-a)/Nex,0)    # discretization of x1 variables
x1i,w1i,Pi = PortHamiltonian.lglnodes(2*N1,(b-a)/Nex,0)  # these points are used to integrate the mass matrix using quadrature
if N2 > 1
    x2,w2 = PortHamiltonian.lglnodes(N2-1,(b-a)/Nex,0)
else
    x2,w2 = [(a+b)/2], [(b-a)/Nex];
end
x2i,w2i,Pi = PortHamiltonian.lglnodes(2*N2,(b-a)/Nex,0)  # these points are used to integrate the mass matrix using quadrature

if N3 > 1
    x3,w3 = PortHamiltonian.lglnodes(N3-1,(b-a)/Nex,0)
else
    x3,w3 = [(a+b)/2], [(b-a)/Nex];
end
x3i,w3i,Pi = PortHamiltonian.lglnodes(2*N3,(b-a)/Nex,0)  # these points are used to integrate the mass matrix using quadrature


xi,wi,Pi = PortHamiltonian.lglnodes(2*N1,(b-a)/Nex,0)
#M = massmatrix(ei,xi,we);
#xquad, wquad = lgwt(Ne,a,b)\
#M = massmatrix(ei,xi, xquad, wquad);

# basis functions for each element
phi1 = x-> map(i-> PortHamiltonian.leg_pol(x, x1,i), 1:length(x1))
phi1xy = (x,y) -> (phi1(x)*phi1(y)')'[:]


phi2 = x-> map(i-> PortHamiltonian.leg_pol(x, x2,i), 1:length(x2))
phi2xy = (x,y) -> (phi2(x)*phi2(y)')'[:]


phi3 = x-> map(i-> PortHamiltonian.leg_pol(x, x3,i), 1:length(x3))
phi3xy = (x,y) -> (phi3(x)*phi3(y)')'[:]


M1 = zeros(N1*N1,N1*N1);
for i = 1:length(x1i)
    for j = 1:length(x1i)
        M1 = M1+w1i[i]*w1i[j]*phi1xy(x1i[i],x1i[j])*phi1xy(x1i[i],x1i[j])';
    end
end


M2 = zeros(N2*N2,N2*N2);
for i = 1:length(x2i)
    for j = 1:length(x2i)
        M2 = M2+w2i[i]*w2i[j]*phi2xy(x2i[i],x2i[j])*phi2xy(x2i[i],x2i[j])';
    end
end



M3 = zeros(N3*N3,N3*N3);
for i = 1:length(x3i)
    for j = 1:length(x3i)
        M3 = M3+w3i[i]*w3i[j]*phi3xy(x3i[i],x3i[j])*phi3xy(x3i[i],x3i[j])';
    end
end

phi1xydx = (x,y)-> (ForwardDiff.derivative(x->phi1xy(x,y),x));
Dx = zeros(N1*N1,N3*N3);
for i = 1:length(x1i)
    for j = 1:length(x1i)
        Dx = Dx+w1i[i]*w1i[j]*phi1xydx(x1i[i],x1i[j])*phi2xy(x1i[i],x1i[j])';
    end
end

phi1xydy = (x,y)-> (ForwardDiff.derivative(y->phi1xy(x,y),y));
Dy = zeros(N1*N1,N3*N3);
for i = 1:length(x1i)
    for j = 1:length(x1i)
        Dy = Dy+w1i[i]*w1i[j]*phi1xydy(x1i[i],x1i[j])*phi2xy(x1i[i],x1i[j])';
    end
end

# boundary input left
phix = phi1;
phiy = phi1;

Bpartial_left = zeros(N1,1);
for i = 1:length(x1i)
    Bpartial_left += w1i[i]*phiy(x1i[i]);
end
Bpartial_down = zeros(N1,1);
for i = 1:length(x1i)
    Bpartial_down += w1i[i]*phix(x1i[i]);
end


# assemblage of matrices
vectorinc = collect(1:(Pe1*Nex+1)*(Pe1*Ney+1));
matrixindex = reshape(vectorinc,Pe1*Nex+1,Pe1*Ney+1)';
elmatrix = (nelx, nely) -> matrixindex[((nely-1)*Pe1+1:Pe1*nely+1), ((nelx-1)*Pe1+1:Pe1*nelx+1)];

if Pe23 == 0
    vectorinc2 = collect(1:(Nex)*(Ney));
    matrixindex2 = reshape(vectorinc2,Nex,Ney)';
    Massmatrix2 = zeros(Nex*Ney,Nex*Ney);
    Massmatrix3 = zeros(Nex*Ney,Nex*Ney);
    Dxmatrix = zeros((Pe1*Nex+1)*(Pe1*Ney+1), Nex*Ney);
    Dymatrix = zeros((Pe1*Nex+1)*(Pe1*Ney+1), Nex*Ney);
    elmatrix2 = (nelx, nely) -> [matrixindex2[nelx, nely]];
else
    vectorinc2 = collect(1:(Pe23*Nex+1)*(Pe23*Ney+1));
    matrixindex2 = reshape(vectorinc2,(Pe23*Nex+1),(Pe23*Ney+1))';
    elmatrix2 = (nelx, nely) -> matrixindex2[((nely-1)*Pe23+1:Pe23*nely+1), ((nelx-1)*Pe23+1:Pe23*nelx+1)];
    
    Massmatrix2 = zeros((Pe23*Nex+1)*(Pe23*Ney+1), (Pe23*Nex+1)*(Pe23*Ney+1));
    Massmatrix3 = zeros((Pe23*Nex+1)*(Pe23*Ney+1),(Pe23*Nex+1)*(Pe23*Ney+1));
    Dxmatrix = zeros((Pe1*Nex+1)*(Pe1*Ney+1), (Pe23*Nex+1)*(Pe23*Ney+1));
    Dymatrix = zeros((Pe1*Nex+1)*(Pe1*Ney+1), (Pe23*Nex+1)*(Pe23*Ney+1));
end


Massmatrix1 = zeros((Pe1*Nex+1)*(Pe1*Ney+1),(Pe1*Nex+1)*(Pe1*Ney+1));

Bpartial_left_full = zeros((Pe1*Nex+1)*(Pe1*Ney+1),Ney);
Bpartial_right_full = zeros((Pe1*Nex+1)*(Pe1*Ney+1),Ney);
Bpartial_up_full = zeros((Pe1*Nex+1)*(Pe1*Ney+1),Nex);
Bpartial_down_full = zeros((Pe1*Nex+1)*(Pe1*Ney+1),Nex);

### assemblage

for ix = 1:Nex
    for iy = 1:Ney
        elMatrix = elmatrix(ix,iy);
        elMatrix = elMatrix[:];
        elMatrix2 = elmatrix2(ix,iy);
        elMatrix2 = elMatrix2[:];
        for i = 1:length(elMatrix)
            for j = 1:length(elMatrix)
            Massmatrix1[elMatrix[i], elMatrix[j]] = Massmatrix1[elMatrix[i], elMatrix[j]] +  M1[i,j];
            end
            for j = 1:length(elMatrix2)
                Dxmatrix[elMatrix[i], elMatrix2[j]] = Dxmatrix[elMatrix[i], elMatrix2[j]] + Dx[i,j];
                Dymatrix[elMatrix[i], elMatrix2[j]] = Dymatrix[elMatrix[i], elMatrix2[j]] + Dy[i,j];
            end
        end
        if Pe23 == 0
            Massmatrix2[matrixindex2[iy,ix],matrixindex2[iy,ix]] = M2[1];
            Massmatrix3[matrixindex2[iy,ix],matrixindex2[iy,ix]] = M3[1];
            else
                 for i = 1:length(elMatrix2)
                     for j = 1:length(elMatrix2)
                         Massmatrix2[elMatrix2[i], elMatrix2[j]] = Massmatrix2[elMatrix2[i], elMatrix2[j]] +  M2[i,j];
                         Massmatrix3[elMatrix2[i], elMatrix2[j]] = Massmatrix3[elMatrix2[i], elMatrix2[j]] +  M3[i,j];
                     end
                 end
                
            end
       if ix == 1
           [iy, ix]
           Bpartial_left_full[elMatrix[1], iy] = Bpartial_left_full[elMatrix[1], iy]+Bpartial_left[1];
           Bpartial_left_full[elMatrix[2], iy] = Bpartial_left_full[elMatrix[2], iy] +Bpartial_left[2];
       end
       if ix == Nex
           [iy, ix]
           Bpartial_right_full[elMatrix[3], iy] = Bpartial_right_full[elMatrix[3], iy]+Bpartial_left[1];
           Bpartial_right_full[elMatrix[4], iy] = Bpartial_right_full[elMatrix[4], iy]+Bpartial_left[2];
       end
       if iy == 1
           Bpartial_up_full[elMatrix[1], ix] = Bpartial_up_full[elMatrix[1], ix]+Bpartial_down[1];
           Bpartial_up_full[elMatrix[3], ix] = Bpartial_up_full[elMatrix[3], ix]+Bpartial_down[1];
       end
       if iy == Ney
           Bpartial_down_full[elMatrix[2], ix] = Bpartial_down_full[elMatrix[2], ix]+Bpartial_down[1];
           Bpartial_down_full[elMatrix[4], ix] = Bpartial_down_full[elMatrix[4], ix]+Bpartial_down[1];
       end
    end
end


## saint-venant equations
N1 = size(Dxmatrix,1);
N2 = size(Dxmatrix,2);
J = [zeros(size(Massmatrix1)) Dxmatrix Dymatrix; -Dxmatrix' zeros(N2,2*N2); -Dymatrix' zeros(N2,2*N2)];
Q = blkdiag(blkdiag(inv(Massmatrix1), inv(Massmatrix2)), inv(Massmatrix3));


Mpp = blkdiag(Massmatrix2,Massmatrix3);
DD = [Dxmatrix Dymatrix];

#a,v = eig(J*Q);



A = J*Q;

xpos = linspace(0,Lx, Nex+1);
ypos = linspace(0,Ly, Ney+1);

#[XX, YY] = meshgrid(xpos, ypos)
XX = collect(xpos) .+ xpos'*0;
YY = collect(ypos)' .+ ypos*0;

x0 = zeros(size(A,1),1);
x0[1:N1] = Massmatrix1 * (1+XX[:]*0.0+YY[:]*0.0 - XX[:].*YY[:]*0.5*0);

B = [Bpartial_up_full*ones(Nex,1) -Bpartial_left_full*ones(Ney,1); zeros(size(A,1)-(Nex+1)*(Ney+1),2)]*0.2;
C = eye(size(A,1));


function Hamiltonianintegrand(alpha1, alpha2, alpha3)
	0.5*(alpha1^2 + alpha1*(alpha2^2+alpha3^2))
end

function DiscreteHamiltonian(X)
	alpha1 = X[1:N1];
    alpha2 = X[N1+1:N1+N2];
    alpha3 = X[N1+N2+1:N1+N2*2];
	Hd = 0;	
    for ix = 1:Nex
        for iy = 1:Ney
            alpha1_el = alpha1[elmatrix(ix,iy)'[:]];
            alpha2_el = alpha2[elmatrix2(ix,iy)'[:]];
            alpha3_el = alpha3[elmatrix2(ix,iy)'[:]];

            alpha1_elxy = (x,y) -> (alpha1_el'*phi1xy(x,y))[1];
            alpha2_elxy = (x,y)-> (alpha2_el'*phi2xy(x,y))[1];
            alpha3_elxy = (x,y)-> (alpha3_el'*phi3xy(x,y))[1];

            for i = 1:length(xi)
                for j = 1:length(xi)
                Hd = Hd + wi[i]*wi[j]'*Hamiltonianintegrand(alpha1_elxy(xi[i],xi[j]), alpha2_elxy(xi[i],xi[j]), alpha3_elxy(xi[i],xi[j]));
                end
            end
        end
	end
	Hd
end

grad = x-> ForwardDiff.gradient(DiscreteHamiltonian, Q*x)
hessian = x-> ForwardDiff.hessian(DiscreteHamiltonian, Q*x)


X = randn(size(Q,1));

println(DiscreteHamiltonian(Q*X)-X'*Q*X/2)


using DifferentialEquations
amp = 2;
function dynamics_linear(xd, x, p, t)
    xd[:] = A*x + amp*B*[1, 1]*(t/(1+t))*sin(pi*t)*(t<1)
    println(t)
end

function dynamicsNL(xd, x, p, t)
    xd[:] = A*grad(x) + amp*B*[1, 1]*(t/(1+t))*sin(pi*t)* (t<1)
end

prob = ODEProblem(dynamicsNL,x0,(0.,2))
sol = solve(prob,Trapezoid());

problin = ODEProblem(dynamics_linear,x0*0,(0.,2))
sollin = solve(problin,Trapezoid());


function plotx(x)
    surf(XX,YY,reshape((Massmatrix1 \ x[1:N1])',  Nex+1,Ney+1)[:,:], cmap = "winter")
end

function plotxlin(x)
    surf(XX,YY,reshape((Massmatrix1 \ x[1:N1])',  Nex+1,Ney+1)[:,:]+1, alpha=0.5, color = "k")
end
times = 0.4:0.4:1.6;
using PyPlot
for i = times
    figure()
    plotx(sol(i));    
    plotxlin(sollin(i));
    ax = gca();
    ax[:view_init](elev = 45, azim = -60)
    ax[:set_zlim]([0.7,1.4])
end
for i = 1:4
    figure(i)
    legend(["nonlinear", "linear"])
    filename = string("simulation2DNL_", Int(times[i]*10));
    savefig(string(filename,".jpg"));
end



for i = times
    figure();
    plotx(sol(i));    
    plotxlin(sollin(i));
    ax = gca();
    ax[:view_init](elev = 45, azim = -60)
    ax[:set_zlim]([0.9995,1.0004])
end
for i = 1:4
    figure(i)
    legend(["nonlinear", "linear"])
    filename = string("simulation2DNL_large_", Int(times[i]*10));
    savefig(string(filename,".jpg"));
end


times = 0:0.1:2;
Hd = zeros(length(times))
Vd = zeros(length(times))
for i = 1:length(times)
    Hd[i] = DiscreteHamiltonian(Q*sol(times[i]))
    Vd[i] = sum(sol(times[i])[1:N1])
end
figure();
plot(times,Hd);
grid();
xlabel("time (s)");
ylabel("Hamiltonian");
savefig(string("Hamiltonian2Dtime",".jpg"));

figure();
plot(times,Vd);
grid();
xlabel("time (s)");
ylabel("Total volume of fluid");
savefig(string("Volume2Dtime",".jpg"));




figure(14);
subplot(2,2,1);
solti = ti-> reshape((Massmatrix1 \ sol(ti)[1:N1])',  Nex+1,Ney+1)[:,:][6,:]
sollinti = ti-> reshape((Massmatrix1 \ sollin(ti)[1:N1])',  Nex+1,Ney+1)[:,:][6,:]
plot(solti(0.4));
plot(sollinti(0.4)+1,  linestyle = "dashed");
xlabel("x position"); ylabel("fluid height");
subplot(2,2,2);
plot(solti(0.8));
plot(sollinti(0.8)+1,  linestyle = "dashed");
xlabel("x position"); ylabel("fluid height");
subplot(2,2,3);
plot(solti(1.2));
plot(sollinti(1.2)+1,  linestyle = "dashed");
xlabel("x position"); ylabel("fluid height");
subplot(2,2,4);
plot(solti(1.6));
plot(sollinti(1.6)+1,  linestyle = "dashed");
xlabel("x position"); ylabel("fluid height");
legend(["nonlinear", "linear"])
tight_layout();
filename = "simulations2d_snapshot_large";
savefig(string(filename,".jpg"));


figure();
spy(Massmatrix1);
savefig("sparsityMq.pdf");

figure();
spy(Mpp);
savefig("sparsityMp.pdf");

figure();
spy(DD);
savefig("sparsityD.pdf");

