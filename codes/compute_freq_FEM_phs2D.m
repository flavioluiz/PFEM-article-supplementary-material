function [ freqcomp, DOF ] = compute_freq_FEM_phs2D( Nex, Ney, Pe1, Pe23 )
    %to do: the element matrices seem OK, now we need to generalize the way the
    %assemblage is made!


    % domain dimensions
    Lx = 1;
    Ly = 1;


    % element size
    dx = Lx / Nex;
    dy = Ly /Ney;
    % basis functions
    syms x y
    phix = [(dx/2-x)/dx; (x+dx/2)/dx]; % first order Lag. polynomial:

    phiy = [(y+dy/2)/dy; (dy/2-y)/dy];

    phixy1 = phix * phiy.';
    phixy1 = phixy1.';
    phixy1 = phixy1(:);

    phix2 = [(x-0)*(x-dx/2)/(-dx/2 - 0)/(-dx/2 -dx/2) ; (x+dx/2)*(x-dx/2)/(x + dx/2)/(x-dx/2); (x-0)*(x+dx/2)/(dx/2 - 0)/(dx/2 +dx/2)]; % second order Lag. polynomial:
    phiy2 = [(y-0)*(y-dy/2)/(-dy/2 - 0)/(-dy/2 -dy/2) ; (y+dy/2)*(y-dy/2)/(y + dy/2)/(y-dy/2); (y-0)*(y+dy/2)/(dy/2 - 0)/(dy/2 +dy/2)]; % second order Lag. polynomial:
    phixy2 = phix2 * phiy2.';
    phixy2 = phixy2.';
    phixy2 = phixy2(:);


    phixy0 = 1*(x/x);

    if Pe1 == 1
        phi1 = phixy1;
        elseif Pe1 == 2
            phi1 = phixy2;
    else
        display("invalid option for polynomial approximation of space 1");
    end

    if Pe23 == 0
        phi2 = phixy0;
        phi3 = phixy0;
        elseif Pe23 == 1;
        phi2 = phixy1;
        phi3 = phixy1;
        elseif Pe23 == 2;
        phi2 = phixy2;
        phi3 = phixy2;
        else
        display("invalid option for polynomial approximation of space 23");
    end

    % element matrices
    M1 = int(int(phi1 * phi1.',-dx/2,dx/2),-dy/2,dy/2);
    M2 = int(int((phi2 * phi2.'),-dx/2,dx/2),-dy/2,dy/2);
    M3 = int(int((phi3 * phi3.'),-dx/2,dx/2),-dy/2,dy/2);
    Dx = int(int( diff(phi1,x) * phi2.',-dy/2,dy/2),-dx/2,dx/2);
    Dy = int(int( diff(phi1,y) * phi3.',-dx/2,dx/2),-dy/2,dy/2);
    %%
    % boundary input left
    Bpartial_left = int(phiy, -dy/2,dy/2);
    Bpartial_down = int(phix, -dx/2,dx/2);

    % assemblage of matrices
    vectorinc = [1:(Pe1*Nex+1)*(Pe1*Ney+1)];
    matrixindex = reshape(vectorinc,Pe1*Nex+1,Pe1*Ney+1)';
    elmatrix = @(nelx, nely) matrixindex(((nely-1)*Pe1+1:Pe1*nely+1), ((nelx-1)*Pe1+1:Pe1*nelx+1));
    %%
    if Pe23 == 0
        vectorinc2 = [1:(Nex)*(Ney)];
        matrixindex2 = reshape(vectorinc2,Nex,Ney)';
        Massmatrix2 = zeros(Nex*Ney);
        Massmatrix3 = zeros(Nex*Ney);
        Dxmatrix = zeros((Pe1*Nex+1)*(Pe1*Ney+1), Nex*Ney);
        Dymatrix = zeros((Pe1*Nex+1)*(Pe1*Ney+1), Nex*Ney);
        elmatrix2 = @(nelx, nely) matrixindex2(nelx, nely);
    else
        vectorinc2 = [1:(Pe23*Nex+1)*(Pe23*Ney+1)];
        matrixindex2 = reshape(vectorinc2,(Pe23*Nex+1),(Pe23*Ney+1))';
        elmatrix2 = @(nelx, nely) matrixindex2(((nely-1)*Pe23+1:Pe23*nely+1), ((nelx-1)*Pe23+1:Pe23*nelx+1));

        Massmatrix2 = zeros((Pe23*Nex+1)*(Pe23*Ney+1));
        Massmatrix3 = zeros((Pe23*Nex+1)*(Pe23*Ney+1));
        Dxmatrix = zeros((Pe1*Nex+1)*(Pe1*Ney+1), (Pe23*Nex+1)*(Pe23*Ney+1));
        Dymatrix = zeros((Pe1*Nex+1)*(Pe1*Ney+1), (Pe23*Nex+1)*(Pe23*Ney+1));
    end

    %%
    Massmatrix1 = zeros((Pe1*Nex+1)*(Pe1*Ney+1));

    Bpartial_left_full = zeros((Pe1*Nex+1)*(Pe1*Ney+1),Ney);
    Bpartial_right_full = zeros((Pe1*Nex+1)*(Pe1*Ney+1),Ney);
    Bpartial_up_full = zeros((Pe1*Nex+1)*(Pe1*Ney+1),Nex);
    Bpartial_down_full = zeros((Pe1*Nex+1)*(Pe1*Ney+1),Nex);
    %
    for ix = 1:Nex
        for iy = 1:Ney
            elMatrix = elmatrix(ix,iy);
            elMatrix = elMatrix(:);
            elMatrix2 = elmatrix2(ix,iy);
            elMatrix2 = elMatrix2(:);
            for i = 1:length(elMatrix)
                for j = 1:length(elMatrix)
                Massmatrix1(elMatrix(i), elMatrix(j)) = Massmatrix1(elMatrix(i), elMatrix(j)) +  M1(i,j);
                end
                for j = 1:length(elMatrix2)
                    Dxmatrix(elMatrix(i), elMatrix2(j)) = Dxmatrix(elMatrix(i), elMatrix2(j)) + Dx(i,j);
                    Dymatrix(elMatrix(i), elMatrix2(j)) = Dymatrix(elMatrix(i), elMatrix2(j)) + Dy(i,j);
                end
            end
           if Pe23 == 0
           Massmatrix2(matrixindex2(iy,ix),matrixindex2(iy,ix)) = M2;
           Massmatrix3(matrixindex2(iy,ix),matrixindex2(iy,ix)) = M3;
           else
                for i = 1:length(elMatrix2)
                    for j = 1:length(elMatrix2)
                        Massmatrix2(elMatrix2(i), elMatrix2(j)) = Massmatrix2(elMatrix2(i), elMatrix2(j)) +  M2(i,j);
                        Massmatrix3(elMatrix2(i), elMatrix2(j)) = Massmatrix3(elMatrix2(i), elMatrix2(j)) +  M3(i,j);
                    end
                end

           end
           if ix == 1
               [iy, ix];
               Bpartial_left_full(elMatrix(1), iy) = Bpartial_left_full(elMatrix(1), iy)+Bpartial_left(1);
               Bpartial_left_full(elMatrix(2), iy) = Bpartial_left_full(elMatrix(2), iy )+Bpartial_left(2);
           end
           if ix == Nex
               [iy, ix];
               Bpartial_right_full(elMatrix(3), iy) = Bpartial_right_full(elMatrix(3), iy)+Bpartial_left(1);
               Bpartial_right_full(elMatrix(4), iy) = Bpartial_right_full(elMatrix(4), iy)+Bpartial_left(2);
           end
           if iy == 1
               Bpartial_up_full(elMatrix(1), ix) = Bpartial_up_full(elMatrix(1), ix)+Bpartial_down(1);
               Bpartial_up_full(elMatrix(3), ix) = Bpartial_up_full(elMatrix(3), ix)+Bpartial_down(1);
           end
           if iy == Ney
               Bpartial_down_full(elMatrix(2), ix) = Bpartial_down_full(elMatrix(2), ix)+Bpartial_down(1);
               Bpartial_down_full(elMatrix(4), ix) = Bpartial_down_full(elMatrix(4), ix)+Bpartial_down(1);
           end
        end
    end

    %% saint-venant equations
    N1 = size(Dxmatrix,1);
    N2 = size(Dxmatrix,2);
    J = [zeros(size(Massmatrix1)) Dxmatrix Dymatrix; -Dxmatrix' zeros(N2,2*N2); -Dymatrix' zeros(N2,2*N2)];
    Q = blkdiag(inv(Massmatrix1), inv(Massmatrix2), inv(Massmatrix3));
    %
    [v,a] = eig(J*Q);
    freq = imag(diag(a));
    positivefreqs = freq>.001;
    freq = freq(positivefreqs);
    [freqord freqindex] = sort(freq)
    vpositive = v(:, positivefreqs);
    %
    modalshapes = vpositive(1:N1,:);
    Nfreq  =  2;
    modalshapen = Massmatrix1 \ modalshapes(:,freqindex(Nfreq));
    freqcomp = freq(freqindex(Nfreq))
    DOF = length(J);
end

