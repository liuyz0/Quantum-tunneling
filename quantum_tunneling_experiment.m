%% On quantum simulation for non-convex optimization
% This script runs the 3 numerical experiments presented in the paper
% This is the code for paper
% "﻿On Quantum Speedups for Nonconvex Optimization via
% Quantum Tunneling Walks". Do not distribute.
% The computing environment is MATLAB 2020b.

% we use the finite difference method (FDM) and the 
% leapfrog scheme (LFS) to solve time-dependent Schrodinger equation. 
% The integrator 'leapfrogsolver' (available on MathWorks) was developed 
% by Mauger Franois. 


%% Test 1: 1-dim examples
% including tow parts: QTW and SGD
% see Appendix E.1 for details

%% Initial state illustration

% settings
n = 512; %number of grids in 1-dim
n1 = n + 1;

% % Example 1, the critical case (parameters):
% a = 5; %length of a well
% b = 0.243398682158585; % length of a barrier
% k = 0.014603920929515; % quadratic
% sig = 1; % Gassian: standard deviation
% eps = 0.15; % Gassian: constant

% % Example 2, flatness (parameters):
% a = 5; %length of a well
% b = 1.694890449331924; % length of a barrier
% k = 0.006101605617595; % quadratic
% sig = 1; % Gassian: standard deviation
% eps = 0.009; % Gassian: constant


% % Example 3, sharpness (parameters):
a = 5; %length of a well
b = 1.034579727740753; % length of a barrier
k = 0.014566882566590; % quadratic
sig = 0.5; % Gassian: standard deviation
eps = 0.0088; % Gassian: constant


N = 3; %number of wells
L = 2*(N*a + (N-1)*b); % length of the box
Period = 2*a + 2*b;
dx = L/n1;
X = dx:dx:L-dx;


% set the potential
f0 = zeros(n,1);
for i = 1:n
    for j = 1:N-1
        if (j-1)*Period <= X(i) && X(i) < (j-1)*Period + 2*a
            f0(i) = funcwell(X(i)-(j-1)*Period-a,k);
        end
        if (j-1)*Period + 2*a <= X(i) && X(i) < (j)*Period
            f0(i) = funcbarr(X(i)-(j)*Period+b,sig) +k*a^2/2 - eps;
        end
    end
    if (N-1)*Period <= X(i)
        f0(i) = funcwell(X(i)-(N-1)*Period-a,k);
    end
end
U0 = spdiags(f0,0,n,n);

% build the Laplacian
Laplacian = lap1d(n)./dx^2;

nh = 100;
h = linspace(0.1,1,nh);
Proj = zeros(nh,1);

for i = 1:nh
    % build the Hamiltonian
    H = - Laplacian.*h(i)^2 + U0;
    % eigenstates preparation
    [Gg,Eg] = eigs(H,3,'smallestabs');
    Gg1 =  Gg(:,1)/sqrt(sum(Gg(:,1).^2)*dx);
    Gg2 =  Gg(:,2)/sqrt(sum(Gg(:,2).^2)*dx);
    Gg3 =  Gg(:,3)/sqrt(sum(Gg(:,3).^2)*dx);
    
    % local ground state preparation
    eta = 2*a+2*b; % right boundary of the ground state
    ng = fix(eta./dx);
    Hl = full(H(1:ng,1:ng));
    [Lg,El0] = eigs(Hl,1,'smallestabs');
    Lg = Lg/sqrt(sum(Lg.^2)*dx);
    Lg = [Lg; zeros(n-ng,1)];
    Proj(i) = (sum(Lg.*Gg1)*dx)^2+ (sum(Lg.*Gg2)*dx)^2 + (sum(Lg.*Gg3)*dx)^2;
end

%% plot
nh = 100;
h = linspace(0.1,1,nh);
figure(1)
set(gcf, 'units','points','position',[0 0 500 350]);
load('overlap_1_2a.mat')
plot(h,Proj,'LineWidth',3)
hold on
load('overlap_1_2a1b.mat')
plot(h,Proj,'LineWidth',3)
load('overlap_1_2a2b.mat')
plot(h,Proj,'LineWidth',3)
xlabel('h')
ylabel('Overlap with the subspace')
title('Example 1');
le = legend('Small local region','Middle local region','Large local region',...
    'location','SouthWest');
set(le,'FontName','Times','FontSize',25)
set(gca,'FontSize',25,'FontName','Times')
xlim([0.1,1])
ylim([0.8,1])
%xticks(0.1:0.2:1);


figure(2)
set(gcf, 'units','points','position',[0 0 500 350]);
load('overlap_2_2a.mat');
plot(h,Proj,'LineWidth',3)
hold on
load('overlap_2_2a1b.mat');
plot(h,Proj,'LineWidth',3)
load('overlap_2_2a2b.mat');
plot(h,Proj,'LineWidth',3)
hold off
xlabel('h')
%ylabel('Overlap with the subspace')
title('Example 2');
le = legend('Small local region','Middle local region','Large local region',...
    'location','SouthWest');
set(le,'FontName','Times','FontSize',25)
set(gca,'FontSize',25,'FontName','Times')
xlim([0.1,1])
%xticks(0.1:0.2:1);

figure(3)
set(gcf, 'units','points','position',[0 0 500 350]);
load('overlap_3_2a.mat');
plot(h,Proj,'LineWidth',3)
hold on
load('overlap_3_2a1b.mat');
plot(h,Proj,'LineWidth',3)
load('overlap_3_2a2b.mat');
plot(h,Proj,'LineWidth',3)
hold off
xlabel('h')
%ylabel('Overlap with the subspace')
title('Example 3');
le = legend('Small local region','Middle local region','Large local region',...
    'location','SouthWest');
set(le,'FontName','Times','FontSize',25)
set(gca,'FontSize',25,'FontName','Times')
xlim([0.1,1])
ylim([0.8,1])
%xticks(0.1:0.2:1);
%% Quantum tunneling walks
% settings
n = 512; %number of grids in 1-dim
n1 = n + 1;

% Example 1, the critical case (parameters):
a = 5; %length of a well
b = 0.243398682158585; % length of a barrier
k = 0.014603920929515; % quadratic
sig = 1; % Gassian: standard deviation
eps = 0.15; % Gassian: constant
h = 0.8; % Planck constant which can be varied!


% Example 2, flatness (parameters):
% a = 5; %length of a well
% b = 1.694890449331924; % length of a barrier
% k = 0.006101605617595; % quadratic
% sig = 1; % Gassian: standard deviation
% eps = 0.009; % Gassian: constant
% h = 0.8; % Planck constant which can be varied!

% % Example 3, sharpness (parameters):
% a = 5; %length of a well
% b = 1.034579727740753; % length of a barrier
% k = 0.014566882566590; % quadratic
% sig = 0.5; % Gassian: standard deviation
% eps = 0.0088; % Gassian: constant
% h = 1; % Planck constant which can be varied!


N = 3; %number of wells
L = 2*(N*a + (N-1)*b); % length of the box
Period = 2*a + 2*b;
dx = L/n1;
X = dx:dx:L-dx;


% set the potential
f0 = zeros(n,1);
for i = 1:n
    for j = 1:N-1
        if (j-1)*Period <= X(i) && X(i) < (j-1)*Period + 2*a
            f0(i) = funcwell(X(i)-(j-1)*Period-a,k);
        end
        if (j-1)*Period + 2*a <= X(i) && X(i) < (j)*Period
            f0(i) = funcbarr(X(i)-(j)*Period+b,sig) +k*a^2/2 - eps;
        end
    end
    if (N-1)*Period <= X(i)
        f0(i) = funcwell(X(i)-(N-1)*Period-a,k);
    end
end
U0 = spdiags(f0,0,n,n);

% build the Laplacian
Laplacian = lap1d(n)./dx^2;

% build the Hamiltonian 
H = - Laplacian.*h^2 + U0;

% local ground state preparation
eta = 2*a; % right boundary of the ground state
ng = fix(eta./dx);
Hg = full(H(1:ng,1:ng));
[Eigenfs,Eigenvs]=eig(Hg);
Ground = abs(Eigenfs(:,1))/sqrt(sum(Eigenfs(:,1).^2)*dx);
Ground = [Ground; zeros(n-ng,1)];

%% run the simulation by LFS to see time evolution
tic
dt = 1/Eigenvs(ng,ng);
tspan = linspace(0,1000);
Phi = leapfrogsolver(H,Ground,tspan,dt);
tEnd=toc;
fprintf('Test 1, running time = %d\n',tEnd);
%% plot
% exprisk = zeros(100,1);
% for i = 1:100
%     exprisk(i) = sum(abs(Phi(i,:)).^2.*f0'*dx);
% end

% figure(4)
% set(gcf, 'units','points','position',[0 0 1000 1000]);
% for i = 1:9
%     if mod(i,1)==0
%         subplot(3,3,i)
%         yyaxis left
%         plot(X,abs(Phi(i,:)).^2,'LineWidth',2);
%         ylim([0,0.25])
%         ylabel('Distribution')
%         yyaxis right
%         plot(X,f0,'LineWidth',2)
%         xticks([0,10,20,30]);
%         %yticks([-2,-1,0,1,2]);
%         xlim([0,L])
%         xlabel('X')
%         ylabel('Objective function')
%         title(['t =',num2str((i-1)*36)]);
%     end
% end
% set(gca,...
%             'FontSize',14,'FontName','Times');

mean(exprisk)
std(exprisk)
figure(5)
plot(1:100,exprisk)

%% run the simulation by LFS and then measuring
experiment = 1000;
Thit = [];
Xhit = [];

tic

for ii = 1:experiment
tau = 288;
trial = 1000; % number of trials
totalt = 0;
for i = 1:trial 
    t = rand*tau;
    totalt = t + totalt;
    dt = 1/Eigenvs(ng,ng);
    Phi = leapfrogsolver(H,Ground,[0,t],dt);
    fprintf('Trial = %d',i);
    PDF = abs(Phi(2,:)).^2*dx;
    measurement = rand;
    CDF = 0;
    for j = 1:n
        if CDF < measurement && measurement <= CDF + PDF(j)
            measuredx = j*dx;
            break
        else
            CDF = CDF + PDF(j);
        end
    end
    fprintf('x = %d\n', measuredx);
    if L-2*a <= measuredx
        fprintf('Success!\n');
        Thit = [Thit,totalt];
        Xhit = [Xhit,measuredx];
        break
    end
end
end
totalEnd=toc;
fprintf('total running time = %d\n',totalEnd);

%% plot
mean(Thit)
sqrt(var(Thit))

histogram(Thit,50)

%% Stoachastic gradient descent

% Example 1, the critical case (parameters):
% a = 5; %length of a well
% b = 0.243398682158585; % length of a barrier
% k = 0.014603920929515; % quadratic
% sig = 1; % Gassian: standard deviation
% eps = 0.15; % Gassian: constant
% s = 0.1525; % learning rate which can be varied!


% 
% Example 2, flatness (parameters):
% a = 5; %length of a well
% b = 1.694890449331924; % length of a barrier
% k = 0.006101605617595; % quadratic
% sig = 1; % Gassian: standard deviation
% eps = 0.009; % Gassian: constant
% s = 0.129; % learning rate which can be varied!

% % Example 3, sharpness (parameters):
a = 5; %length of a well
b = 1.034579727740753; % length of a barrier
k = 0.014566882566590; % quadratic
sig = 0.5; % Gassian: standard deviation
eps = 0.0088; % Gassian: constant
s = 0.189; % experiment determined s; do not change!

N = 3; %number of wells
L = 2*(N*a + (N-1)*b); % length of the box
Period = 2*a + 2*b;



% Gibbs distribution 
Z = sum(exp(-2*f0/s)*dx);
Expectedf = sum(f0.*exp(-2*f0/s)*dx)/Z;


%% run SGD
% ns = 10;
% sspan = linspace(0.1,0.5,ns);
% Tsgd = zeros(1000,ns);
% Expecf = zeros(1,ns);
% for nn = 1:ns
% s = sspan(nn);
% 
% % Gibbs distribution
% Z = sum(exp(-2*f0/s)*dx);
% Expecf(nn) = sum(f0.*exp(-2*f0/s)*dx)/Z;

experiment = 1000;
Tsgd = [];

for ii = 1:experiment
    Tupper = 1e+7;
    Stepupper = Tupper/s;
    Xs = a;
    %Trajectory = Xs;
    for i = 1:Stepupper
        for j = 1:N-1
            if (j-1)*Period <= Xs && Xs < (j-1)*Period + 2*a
                gradXs = diffwell(Xs-(j-1)*Period-a,k);
            end
            if (j-1)*Period + 2*a <= Xs && Xs < (j)*Period
                gradXs = diffbarr(Xs-(j)*Period+b,sig);
            end
        end
        if (N-1)*Period <= Xs
            gradXs = diffwell(Xs-(N-1)*Period-a,k);
        end
        Xs1 = Xs -s*gradXs - s*randn;
        if Xs1 < 0 || Xs1 > L
            Xs = Xs+0;
        else
            Xs = Xs1;
        end
        %Trajectory = [Trajectory,Xs];
        if Xs >= L-2*a
            %fprintf('Success!\n');
            thit = i*s;
            Tsgd = [Tsgd,thit];
            %Tsgd(ii,nn) = thit;
            break
        end
    end
end


%end

%% plot
figure(5)
set(gcf, 'units','points','position',[0 10 1000 300])
subplot(1,3,1)
load('Tsgd1000Exp1.mat')
histogram(Tsgd/10,'BinWidth',50)
hold on 
load('Thit1000Exp1.mat')
histogram(Thit,'BinWidth',50)
ylabel('Frequency')
xlabel('Hitting time')
ylim([0,250])
legend('T_{hit}^{SGD}/10','T_{hit}^{QTW}')
title('Example 1')
set(gca,'FontSize',15,'FontName','Times','units','points','position',[55 45 285 235])

subplot(1,3,2)
load('Tsgd1000Exp2.mat')
histogram(Tsgd/10,'BinWidth',250)
hold on
load('Thit1000Exp2.mat')
histogram(Thit,'BinWidth',250)
xlabel('Hitting time')
legend('T_{hit}^{SGD}/10','T_{hit}^{QTW}')
title('Example 2')
set(gca,'FontSize',15,'FontName','Times','units','points','position',[380 45 285 235])

subplot(1,3,3)
load('Tsgd1000Exp3.mat')
histogram(Tsgd/10,'BinWidth',2000)
hold on
load('Thit1000Exp3.mat')
histogram(Thit,'BinWidth',2000)
xlabel('Hitting time')
legend('T_{hit}^{SGD}/10','T_{hit}^{QTW}')
title('Example 3')
set(gca,'FontSize',15,'FontName','Times','units','points','position',[710 45 285 235])

% mean(Thit)
% mean(Tsgd)


%Tsgdmean = mean(Tsgd);
% ypos = log10(max(Tsgd)) - log10(Tsgdmean);
% yneg = log10(Tsgdmean) - log10(min(Tsgd));
% 
% 
% figure(5)
% set(gcf, 'units','points','position',[0 0 600 450]);
% yyaxis left
% errorbar(sspan(3:10),log10(Tsgdmean(3:10)),yneg(3:10),ypos(3:10),'-o','LineWidth',2,'markersize',10)
% xlabel('s')
% ylabel('log_{10} (T^{SGD}_{hit})')
% yyaxis right
% plot(sspan,Expecf,'-s','LineWidth',2,'markersize',10)
% ylabel('Expected risk')
% title('Example 3')
% set(gca,'looseInset',[0 0 0 0],'FontSize',25,'FontName','Times') %去掉多余白边

%% Test 2: 2-dim example
% including tow parts: QTW and SGD
% we use the constructed landscape (dimension = 2) in Appendix E.3

%% initial state preparation
n = 399; %number of grids in 1-dim
n1 = n + 1;
a = 1;
L = 4*a;
dx = 2*L/n1;
xaxis = -L+dx:dx:L-dx;
yaxis = xaxis;
[X,Y]=meshgrid(xaxis,yaxis);

omega = 0.5;
b = 1.4;
H1 = omega^2*a^2/2;
H2 = 20*H1;

% objective function
f0 = ones(n,n)*H2;
for i = 1:n
    for j = 1:n
        if X(i,j) <= 2*b - a/sqrt(2) && a/sqrt(2) <= X(i,j) && -a/sqrt(2) <= Y(i,j) && Y(i,j)<=a/sqrt(2)
            f0(i,j) = H1;
        end
        if X(i,j)^2 + Y(i,j)^2 <= a^2
            f0(i,j) = funcWell(X(i,j),Y(i,j),omega);
        end
        if (X(i,j)-2*b)^2 + Y(i,j)^2 <= a^2
            f0(i,j) = funcWell(X(i,j)-2*b,Y(i,j),omega);
        end
    end
end
v0 = f0(:);
V0 = spdiags(v0,0,n^2,n^2);


% build the global Laplacian
Laplacian = lap2d(n)./dx^2;

%%
nh = 100;
hspan = linspace(0.05,0.3,nh);
overl = zeros(nh,1);
for ii = 1:nh
    h = hspan(ii);
    % build the global Hamiltonian
    H = - Laplacian.*h^2 + V0;
    
    % local Hamiltonian for the local ground state
    vecX = X(:);
    vecY = Y(:);
    pointer = [];
    for i = 1:n^2
        if vecX(i)<= a && -a <= vecX(i) && -a <= vecY(i) && vecY(i)<= a
            pointer = [pointer,i];
        end
    end
    nl = sqrt(length(pointer));
    vl = v0(pointer);
    Vl = spdiags(vl,0,nl^2,nl^2);
    % build the global Laplacian
    Ll = lap2d(nl)./dx^2;
    Hlocal = - Ll.*h^2 + Vl;
    [Phil,~] = eigs(Hlocal,1,'smallestabs');
    
    Phi0 = zeros(n^2,1);
    for i = 1:nl^2
        Phi0(pointer(i)) = Phil(i);
    end
    Phi0 = Phi0/sqrt(sum(abs(Phi0).^2*dx^2));
    
    % Global ground state and first excited state
    [Phigl,E00] = eigs(H,2,'smallestabs');
    Gg1 =  Phigl(:,1)/sqrt(sum(Phigl(:,1).^2)*dx^2);
    Gg2 =  Phigl(:,2)/sqrt(sum(Phigl(:,2).^2)*dx^2);
    overl(ii) = sum(Phi0.*Gg1*dx^2)^2 + sum(Phi0.*Gg2*dx^2)^2;
end

%% plot
figure()
set(gcf, 'units','points','position',[0 0 500 350]);
plot(hspan,overl,'LineWidth',3)
xlabel('h')
ylabel('Overlap with the subspace')
xlim([0.1,0.3])
ylim([0.8,1])
set(gca,'FontSize',25,'FontName','Times')

%% Quantum tunneling walks
n = 399; %number of grids in 1-dim
n1 = n + 1;
a = 1;
L = 4*a;
dx = 2*L/n1;
xaxis = -L+dx:dx:L-dx;
yaxis = xaxis;
[X,Y]=meshgrid(xaxis,yaxis);

omega = 0.5;
b = 1.4;
H1 = omega^2*a^2/2;
H2 = 20*H1;
h = 0.26; %0.2 is the special one!
S0 = omega*a^2/sqrt(2) + 2*(b-a)*sqrt(H1);

% build the global Laplacian
Laplacian = lap2d(n)./dx^2;

% build the global Hamiltonian 
H = - Laplacian.*h^2 + V0;

% local Hamiltonian for the local ground state
vecX = X(:);
vecY = Y(:);
pointer = [];
for i = 1:n^2
    if vecX(i)<= a && -a <= vecX(i) && -a <= vecY(i) && vecY(i)<= a
        pointer = [pointer,i];
    end
end
nl = sqrt(length(pointer));
vl = v0(pointer);
Vl = spdiags(vl,0,nl^2,nl^2);
% build the global Laplacian
Ll = lap2d(nl)./dx^2;
Hlocal = - Ll.*h^2 + Vl;
[Phil,E0] = eigs(Hlocal,1,'smallestabs');

Phi0 = zeros(n^2,1);
for i = 1:nl^2
    Phi0(pointer(i)) = Phil(i);
end
Phi0 = Phi0/sqrt(sum(abs(Phi0).^2*dx^2));

%sum(v0.*Phi0.^2*dx^2)% accuracy

%% run the simulation by LFS to see time evolution
tic
dt = 0.2*dx^2/h^2;
tspan = 0:25:200;
Phi = leapfrogsolver(H,Phi0,tspan,dt);
tEnd=toc;
fprintf('Test 2, running time = %d\n',tEnd);

sum(abs(Phi(2,:)).^2*dx^2)


%% run the simulation by LFS to find period
pointer = [];
for i = 1:n^2
    if (vecX(i)-2*b)^2 + vecY(i)^2 <= a^2
        pointer = [pointer,i];
    end
end

dt = 0.2*dx^2/h^2;
tspan = [0 5];
Phi = leapfrogsolver(H,Phi0,tspan,dt);

iter = 100;
for i = 1:iter
    tic
    dt = 0.2*dx^2/h^2;
    Phi = leapfrogsolver(H,Phi(2,:),[0 5],dt);
    tEnd=toc;
    fprintf('Test find period, running time = %d\n',tEnd);
    phit = sum(abs(Phi(2,pointer)).^2*dx^2);
    totalp = sum(abs(Phi(2,:)).^2*dx^2);
    if 1.1 <= totalp
        fprintf('fail');
        break
    end
    if 0.90 <= phit
        fprintf('h = %d, phit = %d, half period = %d\n',h,phit,(i+1)*5);
        break
    end
end

%% plot
% figure(3)
% surf(X,Y,reshape(Phigl(:,1)-Phigl(:,2),[n,n]),'FaceAlpha',0.85,'EdgeColor','none')


% exprisk = zeros(100,1);
% for i = 1:100
%     exprisk(i) = sum(abs(Phi(i,:)).^2.*v0'*dx^2);
% end
% mean(exprisk)
% 
% figure(6)
% set(gcf, 'units','points','position',[0 0 600 450]);
% 
% plot(tspan,exprisk,'LineWidth',2)
% ylim([0.03,0.09])
% hold on
% histogram(exprisk,10,'Orientation','horizontal')
% ylim([0.03,0.09])
% xlabel('Time')
% ylabel('Expected risk')
% set(gca,'looseInset',[0 0 0 0],'FontSize',25,'FontName','Times') %去掉多余白边

% figure(7)
% % imagesc(xaxis,yaxis,abs(reshape(Phi(2,:),[n,n])).^2,'AlphaData',.75,[0 1]);
% % hold on
% % contour(X,Y,f0,0:0.025:0.125,'w-','ShowText','on');
% set(gcf, 'units','points','position',[0 0 1000 1000]);
% for i = 1:9
%     subplot(3,3,i)
%     set(gca,...
%         'FontSize',12,'FontName','Times');
%     surf(X,Y,abs(reshape(Phi(i,:),[n,n])).^2+1,'FaceAlpha',0.85,'EdgeColor','none');
%     caxis([1,1.5]);
%     hold on
%     %caxis('manual');
%     contourf(X,Y,f0*1.2);
%     xlabel('X')
%     ylabel('Y')
%     %zlabel('Distribution')
%     title(['t =',num2str((i-1)*25)]);
% end


hspan = 0.16:0.02:0.28;
Tspan = [495,285,205,145,125,85,70];
Tmin = Tspan - 5;
Tmax = Tspan + 5;
yneg = log(Tspan) - log(Tmin);
ypos =  -log(Tspan) + log(Tmax);
hs = linspace(0.16,0.28);
lnT = log(pi/2) - log(2.*hs./pi)/2 - log(H1*2*b*omega^2/4/sqrt(H1))/2 + S0./hs - 2*omega*(b-a)/sqrt(2*H1) + 4*log(b/a);

figure(8)
set(gcf, 'units','points','position',[0 0 600 450]);
errorbar(hspan,log(Tspan),yneg,ypos,'o','LineWidth',2,'markersize',10)
hold on 
plot(hs, lnT,'LineWidth',2)
hold off
xlabel('h')
ylabel('ln(T_{half})')
xlim([0.16,0.28])
legend('Experimental', 'Theoretical')
set(gca,'looseInset',[0 0 0 0],'FontSize',25,'FontName','Times') %去掉多余白边

%% run the simulation by LFS and then measuring
experiment = 100;
Thit = zeros(experiment,1);
%Xhit = zeros(experiment,1);

for ii = 1:experiment
tau = 200;
trial = 1000; % number of trials
totalt = 0;
for i = 1:trial 
    t = rand*tau;
    totalt = t + totalt;
    tic
    dt = 0.22*dx^2/h^2;
    Phi = leapfrogsolver(H,Phi0,[0,t],dt);
    tEnd=toc;
    fprintf('Trial = %d, running time = %d\n',i,tEnd);
    PDF = abs(Phi(2,:)).^2*dx^2;
    measurement = rand;
    CDF = 0;
    for j = 1:n^2
        if CDF < measurement && measurement <= CDF + PDF(j)
            measuredx = [vecX(j),vecY(j)];
            break
        else
            CDF = CDF + PDF(j);
        end
    end
    fprintf('x = %d, y = %d\n', measuredx(1),measuredx(2));
    if (measuredx(1)-2*b)^2 + measuredx(2)^2 <= a^2
        fprintf('Success!\n');
        Thit(ii) = totalt;
        %Xhit = measuredx;
        break
    end
end
end

%% plot
mean(Thit(:))
histogram(Thit,10)


%% Stoachastic gradient descent

a = 1;
omega = 0.5;
b = 1.4;

s = 0.37;% learning rate unchange!! 0.37
experiment = 500;
Tsgd = [];

for ii = 1:experiment
    Tupper = 1e+7;
    Stepupper = Tupper/s;
    Xs = [0,0];
    %Trajectory = Xs;
    for i = 1:Stepupper
        gradXs = [0,0];
        if Xs(1)^2 + Xs(2)^2 <= a^2
            gradXs = [Xs(1)*omega^2, Xs(2)*omega^2];
        end
        if (Xs(1)-2*b)^2 + Xs(2)^2 <= a^2
            gradXs = [(Xs(1)-2*b)*omega^2,Xs(2)*omega^2];
        end
        if Xs(1)^2 + Xs(2)^2 < 32*a^2
            Xs1 = Xs -s*gradXs - s*randn(1,2);
        end
        if Xs1(1)^2 + Xs1(2)^2 < 32*a^2
            Xs = Xs1;
        else
            Xs = Xs+0;
        end
        %Trajectory = [Trajectory,Xs];
        if (Xs(1)-2*b)^2 + Xs(2)^2 <= a^2
            %fprintf('Success!\n');
            thit = i*s;
            Tsgd = [Tsgd,thit];
            break
        end
    end
end

%% plot
figure(4)
plot(1:experiment, Tsgd)

mean(Tsgd(:))

%% Test 3: dimension dependence
% see Appendix F.2 for details

a = 1;
omega = 0.5;
b = 1.4;
sspan = [0.1 0.2 0.3 0.4 0.5 0.6]*a;
ns = 6;
nd = 18;
Qmeans = zeros(ns,nd);


for iiii = 1:ns
s = sspan(iiii);
dspan = zeros(nd,1);
for i = 1:nd
    dspan(i) = 10 + (i-1)*5;
end

experiment = 1000;
Qsuc = zeros(experiment,nd);
tic
for iii = 1:nd
d = dspan(iii);
for ii = 1:experiment
    Tupper = 1e+7;
    Xs = zeros(d,1);
    Stepupper = fix(Tupper/s);
    for i = 1:Stepupper
        if norm(Xs) <= a
            gradXs = Xs*omega^2;
        else
            gradXs = 0;
        end
        Xs1 = Xs -s*gradXs - s*randn(d,1);
        if norm(Xs1) <= 4*sqrt(2)*a
            Xs = Xs1;
        else
            Xs = Xs + 0;
        end
        if abs(Xs(1)) > a/2
            %fprintf('Exp = %d, Time=%d, Success!\n',ii,i*s);
            Qsuc(ii,iii) = i;
            break
        end
    end
end

end
tEnd = toc;
fprintf('DimDepen, Runtime=%d \n',tEnd);
%Qmean = mean(Qsuc);
Qmeans(iiii,:) = mean(Qsuc);
end
%% plot
% figure()
% set(gcf, 'units','points','position',[0 0 600 450])
% for i = 1:ns
%     plot(10:5:95, log(Qmeans(i,:)),'-o','LineWidth',2)
%     hold on
% end
% plot(10:5:95,log(min(Qmeans))/1.05,'--k','LineWidth',2)
% line([30,90],[1,1],'LineWidth',2,'Color','black')
% line([90,90],[1,1+60/256],'LineWidth',2,'Color','black')
% line([90,30],[1+60/256,1],'LineWidth',2,'Color','black')
% text(60,1.8,'d/256','FontSize',20,'FontName','Times')
% hold off
% legend('s = 0.1','s = 0.2','s = 0.3', 's = 0.4', 's = 0.5', 's = 0.6','')
% xlim([10,95])
% ylim([0,15])
% xlabel('d')
% ylabel('ln T_{esc}')
% set(gca,'looseInset',[0 0 0.01 0],'FontSize',25,'FontName','Times') 

figure()
subplot(1,3,1)
load('Qsuc_s1_d10_55.mat')
g = 10:5:55;
boxplot(log(Qsuc),g)
ylim([-1,16])
xlabel('d')
ylabel('lnQ_{esc}')
title('s = 1')
set(gcf, 'units','points','position',[0 10 1000 300])
set(gca,'FontSize',15,'FontName','Times','units','points','position',[50 40 290 240])

subplot(1,3,2)
load('Qsuc_s0.5_d10_100.mat')
g = 10:10:100;
boxplot(log(Qsuc),g)
ylim([-1,16])
xlabel('d')
%ylabel('lnT_{esc}')
title('s = 0.5')
set(gca,'FontSize',15,'FontName','Times','units','points','position',[375 40 290 240])

subplot(1,3,3)
load('Qsuc_s0.25_d10_190.mat')
g = 10:20:190;
boxplot(log(Qsuc),g)
ylim([-1,16])
xlabel('d')
%ylabel('lnT_{esc}')
title('s = 0.25')
set(gca,'FontSize',15,'FontName','Times','units','points','position',[700 40 290 240])


%% Definitions of functions

%functions of one well and one barrier
function y = funcwell(x,k)
y = k.*x.^2./2;
end

function y = funcbarr(x,sig)
y = 1/2/pi/sig^2.*exp(-x.^2./sig^2);
end

function y = diffwell(x,k)
y = k.*x;
end

function y = diffbarr(x,sig)
y = -x./pi/sig^4.*exp(-x.^2./sig^2);
end

function L = lap1d(n) % Laplacian on 1-d grid
e = ones(n,1);
L = spdiags([e -2*e e],-1:1,n,n);
end

function y = funcWell(x,y,omega)
y = (x.^2 + y.^2)*omega^2./2;
end


function L = lap2d(n) % Laplacian on 2-d grid
e = ones(n,1);
B = spdiags([e -4*e e],-1:1,n,n);
A = spdiags([e 0*e e],-1:1,n,n);
I = speye(n);
L = kron(I,B) + kron(A,I);
end


function [v1,lambda,iter] = Rayleigh(A)
n = length(A);
I = speye(n);
v0 = ones(n,1)/sqrt(n);
u0 = A\v0;
mu0 = dot(u0,v0)/dot(v0,v0);
y0 = (A-I*mu0^(-1))\v0;
v1 = y0/norm(y0);
iter = 1;
while norm(v1-v0) >= 1e-6
    v0 = v1;
    u0 = A\v0;
    mu0 = dot(u0,v0)/dot(v0,v0);
    lambda = mu0^(-1);
    y0 = (A-I*mu0^(-1))\v0;
    v1 = y0/norm(y0);
    iter = iter + 1;
    if  iter >= 1000
        disp('Warning: Rayleigh unable to converge!');
        return
    end
end
end

function y = leapfrogsolver(H,y0,tspan,dt)
    % Mauger Franois (2020). Symplectic Leap Frog Scheme (https://www.mathworks.com/matlabcentral/fileexchange/38652-symplectic-leap-frog-scheme), MATLAB Central File Exchange. Retrieved May 28, 2020.
    %N = length(y0);
    Q0 = real(y0);
    P0 = imag(y0);
    dV = @(t,Q) H*Q;
    dT = @(P) H*P;

    [~,Q1,P1] = sympleapfrog(dV, dT, tspan, Q0, P0, dt);
    y = Q1+1i.*P1;
end
