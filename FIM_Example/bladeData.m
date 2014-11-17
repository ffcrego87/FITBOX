if ~exist('T_init')
    T_init = 0;
end
if ~exist('T_limit')
    T_limit = 10;
end
if ~exist('Ts')
    Ts = 0.01; %sample time [s]
end
if ~exist('fsgn1') || ~exist('fsgn2') || ~exist('fsgn3')
    %set default fault signals
    time = (T_init:Ts:T_limit)';
    fsgn1 = ones(size(time));
    fsgn2 = ones(size(time));
    fsgn3 = zeros(size(time));
end
%Test scenario for SVO class
%clear allx
seed = randi(10^6,[1,2]);
%fault data
Gain_Beta_2_m2=1.2;
%Continuous Time Dynamic System (Nominal)
omega_n=11.11;
xi=0.6;
A = [0 1;-omega_n^2 -2*xi*omega_n];
B = omega_n^2*[0,0,0;1,-0.5,-0.5];
C = [1,0;1,0];
D = [zeros(2,1),eye(2)];
%Discrete Time Dynamic System (Nominal)
[Ad,Bn,Cd,Dn] = ssdata(c2d(ss(A,B,C,D),Ts));
Bd = Bn(:,1);
Ld = Bn(:,2:3);
Dd = Dn(:,1);
Nd = Dn(:,2:3);
nx = size(Ad,1);
ny = size(Cd,1);
nu = size(Bd,2);
nw = size(Ld,2);
%input data
Np_beta=1.5e-3; %Noise power [deg^2/(rad/s)]
u = 0;
wp = 5*sqrt(Np_beta)*ones(nw,1)./sqrt(Ts); %noise upper bound
wm = -5*sqrt(Np_beta)*ones(nw,1)./sqrt(Ts); %noise lower bound
%fault idenfication module
Treset = T_limit/Ts; %reset after 100 iterations
ff = fim(Treset,nx);
hat_x0 = zeros(2,1);
uncertainty = [1;1];
options.approximation = 'hypercube';
options.solver = 'matlab';
options.angTol = 1e-2;
options.numTol = 1e-1;
options.accelerator = 'on';
if ~exist('horizon','var')
    horizon = 1;
end
%nominal SVO
nml= SVO(nx,nu,ny,nw,hat_x0,uncertainty,horizon,options);
classNml = 1; %class of the nominal SVO
ff.addSVO('nominal',classNml,nml); %add nominal SVO to the bank
%global SVO -- constant set
Aglb = eye(nx);
Bglb = zeros(nx,0); 
Lglb = zeros(nx,0);
Cglb = zeros(0,nx);
Dglb = [];
Nglb = [];
glb = SVO(nx,0,0,0,zeros(nx,1),uncertainty*100,horizon,options);
classGlb = 1;
ff.addSVO('global',classGlb,glb);
%Continuous Time Dynamic System (Fault 7)
xi3=0.9;
omega_n3=3.42;
[Apb,Bpb,Cpb,Dpb]=tf2ss([omega_n^2],[1 2*xi*omega_n omega_n^2]);
Bpb = Bpb*sum(Cpb);
Cpb = Cpb/sum(Cpb);
[Apb2,Bpb2,Cpb2,Dpb2]=tf2ss([omega_n3^2],[1 2*xi3*omega_n3 omega_n3^2]);
Bpb2 = Bpb2*sum(Cpb2);
Cpb2 = Cpb2/sum(Cpb2);

A3c = [0 1;-omega_n3^2 -2*xi3*omega_n3];
B3c = omega_n3^2*[0,0,0;1,-0.5,-0.5];
C3c = [1,0;1,0];
D3c = [zeros(2,1),eye(2)];
[A3,Bn,C3,Dn] = ssdata(c2d(ss(A3c,B3c,C3c,D3c),Ts));
B3 = Bn(:,1);
L3 = Bn(:,2:3);
D3 = Dn(:,1);
N3 = Dn(:,2:3);
wp3= wp;
wm3= wm;
%faulty SVO
fault3 = SVO(nx,nu,ny,nw,hat_x0,uncertainty,horizon,options);
fault3.addDelta('A',(Ad-A3)/2);
fault3.addDelta('B',(Bd-B3)/2);
fault3.addDelta('L',(Ld-L3)/2);
classFaulty = 1;
ff.addSVO('fault3',classFaulty,fault3);
%Continuous Time Dynamic System (Fault 2)
A = ([0 1;-omega_n^2 -2*xi*omega_n]-0.5*[0;omega_n^2]*[Gain_Beta_2_m2-1 0]);
B = omega_n^2*[0,0,0;1,-0.5,-0.5*Gain_Beta_2_m2];
C = [1,0;Gain_Beta_2_m2,0];
D = [zeros(2,1),[1,0;0,Gain_Beta_2_m2]];
[A2,Bn,C2,Dn] = ssdata(c2d(ss(A,B,C,D),Ts));
B2 = Bn(:,1);
L2 = Bn(:,2:3);
D2 = Dn(:,1);
N2 = Dn(:,2:3);
wp2= wp;
wm2= wm;
fault2 = SVO(nx,nu,ny,nw,hat_x0,uncertainty,horizon,options);
classFaulty = 2;
ff.addSVO('fault2',classFaulty,fault2);
%Continuous Time Dynamic System (Fault 1)
Constant_Beta_1_m1=5;
A = ([0 1;-omega_n^2 -2*xi*omega_n] - 0.5*[0;omega_n^2]*[1,0]);
B = omega_n^2*[0,0,0;1,-0.5,-0.5];
C = [0,0;1,0];
D = [zeros(2,1),eye(2)];
[A1,Bn,C1,Dn] = ssdata(c2d(ss(A,B,C,D),Ts));
B1 = Bn(:,1);
L1 = Bn(:,2:3);
D1 = Dn(:,1);
N1 = Dn(:,2:3);
wp1= [Constant_Beta_1_m1;wp(2)];
wm1= [Constant_Beta_1_m1;wm(2)];
fault1 = SVO(nx,nu,ny,nw,hat_x0,uncertainty,horizon,options);
classFaulty = 2;
ff.addSVO('fault1',classFaulty,fault1);