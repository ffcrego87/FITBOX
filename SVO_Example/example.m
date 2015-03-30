set(0,'defaulttextinterpreter','latex')
set(0,'DefaultAxesFontSize',10)
%Test scenario for SVO class
%clear all
Ts = 0.01; %sample time [s]
T  = 100*Ts; %test length [s]
s = RandStream('mt19937ar','Seed',1);
RandStream.setGlobalStream(s);
%Continuous Time Dynamic System (Nominal)
omega_n=11.11;
xi=0.6;
A = [0 1;-omega_n^2 -2*xi*omega_n];
B = omega_n^2*[0,0,0;1,-0.5,-0.5];
C = [1,0;1,0];
D = [zeros(2,1),eye(2)];
[An,B,Cn,D] = ssdata(c2d(ss(A,B,C,D),Ts));
Bn = B(:,1);
Ln = B(:,2:3);
Dn = D(:,1);
Nn = D(:,2:3);
%noise vector
Np_beta=1.5e-3; %Noise power [deg^2/(rad/s)]
w = randn(2,T/Ts)*sqrt(Np_beta)/sqrt(Ts);
%other sim data
nx = size(An,1);
ny = size(Cn,1);
nu = size(Bn,2);
nw = size(Ln,2);
x = zeros(nx,T/Ts);
u = ones(nu,T/Ts);
y = zeros(ny,T/Ts);
%SVO
uncertainty = [1;1];
options.approximation = 'hypercube';
options.solver = 'matlab';
options.angTol = 1e-2;
options.numTol = 1e-1;
options.accelerator = 'off';
horizon = 1;
nml= SVO(nx,nu,ny,nw,x(:,1),uncertainty,horizon,options);
xplus= zeros(nx,T/Ts);
xminus= zeros(nx,T/Ts);
wplus = 5*sqrt(Np_beta)*ones(nw,1)./sqrt(Ts);
wminus= -5*sqrt(Np_beta)*ones(nw,1)./sqrt(Ts);
for ii = 1:T/Ts-1
    y(:,ii) = Cn*x(:,ii)+Dn*u(:,ii)+Nn*w(:,ii);
    nml.update(An,Bn,Ln,Cn,Dn,Nn,u(:,ii),y(:,ii),...
        wplus,wminus);
    [xhat,d] = nml.getStateEstimate();
    xplus(:,ii) = xhat+d/2;
    xminus(:,ii) = xhat-d/2;
    x(:,ii+1) = An*x(:,ii)+Bn*u(:,ii)+Ln*w(:,ii);
end
%%
figure('units','centimeters','position',[0 0 12 7])
h = create_axis([1 0.1 0.35 0.1],[1 0.2 0.1 0.1]);
plot(1:T/Ts-1,[x(1,1:end-1)' xminus(1,1:end-1)' xplus(1,1:end-1)' y(:,1:end-1)'])
grid on
xlabel('Time [s]')
ylabel('Blade Pitch [deg]')
legend({'$x(kT_s)$','$x^-(kT_s)$','$x^+(kT_s)$','$y_1(kT_s)$','$y_2(kT_s)$'},'position',[0.7, 0.2 0.25 0.7])
matlabfrag('../../blade')