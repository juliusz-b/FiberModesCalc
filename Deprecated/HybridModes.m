clear all;clc;close all;
a = 8.3/2;
b = a;
neff = 1.447;
nclad = 1.4440;
ncore = 1.4513;DELTA=(ncore-nclad)/ncore*100;
lambda = 0.85;
%omega = 2*pi*f;
e0 = 8.854187e-12;
m0 = 4*pi*1e-7;
LP = [1,2];


k = 2*pi/(lambda);LS=physconst('lightspeed');
beta_range = [k*nclad k*ncore];
c=LS;
omega = k * LS;
V = a*k*(ncore^2-nclad^2)^0.5;
m=LP(1);l=LP(2);
if m==0
    Xm = @(m,w) (-besselk(1,w))./(w.*besselk(0,w));
    Ym = @(m,u) (-besselj(1,u))./(u.*besselj(0,u));
    besselkDerivative = @(m, x) -besselk(1,x);
    besseljDerivative = @(m, x) -besselj(1,x);
else
    Xm = @(m,w) (-besselk(1,w))./(w.*besselk(0,w));
    Ym = @(m,u) (-besselj(1,u))./(u.*besselj(0,u));
    besselkDerivative = @(m, x) (-0.5*(besselk(m-1,x)+besselk(m+1,x)));
    besseljDerivative = @(m, x) (0.5*(besselj(m-1,x)-besselj(m+1,x)));
end
u = @(beta) a*(k^2*ncore^2-beta.^2).^0.5;
w = @(beta) a*(beta.^2 - k^2*nclad^2).^0.5;
B = @(beta) (w(beta).^2)./(u(beta).^2+w(beta).^2);
equation = @(beta) (Xm(m,w(beta))+Ym(m,u(beta))).*(nclad^2*Xm(m,w(beta))+ncore^2*Ym(m,u(beta)))-((m*beta/k)./((u(beta)).^1.*B(beta).^1)).^2;

%%% mode search
mode_beta=[];
for bet=beta_range(1)+0.01:0.001:beta_range(2)-0.01
    mode_beta((beta_range(1)+0.01:0.001:beta_range(2)-0.01)==bet) = fsolve(equation,[bet], optimoptions('fsolve', 'Display', 'off'));
end
mode_beta(mode_beta<beta_range(1))=[];
mode_beta(mode_beta>beta_range(2))=[];

if isempty(mode_beta)
   error('Cannot found appropriate beta for given parameters'); 
end

mode_beta(diff(mode_beta)<0.0001)=[];
mode_beta = flip(sort(mode_beta));
%%% plotting dispersion equation vs beta
figure();
betas = [beta_range(1):.001:beta_range(2)];
plot(betas, equation(betas), 'linewidth', 3);
hold on;
plot(mode_beta, equation(mode_beta),'o', 'linewidth', 3);
plot(mode_beta(l), equation(mode_beta(l)),'+', 'linewidth', 1.6, 'markersize', 10, 'color', 'black');
ylim([-1 1]);xlim(beta_range);
xlabel('\beta');ylabel('Dispersion equation value');grid on;
set(gca, 'fontsize', 15);
legend({'Disp. eq.', 'Found betas', 'Chosen beta'}, 'location', 'best');
%%% choosing beta
beta = mode_beta(l);

M = zeros(4,4);
M(1,1) = besselj(m,u(beta));
M(1,2) = 0;
M(1,3) = -besselk(m,w(beta));
M(1,4) = 0;

M(2,1) = 0;
M(2,2) = M(1,1);
M(2,3) = 0;
M(2,4) = M(1,3);

M(3,1) = 1i*beta*m/(a*(u(beta)/a)^2)*besselj(m,u(beta));
M(3,2) = -(k*m0)/(u(beta)/a)*( besseljDerivative(m, u(beta)) );
M(3,3) = 1i*beta*m/(a*(w(beta)/a)^2)*besselk(m,w(beta));
M(3,4) = -(k*m0)/(w(beta)/a)*( besselkDerivative(m, w(beta)) );

M(4,1) = (k*e0)/(u(beta)/a)*ncore^2*( besseljDerivative(m, u(beta)) );
M(4,2) = M(3,1);
M(4,3) = (k*e0)/(w(beta)/a)*nclad^2*( besselkDerivative(m, w(beta)) );
M(4,4) = M(3,3);
Y = linsolve(M,[1;0;0;0]);
Be = Y(3);Bh = Y(4);Ae = Y(1);Ah = Y(2);

phase = pi/4;
r_s = linspace(0.0001, b, 250);fi_s = linspace(0.0001, 2*pi, 250);
P = [];Ez = [];
for r=r_s
    
    for fi=fi_s
        
        gamma = w(beta)/a;
        chi = u(beta)/a;
        
        if r>=0 && r<=a
            %Core region
            %E = [Er, Ephi, Ez];%H = [Hr, Hphi, Hz];
            E(1,1) = -1i/chi^2*( Ae*beta*chi*besseljDerivative(m,chi*r) + 1i*omega*m0*m/r*Ah*besselj(m, chi*r) )*cos(m*fi+phase);
            E(1,2) = -1i/chi^2*( Ae*beta*m/r*besselj(m,chi*r) - omega*m0*chi*Ah*besseljDerivative(m, chi*r) )*sin(m*fi+phase);
            E(1,3) = Ae*besselj(m,chi*r)*cos(m*fi+phase);
            
            H(1,1) = -1i/chi^2*( Ah*beta*chi*besseljDerivative(m,chi*r) - 1i*omega*e0*ncore^2*m/r*Ae*besselj(m, chi*r) )*sin(m*fi+phase);
            H(1,2) = -1i/chi^2*( Ah*1i*beta*m/r*chi*besselj(m,chi*r) + omega*e0*ncore^2*chi*Ae*besseljDerivative(m, chi*r) )*cos(m*fi+phase);
            H(1,3) = Ah*besselj(m, chi*r)*sin(m*fi+phase);
        else
            %Cladding region
            %E = [Er, Ephi, Ez];%H = [Hr, Hphi, Hz];
            E(1,1) = -1i/gamma^2*( Be*beta*gamma*besselkDerivative(m,gamma*r) + 1i*omega*m0*m/r*Bh*besselk(m, gamma*r) )*cos(m*fi+phase);
            E(1,2) = -1i/gamma^2*( Be*beta*m/r*besselk(m,gamma*r) - omega*m0*gamma*Bh*besselkDerivative(m, gamma*r) )*sin(m*fi+phase);
            E(1,3) = Be*besselk(m,gamma*r)*cos(m*fi+phase);
            
            H(1,1) = -1i/gamma^2*( Bh*beta*gamma*besselkDerivative(m,gamma*r) - 1i*omega*e0*nclad^2*m/r*Be*besselk(m, gamma*r) )*sin(m*fi+phase);
            H(1,2) = -1i/gamma^2*( Bh*1i*beta*m/r*gamma*besselk(m,gamma*r) + omega*e0*nclad^2*gamma*Be*besselkDerivative(m, gamma*r) )*cos(m*fi+phase);
            H(1,3) = Bh*besselk(m, gamma*r)*sin(m*fi+phase);
        end
        
        
        P(r_s==r, fi_s==fi) = E(1,3);%0.5*(E(1,1)*conj(H(1,2)) - E(1,2)*conj(H(1,1)));
        
    end
    
end

figure();
polarplot3d(abs(P), 'RadialRange',[0 b],'plottype','surfn', 'InterpMethod', 'spline','polargrid',{1 1})
view([0 90])