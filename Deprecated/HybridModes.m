clear all;clc;close all;
a = 8.6e-6/2;
neff = 1.447;
nclad = 1.4440;
ncore = 1.4513;DELTA=(ncore-nclad)/ncore*100;
lambda = 650e-9;
%omega = 2*pi*f;
r_s = [0:.001:125]*1e-6; %wektor wszystkich srednic
fi_s = [0:.1:2*pi]; %wektor katow
e0 = 8.854187e-12;
m0 = 4*pi*1e-7;

k = 2*pi/(lambda);LS=physconst('lightspeed');
c=LS;
omega = k * LS;
V = a*k*(ncore^2-nclad^2)^0.5;
m=1;
if m==0
    Xm = @(m,w) (-besselk(1,w))./(w.*besselk(0,w));
    Ym = @(m,u) (-besselj(1,u))./(u.*besselj(0,u));
    besselkDerivative = @(m, x) -besselk(1,x);
    besseljDerivative = @(m, x) -besselj(1,x);
else
    Xm = @(m,w) (-0.5*(besselk(m-1,w)+besselk(m+1,w)))./(w.*besselk(m,w));
    Ym = @(m,u) (0.5*(besselj(m-1,u)-besselj(m+1,u)))./(u.*besselj(m,u));
    besselkDerivative = @(m, x) (-0.5*(besselk(m-1,x)+besselk(m+1,x)));
    besseljDerivative = @(m, x) (0.5*(besselj(m-1,x)-besselj(m+1,x)));
end
u = @(beta) a*(k^2*ncore^2-beta.^2).^0.5;
w = @(beta) a*(beta.^2 - k^2*nclad^2).^0.5;
B = @(beta) (w(beta).^2)./(u(beta).^2+w(beta).^2);
equation = @(beta) real(Xm(m,w(beta))+Ym(m,u(beta))).*(nclad^2*Xm(m,w(beta))+ncore^2*Ym(m,u(beta)))-((m*beta/k)./((u(beta)).^1.*B(beta).^1)).^2;
equation2 = @(beta) (beta/k).^2.*(1./(u(beta).^2)+1./(w(beta).^2)) - (ncore./u(beta)).^2 - (nclad./w(beta)).^2;
equation3 = @(beta) (besseljDerivative(m, u(beta))./(u(beta).*besselj(m, u(beta))) + ( (besselkDerivative(m, w(beta)))./(w(beta).*besselk(m, w(beta))) ) ).*(besseljDerivative(m, u(beta))./(u(beta).*besselj(m, u(beta))) + (nclad/ncore)^2*( (besselkDerivative(m, w(beta)))./(w(beta).*besselk(m, w(beta))) ) )-m^2*(1./(u(beta).^2)+1./(w(beta).^2)).*(1./(u(beta).^2)+(nclad/ncore)^2*1./(w(beta).^2));
opt = optimoptions('fsolve', 'FunctionTolerance', 1e-20, 'OptimalityTolerance', 1e-20);
ourBeta = fsolve(equation3, [1.7e6], opt);
beta = ourBeta;

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

phase = pi/2;
A=1;
C=-A*beta/omega/m0*(m*(1/(u(beta)^2)+1/(w(beta^2)))/((besseljDerivative(m,u(beta))/(u(beta)*besselj(m,u(beta)))) + (besselkDerivative(m,w(beta))/(w(beta)*besselk(m,w(beta))))));
r_s = linspace(0.0001, a, 100)*1e-6;fi_s = linspace(0.0001, 2*pi, 100);
P = [];Ez = [];
for r=r_s
    
    for fi=fi_s
        
        
        
        if r>=0 && r<=a
            %Core region
            %E = [Er, Ephi, Ez];%H = [Hr, Hphi, Hz];
            E(1,1) = -1i*a^2/u(beta)^2*( A*beta*u(beta)/a*besseljDerivative(m,u(beta)*r/a) + C*omega*m0*m/r*besselj(m,u(beta)*r/a) )*cos(m*fi+phase);
            E(1,2) = -1i*a^2/u(beta)^2*( -A*beta*m/r*besselj(m,u(beta)*r/a) - C*omega*m0*u(beta)/a*besseljDerivative(m,u(beta)*r/a) )*sin(m*fi+phase);
            E(1,3) = A*besselj(m,u(beta)*r/a)*cos(m*fi+phase);
            H(1,1) = -1i*a^2/u(beta)^2*( A*omega*e0*ncore^2*m/r*besselj(m, u(beta)*r/a) + C*beta*u(beta)/a*besseljDerivative(m, u(beta)*r/a) )*sin(m*fi+phase);
            H(1,2) = -1i*a^2/u(beta)^2*( A*omega*e0*ncore^2*u(beta)/a*besseljDerivative(m, u(beta)*r/a) + C*beta*m/r*besselj(m, u(beta)*r/a) )*cos(m*fi+phase);
            H(1,3) = C*besselj(m,u(beta)*r/a)*sin(m*fi+phase);
        else
            %Cladding region
            %E = [Er, Ephi, Ez];%H = [Hr, Hphi, Hz];
            E(1,1) = 1i*a^2/w(beta)^2*( A*beta*w(beta)/a*besselkDerivative(m,w(beta)*r/a) + C*omega*m0*m/r*besselk(m,w(beta)*r/a) )*besselj(m, u(beta))/besselk(m,w(beta))*cos(m*fi+phase);
            E(1,2) = 1i*a^2/w(beta)^2*( -A*beta*m/r*besselk(m,w(beta)*r/a) - C*omega*m0*w(beta)/a*besselkDerivative(m,w(beta)*r/a) )*besselj(m, u(beta))/besselk(m,w(beta))*sin(m*fi+phase);
            E(1,3) = A*besselj(m,u(beta))/besselk(m, w(beta))*besselk(m,w(beta)*r/a)*cos(m*fi+phase);
            H(1,1) = 1i*a^2/w(beta)^2*( A*omega*e0*nclad^2*m/r*besselk(m, w(beta)*r/a) + C*beta*w(beta)/a*besselkDerivative(m, w(beta)*r/a) )*besselj(m, u(beta))/besselk(m,w(beta))*sin(m*fi+phase);
            H(1,2) = 1i*a^2/w(beta)^2*( A*omega*e0*nclad^2*w(beta)/a*besselkDerivative(m, w(beta)*r/a) + C*beta*m/r*besselj(m, w(beta)*r/a) )*besselj(m, u(beta))/besselk(m,w(beta))*cos(m*fi+phase);
            H(1,3) = C*besselj(m,u(beta))/besselk(m,w(beta))*besselk(m,w(beta)*r/a)*sin(m*fi+phase);
        end
        
        
        P(r_s==r, fi_s==fi) = 0.5*(E(1,1)*conj(H(1,2)) - E(1,2)*conj(H(1,1)));
        
    end
    
end

figure();
polarplot3d(abs(P), 'RadialRange',[0 a],'plottype','surfn', 'InterpMethod', 'spline','polargrid',{1 1})
