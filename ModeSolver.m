function [E_comp, H_comp, r_s, fi_s] = ModeSolver(lambda, a, b, A, ML, phase, n_s, type)
%Juliusz Bojarczuk 2020. Based on: https://doi.org/10.1016/B978-0-12-525096-2.X5000-4
%lambda - wavelength
%a - core RADIUS
%b - cladding RADIUS
%A - amplitude component related with power carried my mode. More info: 10.1016/B978-0-12-525096-2.X5000-4
%ML - when LP is calculated ML(1) = m and ML(2) = l. In other cases: ML(1)= m and ML(2) = 1;
%phase - initial phase for calculations
%n_s(1) - ncore; n_s(2) - ncladd
%type: 'LP' or 'hybrid'

%%% Ph. cons.
e0 = 8.854187e-12;
m0 = 4*pi*1e-7;
c=physconst('lightspeed');
%%%
ncore = n_s(1);nclad = n_s(2);
m = ML(1);
l = ML(2);
%%%
k = 2*pi/(lambda);
beta_range = [k*nclad k*ncore];
omega = 2*pi/(lambda)* c;
V = a*k*(ncore^2-nclad^2)^0.5;
%%% Bessel's derivatives definition
if m==0
    besselkDerivative = @(m, x) -besselk(1,x);
    besseljDerivative = @(m, x) -besselj(1,x);
else
    besselkDerivative = @(m, x) (-0.5*(besselk(m-1,x)+besselk(m+1,x)));
    besseljDerivative = @(m, x) (0.5*(besselj(m-1,x)-besselj(m+1,x)));
end
%%% u and w definitions
u = @(beta) a*sqrt(k^2*ncore^2-beta.^2);
w = @(beta) a*sqrt(beta.^2 - k^2*nclad^2);
%%% dispersion equation
if type == 'LP'
    if m==0
        equation = @(beta) ((besselj(0,u(beta)))./(u(beta).*(besselj(1,u(beta)))) - (besselk(0,w(beta)))./(w(beta).*(besselk(1,w(beta)))));
    elseif m==1
        equation = @(beta) ((besselj(1,u(beta)))./(u(beta).*(besselj(0,u(beta)))) + (besselk(1,w(beta)))./(w(beta).*(besselk(0,w(beta)))));
    else
        equation = @(beta) ((besselj(m,u(beta)))./(u(beta).*(besselj(m-1,u(beta)))) + (besselk(m,w(beta)))./(w(beta).*(besselk(m-1,w(beta)))));
    end
elseif type == 'hybrid'
    equation = @(beta) (besseljDerivative(m, u(beta))./(u(beta).*besselj(m, u(beta))) + ( (besselkDerivative(m, w(beta)))./(w(beta).*besselk(m, w(beta))) ) ).*(besseljDerivative(m, u(beta))./(u(beta).*besselj(m, u(beta))) + (nclad/ncore)^2*( (besselkDerivative(m, w(beta)))./(w(beta).*besselk(m, w(beta))) ) )-m^2*(1./(u(beta).^2)+1./(w(beta).^2)).*(1./(u(beta).^2)+(nclad/ncore)^2*1./(w(beta).^2));
else
    error('Wrong mode type!');
end
%%% mode search
mode_beta=[];
for bet=beta_range(1)+0.01:0.001:beta_range(2)-0.01
    mode_beta((beta_range(1)+0.01:0.001:beta_range(2)-0.01)==bet) = fsolve(equation,[bet]);
end
mode_beta(mode_beta<beta_range(1))=[];
mode_beta(mode_beta>beta_range(2))=[];

if isempty(mode_beta)
   error('Cannot found appropriate beta for given parameters'); 
end

mode_beta(diff(mode_beta)<0.00001)=[];
mode_beta = flip(sort(mode_beta));
%%% plotting dispersion equation vs beta
figure();
betas = [beta_range(1):.001:beta_range(2)];
plot(betas, equation(betas), 'linewidth', 3);
hold on;
plot(mode_beta, equation(mode_beta),'o', 'linewidth', 3);ylim([-20 20]);xlim(beta_range);
xlabel('\beta');ylabel('Dispersion equation value');
set(gca, 'fontsize', 15);
%%% choosing beta
beta = mode_beta(l);
%%% calculating C as it is shown in 10.1016/B978-0-12-525096-2.X5000-4
C=-A*beta/omega/m0*(m*(1/(u(beta)^2)+1/(w(beta^2)))/((besseljDerivative(m,u(beta))/(u(beta)*besselj(m,u(beta)))) + (besselkDerivative(m,w(beta))/(w(beta)*besselk(m,w(beta))))));
r_s = linspace(0.0001, b, 250);fi_s = linspace(0.0001, 2*pi, 250);
P = [];Ez = [];E=[];H=[];E_comp=[];H_comp=[];
for r=r_s
    
    for fi=fi_s
        
        s = (m*(1/u(beta)^2+1/w(beta)^2))/((besseljDerivative(m, u(beta))/(u(beta)*besselj(m, u(beta))))+(besselkDerivative(m,w(beta))/(w(beta)*besselk(m,u(beta)))) );
        Er = -1i*A*beta*a/u(beta)*((1-s)/2*besselj(m-1,u(beta)*r/a)-(1+s)/2*besselj(m+1, u(beta)*r/2))*cos(m*fi+phase);
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
        
        E_comp(r_s==r, fi_s==fi, :) = E;
        H_comp(r_s==r, fi_s==fi, :) = H;
        
    end
    
end

end

