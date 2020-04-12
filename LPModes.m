clear all;clc;close all;
a = 20.3/2;
neff = 1.447;
nclad = 1.4440;
ncore = 1.4513;DELTA=(ncore-nclad)/ncore*100;%1.4513
lambda = 1;
e0 = 8.854187e-12;
m0 = 4*pi*1e-7;
LP_type = [3,1];
m=LP_type(1);l=LP_type(2);

k = 2*pi/(lambda);
beta_range = [k*nclad k*ncore];
c=physconst('lightspeed');
omega = 2*pi/(lambda)* c;
V = a*k*(ncore^2-nclad^2)^0.5;
if m==0
    besselkDerivative = @(m, x) -besselk(1,x);
    besseljDerivative = @(m, x) -besselj(1,x);
else
    besselkDerivative = @(m, x) (-0.5*(besselk(m-1,x)+besselk(m+1,x)));
    besseljDerivative = @(m, x) (0.5*(besselj(m-1,x)-besselj(m+1,x)));
end
u = @(beta) a*sqrt(k^2*ncore^2-beta.^2);
w = @(beta) a*sqrt(beta.^2 - k^2*nclad^2);
B = @(beta) (w(beta).^2)./(u(beta).^2+w(beta).^2);

if m==0
   equation = @(beta) ((besselj(0,u(beta)))./(u(beta).*(besselj(1,u(beta)))) - (besselk(0,w(beta)))./(w(beta).*(besselk(1,w(beta)))));
elseif m==1
   equation = @(beta) ((besselj(1,u(beta)))./(u(beta).*(besselj(0,u(beta)))) + (besselk(1,w(beta)))./(w(beta).*(besselk(0,w(beta)))));
else
   equation = @(beta) ((besselj(m,u(beta)))./(u(beta).*(besselj(m-1,u(beta)))) + (besselk(m,w(beta)))./(w(beta).*(besselk(m-1,w(beta)))));
end

betas = [beta_range(1):.001:beta_range(2)];
aaaa=[];bbbb=[];
opt = optimoptions('fsolve', 'FunctionTolerance', 1e-20, 'OptimalityTolerance', 1e-20);
for bet=beta_range(1)+0.01:0.001:beta_range(2)-0.01
    aaaa((beta_range(1)+0.01:0.001:beta_range(2)-0.01)==bet) = fsolve(equation,[bet]);
end
aaaa(aaaa<beta_range(1))=[];
aaaa(aaaa>beta_range(2))=[];
aaaa(diff(aaaa)<0.00001)=[];
aaaa = flip(sort(aaaa));

figure();
plot(betas, equation(betas), 'linewidth', 3);
hold on;
plot(aaaa, equation(aaaa),'o', 'linewidth', 3);ylim([-20 20]);xlim(beta_range);
xlabel('\beta');ylabel('Dispersion equation value');

beta = aaaa(l);


phase = 0;
A=1;
C=-A*beta/omega/m0*(m*(1/(u(beta)^2)+1/(w(beta^2)))/((besseljDerivative(m,u(beta))/(u(beta)*besselj(m,u(beta)))) + (besselkDerivative(m,w(beta))/(w(beta)*besselk(m,w(beta))))));
r_s = linspace(0.0001, 125/2, 100);fi_s = linspace(0.0001, 2*pi, 100);
P = [];Ez = [];
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
        
        
        P(r_s==r, fi_s==fi) = E(1,3);
        %P(r_s==r, fi_s==fi) = 0.5*(E(1,1)*conj(H(1,2)) - E(1,2)*conj(H(1,1)));
        Ex(r_s==r, fi_s==fi) = E(1,1)*cos(fi)-E(1,2)*sin(fi);
        Ey(r_s==r, fi_s==fi) = E(1,1)*sin(fi)+E(1,2)*cos(fi);
        Hx(r_s==r, fi_s==fi) = H(1,1)*cos(fi)-H(1,2)*sin(fi);
        Hy(r_s==r, fi_s==fi) = H(1,1)*sin(fi)+H(1,2)*cos(fi);
        Pxy(r_s==r, fi_s==fi) = 0.5*(Ex(r_s==r, fi_s==fi)*conj(Hy(r_s==r, fi_s==fi)) - Ey(r_s==r, fi_s==fi)*conj(Hx(r_s==r, fi_s==fi)));
        Pxy(r_s==r, fi_s==fi) = Ex(r_s==r, fi_s==fi)+Ey(r_s==r, fi_s==fi);
  
        
    end
    
end
figure();
polarplot3d(real(P), 'RadialRange',[0 a],'plottype','surfn', 'InterpMethod', 'spline','polargrid',{1 1});
view([0 90])
colormap jet

figure();
subplot(2,2,1)
polarplot3d(abs(Ex), 'RadialRange',[0 a],'plottype','surfn', 'InterpMethod', 'spline','polargrid',{1 1});
view([0 90])
title('Ex')
colormap jet

subplot(2,2,2)
polarplot3d(abs(Hx), 'RadialRange',[0 a],'plottype','surfn', 'InterpMethod', 'spline','polargrid',{1 1});
view([0 90])
title('Hx')
colormap jet

subplot(2,2,3)
polarplot3d(abs(Ey), 'RadialRange',[0 a],'plottype','surfn', 'InterpMethod', 'spline','polargrid',{1 1});
view([0 90])
title('Ey')
colormap jet

subplot(2,2,4)
polarplot3d(abs(Hy), 'RadialRange',[0 a],'plottype','surfn', 'InterpMethod', 'spline','polargrid',{1 1});
view([0 90])
title('Hy')
colormap jet

