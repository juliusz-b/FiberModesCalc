function [] = PlotMode(plothook,E_comp, H_comp, r_s, fi_s, a, b, plot_type)

%gca = plothook;

if plot_type == "whole"
    polarplot3d(real(E_comp(:,:,3)), 'RadialRange',[0 b],'plottype','surfn', 'InterpMethod', 'spline','polargrid',{1 1});
    hold on;
    max_amp=max(real(E_comp(:,:,3)));
    theta = linspace(0,2*pi,100);
    plot3(a*cos(theta),a*sin(theta),ones(1,length(theta)),'-',LineWidth=2,Color='w');
elseif plot_type == "core"
    polarplot3d(real(E_comp(r_s<=a,:,3)), 'RadialRange',[0 a],'plottype','surfn', 'InterpMethod', 'spline','polargrid',{1 1});
else
    error('Wrong plot type');
end
%view(plothook,[0 90])
xlabel(plothook,'Radius [\mu m]')
axis(plothook,"square");
zlim(plothook,[-1,1]);
colormap(plothook,"jet")
c = colorbar(plothook);
c.Label.String = 'Amplitude';
c.Location="EastOutside";
set(plothook, 'fontsize', 15);
hold off
end

