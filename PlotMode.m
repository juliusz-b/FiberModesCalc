function [] = PlotMode(plothook,E_comp, H_comp, r_s, fi_s, a, b, plot_type)

%gca = plothook;

if plot_type == "whole"
    %figure();
    polarplot3d(real(E_comp(:,:,3)), 'RadialRange',[0 b],'plottype','surfn', 'InterpMethod', 'spline','polargrid',{1 1});
    view([0 90])
    xlabel('Radius [\mu m]')
    colormap jet
    c = colorbar;
    ylabel(c, 'Amplitude')
elseif plot_type == "core"
    %figure();
    polarplot3d(real(E_comp(r_s<=a,:,3)), 'RadialRange',[0 a],'plottype','surfn', 'InterpMethod', 'spline','polargrid',{1 1});
    view(plothook,[0 90])
    xlabel(plothook,'Radius [\mu m]')
    colormap(plothook,"jet")
    c = colorbar(plothook);
    c.Label.String = 'Amplitude';
    c.Location="EastOutside";
    set(plothook, 'fontsize', 15);
else
    
    error('Wrong plot type');
    
end

end

