function [] = PlotMode(E_comp, H_comp, r_s, fi_s, a, b, plot_type)


if plot_type == "whole"
    figure();
    polarplot3d(real(E_comp(:,:,3)), 'RadialRange',[0 b],'plottype','surfn', 'InterpMethod', 'spline','polargrid',{1 1});
    view([0 90])
    colormap jet
elseif plot_type == "core"
    figure();
    polarplot3d(real(E_comp(r_s<=a,:,3)), 'RadialRange',[0 a],'plottype','surfn', 'InterpMethod', 'spline','polargrid',{1 1});
    view([0 90])
else
    
    error('');
    
end

end

