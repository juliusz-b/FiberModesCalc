function [] = PlotMode(E_comp, H_comp, r_s, fi_s, b)



figure();
polarplot3d(real(E_comp(:,:,3)), 'RadialRange',[0 b],'plottype','surfn', 'InterpMethod', 'spline','polargrid',{1 1});
view([0 90])
colormap jet

end

