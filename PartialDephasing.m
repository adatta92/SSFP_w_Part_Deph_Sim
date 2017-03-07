% This calls the constant flow simulation code or displays saved data.

clear all, close all, clc

set(0, 'defaultlinelinewidth', 1);

load_sim = true; % false; %   

d_phi = 2.5;

flow_vect = (0:28)/28;

if ~load_sim % Run the simulations

	tic

	spectral_profs_0 = zeros(length(flow_vect), 360/d_phi);
	spectral_profs_30 = zeros(length(flow_vect), 360/d_phi);
	spectral_profs_60 = zeros(length(flow_vect), 360/d_phi);

	parfor ind = 1:length(flow_vect)
	    spectral_profs_0(ind, :) = part_dep_sim_parallel(flow_vect(ind), 0, d_phi);
	    spectral_profs_30(ind, :) = part_dep_sim_parallel(flow_vect(ind), 30, d_phi);
	    spectral_profs_60(ind, :) = part_dep_sim_parallel(flow_vect(ind), 60, d_phi);
	end %parfor

	toc

	on_res_static = part_dep_sim_parallel(0, 0, 360);
	spectral_profs = cat(3, spectral_profs_0, spectral_profs_30, spectral_profs_60)/on_res_static;

	save('PD_before_RF.mat', 'spectral_profs')
	% save('PD_after_RF.mat', 'spectral_profs') % See part_dep_sim_parallel.m for how to change it to the other partial dephasing implementation.

else % Display the data from previously-run simulations 

	% close all

	load('PD_before_RF.mat')
	% load('PD_after_RF.mat')

	color_param = linspace(0, 1, length(flow_vect))';
	newDefaultColors = [color_param 1-color_param 1-color_param];

	spectral_profs(:, 1, :)

	for i = 1:3
		figure, set(gca, 'ColorOrder', newDefaultColors, 'NextPlot', 'replacechildren');
		    plot(d_phi:d_phi:360, (spectral_profs(:, :, i))')
		    xlabel('Precession per TR (^o)'), ylabel('|M_{xy}| (a.u.)')
		    axis([0 360 0 18])
	end

	spin_replacement_rate = 100 * flow_vect * 3.5/6; % 3.5 ms TR, 6 mm FWHM slice 

	figure, imagesc(round(spin_replacement_rate)), colormap(newDefaultColors), cbr_h = colorbar;
	set(cbr_h, 'YTick', round(spin_replacement_rate(1:2:end)))

end
