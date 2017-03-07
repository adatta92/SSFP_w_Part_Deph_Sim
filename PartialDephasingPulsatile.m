% This generates the pulsatile flow waveforms and calls the simulation function.

d_phi = 10;

tic

r = [cosspace(0, 10, 75) 10*ones(1, 25) -cosspace(-10, 0, 150) zeros(1, 50)]/14;
r = [r r];
r = r(1:end-100);

M_xy_0 = part_dep_sim_pulsatile(r, 0, d_phi);

toc
tic

% M_xy_30 = part_dep_sim_pulsatile(r, 30, d_phi);
M_xy_60 = part_dep_sim_pulsatile(r, 60, d_phi);

save('pulsatile_fast.mat')

r = 0.75*r;

M_xy_0 = part_dep_sim_pulsatile(r, 0, d_phi);
% M_xy_30 = part_dep_sim_pulsatile(r, 30, d_phi);
M_xy_60 = part_dep_sim_pulsatile(r, 60, d_phi);

save('pulsatile_medium.mat')

toc
