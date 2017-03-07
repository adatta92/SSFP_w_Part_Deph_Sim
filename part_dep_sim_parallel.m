function spectral_prof = part_dep_sim_parallel(r, PD_amt, d_phi, TBW)
    % r is the flow rate in m/s (e.g. 0.8 m/s = 80 cm/s).
    % PD_amt is the amount of partial dephasing in degrees (e.g., 0, 30, or 60).
    % d_phi is the spectral resolution (e.g. 2.5 degrees).
    % TBW is the time-bandwidth of the RF pulse whose profile is stored in flip_prof_TBW# (e.g. flip_prof_TBW2 for TBW = 2).
   
    % spectral_prof is the net transverse signal at each precession frequency, including in-slice and out-of-slice spins. 
 
    TR = 3.5; 
    T1 = 1000; T2 = 150; 
    T2star = 75; 
    
    PC_amt = 180;

    if nargin < 3
        d_phi = 2.5;
    end
    if nargin < 4
        TBW = 2;
    end

    load(sprintf('flip_prof_TBW%d', TBW))
    del_z = 6;
    z_res = 0.125;
    
    r_res = z_res/TR;
    r = round(r/r_res) * r_res;
    shift_amt = round(r/r_res);

    % assuming TE = TR/2;
    max_phi = 360; % 720 for a full spectral period (180 degree phase difference between the two half-periods)
    
    [A2_inslice, b2] = relaxation(TR/2, T1, T2); % E_rel, b_rec
    [A2_outslice, ~] = relaxation(TR/2, T1, T2star);

    if abs(r) > eps
        z = floor(-2.5*del_z/2):z_res:min(ceil(2.5*del_z/2) + ceil(4*T2star * r), 150); % Assuming 150 mm coil sensitivity extent.
        numTRs = ceil(length(z)/shift_amt);
    else
        z = floor(-2.5*del_z/2):z_res:ceil(2.5*del_z/2);
        numTRs = 250;
    end
    
    phi_PD = z/del_z * PD_amt;
    
    alpha = interp1(pos, flip_prof, z, 'linear', 0);
    
    M = zeros(3, length(z), max_phi/d_phi);
    M(3, :) = 1;

    for TRnum = 1:numTRs     
        M = cat(2, repmat([0; 0; 1], 1, shift_amt, max_phi/d_phi), M(:, 1:end-shift_amt, :));
        
        phi_PC = mod(PC_amt * TRnum, 360); % FM_amt * TRnum^2/2 +
        
%         (length(z)-shift_amt*(numTRs-TRnum)) 
        
        for z_ind = 1:(length(z)-shift_amt*(numTRs-TRnum)) 
                    
               A3 = throt(alpha(z_ind), phi_PC); % R_exc
               A4 = zrot(phi_PD(z_ind)); % R_pd % Unbalanced gradient before excitation (after readout)
%              A3 = zrot(phi_PD(z_ind)); % Unbalanced gradient after excitation (before readout)
%              A4 = throt(alpha(z_ind), phi_PC);

%            par -- takes longer than regular for
            for phi_ind = 1:(max_phi/d_phi) % phi_offres = d_phi:d_phi:max_phi         
                A1 = zrot(d_phi * phi_ind/2); % R_off

                if abs(alpha) > eps % in-slice
                    A = A2_inslice*A1*A3*A4*A2_inslice*A1; % A2_inslice*A1*A4*A3*A4*A2_inslice*A1; % 
                    b = A2_inslice*A1*A3*A4*b2 + b2; % A2_inslice*A1*A4*A3*A4*b2_inslice + b2_inslice; % 
                else
                    A = A2_outslice*A1*A3*A4*A2_outslice*A1; % A2_outslice*A1*A4*A3*A4*A2_outslice*A1; % 
                    b = A2_outslice*A1*A3*A4*b2 + b2; % A2_outslice*A1*A4*A3*A4*b2_outslice + b2_outslice; % 
                end %if

                M(:, z_ind, phi_ind) = A*M(:, z_ind, phi_ind) + b;
            end %for
        end %for
    end %for
    
    M_xy = squeeze((M(1, :, :) + j*M(2, :, :)).*exp(-j*pi/180*phi_PC));
    
    spectral_prof = abs(sum(M_xy))/(del_z/z_res);
    
    save_evol_flag = false; % true; % Set to true to save m_xy(z), the evolution of the transverse magnetization as it flows through/the spatial distribution of the magnetization, instead of just the total signal.  This can't be run in parallel.
                
    if save_evol_flag             
        save(sprintf('M_xy_%.0f_SR_%d_PD.mat', 10000*r*TR/del_z, PD_amt), 'z', 'M_xy')
    end %if   
end %function
