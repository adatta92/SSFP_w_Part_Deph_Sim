function spectral_prof_evol = part_dep_sim_pulsatile(r, PD_amt, d_phi, TBW, numOutTRs)
    % r is a vector of the flow rate during each TR in m/s (e.g. 0.8 m/s = 80 cm/s).
    % PD_amt is the amount of partial dephasing in degrees (e.g., 0, 30, or 60).
    % d_phi is the spectral resolution (e.g. 2.5 degrees).
    % TBW is the time-bandwidth of the RF pulse whose profile is stored in flip_prof_TBW# (e.g. flip_prof_TBW2 for TBW = 2).
    % numOutTRs is the number of TRs during which the spectral profile is outputted.  The outputted TRs correspond to the last <numOutTRs> in the flow rate vector.  This is because the initial TRs are affected by SSFP transients. 
   
    % spectral_prof_evol is the net transverse signal during each of the outputted TRs at each precession frequency in a numOutTRs by ~360/d_phi matrix.

    TR = 3.5; 
    T1 = 1000; T2 = 150; 
    T2star = 75; 
    
    PC_amt = 180; 

    if nargin < 4
        TBW = 2;
    end

    if nargin < 5
        numOutTRs = 200;
    end
    
    load(sprintf('flip_prof_TBW%d', TBW))
    del_z = 6;
    z_res = 0.125; %  0.25; % 0.5; % 
    
    r_res = z_res/TR;
    r = round(r/r_res) * r_res;
    shift_amt = round(r/r_res);

    % assuming TE = TR/2;

    max_phi = 360; % 720;
    
    [A2_inslice, b2_inslice] = relaxation(TR/2, T1, T2);
    [A2_outslice, b2_outslice] = relaxation(TR/2, T1, T2star);

    z = floor(-2.5*del_z/2):z_res:min(ceil(2.5*del_z/2) + ceil(4*T2star * max(r)), 150); % Using the peak velocity ensures that a long enough distance is simulated to allow the out-of-slice signal to decay.
    numTRs = length(r);
    
    phi_PD = z/del_z * PD_amt;
    
    alpha = interp1(pos, flip_prof, z, 'linear', 0);
    
    M = zeros(3, length(z), max_phi/d_phi);
    M(3, :) = 1;

    on_res_static = part_dep_sim_parallel(0, 0, 360);
    
    spectral_prof_evol = zeros(numOutTRs, length(d_phi:d_phi:max_phi));
%     figure
    for TRnum = 1:numTRs     
        M = cat(2, repmat([0; 0; 1], 1, shift_amt(TRnum), max_phi/d_phi), M(:, 1:end-shift_amt(TRnum), :));
        
        phi_PC = mod(PC_amt * TRnum, 360); 
        
        for z_ind = 1:(length(z)-sum(shift_amt(TRnum+1:(numTRs-numOutTRs)))) % The rest will exit the region before the first plotted profile.
                    
            A3 = throt(alpha(z_ind), phi_PC);
            A4 = zrot(phi_PD(z_ind));

            for phi_offres = d_phi:d_phi:max_phi         
                A1 = zrot(phi_offres/2);

                if alpha % in-slice
                    A = A2_inslice*A1*A3*A4*A2_inslice*A1; 
                    b = A2_inslice*A1*A3*A4*b2_inslice + b2_inslice; 
                else
                    A = A2_outslice*A1*A3*A4*A2_outslice*A1;
                    b = A2_outslice*A1*A3*A4*b2_outslice + b2_outslice;
                end %if

                M(:, z_ind, phi_offres/d_phi) = A*M(:, z_ind, phi_offres/d_phi) + b;
            end %for
        end %for

        if numTRs - TRnum < numOutTRs          
            M_xy = squeeze((M(1, :, :) + j*M(2, :, :)).*exp(-j*pi/180*phi_PC));
            spectral_prof_evol(numOutTRs - (numTRs - TRnum), :) = sum(M_xy)/(del_z/z_res)/on_res_static; % Normalization
        end
    end %for  
end %function
