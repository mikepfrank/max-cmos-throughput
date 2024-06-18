%|=======================================================================|%
%|                  TOP OF FILE:  IRDS2022_HPHD.m                        |%
%|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~|%
%|                                                                       |%
%|      FILE NAME:  IRDS2022_HPHD.m                 [MATLAB script]      |%
%|                                                                       |%
%|      DESCRIPTION:                                                     |%
%|                                                                       |%
%|          This script carries out an analysis of the maximum           |%
%|          raw throughput densities (in switching events per unit       |%
%|          time per unit area) for conventional vs. adiabatic           |%
%|          switching in different CMOS technology generations           |%
%|          based on data from the 2022 edition of the International     |%
%|          Roadmap for Devices and Systems (IRDS).                      |%
%|                                                                       |%
%|      VERSION:    0.15w (beta)                                         |%
%|                                                                       |%
%|      AUTHORS:    Alexander J. Edwards (AJE)                           |%
%|                  Michael P. Frank (MPF)                               |%
%|                                                                       |%
%|vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv|%

close all;      % Close all figure windows.
clear;          % Clear variables in the workspace.

    % Select whether to analyze the high-performance (HP) or 
    % the high-density (HD) design point:

%is_HP = 1;  % Uncomment this line to select the HP design point.
is_HP = 0;  % Uncomment this line to select the HD design point.

    % [cols=Y] List of years of technology nodes being analyzed:

years = [2022, 2025, 2028, 2031, 2034, 2037];
NY = length(years);

    % [cols=P] Range of power density levels to be studied (W/cm^2):

% Manually selected rounded levels; 6 steps per decade; 13 levels.
%power_dens = [1, 1.5, 2.1, 3.2, 4.6, 6.8, 10, 15, 21, 32, 46, 68, 100];
% Computed levels:
power_dens = logspace(-2,2,49);   % 12 steps per decade; 49 levels.
NP = length(power_dens);

%|=========================================================================
%|  Hard-coded input data from tables in the IRDS More Moore chapter
%|  and associated baseline calculations.
%|vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [cols=Y]  Standard Vdd operating voltage level (V):

vdd = [0.7, 0.65, 0.65, 0.6, 0.6, 0.6];     % Baseline for both HP and HD.

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [cols=Y]  Ion per unit width (uA/um):

if(is_HP == 1)
    ion_um = [874, 787, 851, 753, 737, 753]; %HP
else
    ion_um = [644, 602, 656, 551, 532, 547]; %HD
end

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [cols=Y]  Ion per device (uA/dev):

% Commented out b/c we now calculate this ourselves from eff. device width.
%if(is_HP == 1)
%    ion_dev = [88, 170, 184, 157, 118, 115]; %HP
%else
%    ion_dev = [65, 130, 142, 115, 85, 83]; %HD
%end

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [cols=Y]  Ioff by um (pA/um):

if(is_HP == 1)
    ioff_um = [10000, 10000, 10000, 10000, 10000, 10000]; %HP
else
    ioff_um = [100, 100, 100, 100, 100, 100]; %HD
end

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [cols=Y]  Device effective width (nm):

if(is_HP == 1)
    W_dev = [101.0 216.0 216.0 208.0 160.0 152.0]; %HP
else
    W_dev = [101.0 126.0 96.0 128.0 88.0 80.0]; %HD
end

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [cols=Y]  Ioff by device (pA/dev):

ioff_dev = ioff_um .* W_dev * 0.001;

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [cols=Y]  Cell drive at saturation (ohms):

if(is_HP == 1)
    cell_drive = [1982, 1911, 1769, 1914, 2546, 2620]; %HP (HD not given)
else
    % Scale HP cell drive by difference in Ion/W and Weff to get HD.
    cell_drive_HP = [1982, 1911, 1769, 1914, 2546, 2620]; %HP (HD not given)
    ion_um_HP = [874, 787, 851, 753, 737, 753];         % HP Ion/W
    W_dev_HP = [101.0 216.0 216.0 208.0 160.0 152.0];   % HP Weff/dev
    cell_drive = cell_drive_HP .* (ion_um_HP ./ ion_um) .* (W_dev_HP ./ W_dev);
end

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [cols=Y]  Device to cell width multiplier (unitless):

d2c_mul = (vdd ./ cell_drive) ./ (ion_um * 1e-6 .* W_dev * 1e-3);
    % NOTE: This gives the ratio for the overall on-conductance of
    % the charging path through the cell (as a result of device 
    % sizing) relative to the on-conductance of a typical device.
    % We find it comes out as 2x for 2022, 4x for later years.

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [cols=Y]  Device Roff (G ohms):

Roff_dev = vdd ./ (ioff_dev * 1e-12) * 1e-9;

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [cols=Y]  Cell Roff (G ohms):

Roff_cell = Roff_dev ./ d2c_mul;

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [cols=Y]  Leakage power dissipation per cell at standard voltage
    %|      levels (W/cell):

Pleak_cell_std = vdd .* vdd ./ (Roff_cell * 1e9);

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [cols=Y]  NAND2-eq gate dens (Mgate/mm^2):

gate_dens = [21, 29, 39, 81, 171, 284]; 
    % We assume this is for HP, but it's not specified.
    % We don't know how to adjust this for HD, so we
    % assume it's the same (possibly poor assumption).

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [cols=Y]  Leakage power density at standard voltages (W/cm^2):

Pleak_dens_std = Pleak_cell_std .* gate_dens * 1e8;

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [cols=Y]  Wireload resistance (ohms):

if(is_HP == 1)
    wireload_resistance = [512, 625, 865, 972, 804, 824]; %HP
else
    wireload_resistance = [261, 298, 354, 446, 397, 405]; %HD
end

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [cols=Y]  Total path resistance (ohms):

total_path_resistance = cell_drive + wireload_resistance;

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [cols=Y]  Total load capacitance (fF):

if(is_HP == 1)
    load_cap = [5.9, 5.65, 5.35, 4.76, 4., 3.86]; %HP
else
    load_cap = [1.79, 1.69, 1.49, 1.36, 1.23, 1.20]; %HD
end

    %|~~~~~~~~~~~~~~~~~~~~~~~~
    %|  RC delay (FO3) (ps):

rc_delay = total_path_resistance .* load_cap / 1e3;

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  "Ideal" speed 1/(4RC) (GHz):

ideal_speed = 1./(4 * rc_delay * 1e-12) / 1e9;

% Commented this out because we now calculate this for ourselves.
%Energy per switching (fJ/switch)
%energy_per_sw = [0.65, 0.49, 0.47, 0.4, 0.4, 0.4]; %HP (HD not given)
% (This seems to be an intrinsic value not accounting for parasitics.)

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [cols=Y]  Energy per switching (fJ/switch):

energy_per_sw = load_cap .* vdd .* vdd;

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [cols=Y]  Threshold voltage Vth (mV):

if(is_HP == 1)
    vth = [156, 165, 165, 164, 156, 154]; %HP
else
    vth = [288, 271, 268, 268, 258, 255]; %HD
end

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [cols=Y]  Subthreshold Slope (mV/dec):

if(is_HP == 1)
    ss = [82, 72, 70, 70, 70, 70]; %HP
else
    ss = [75, 67, 67, 65, 65, 65]; %HD
end 

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [cols=Y]  Cell drive transconductance beta (mS/V):

beta = 1 ./ (cell_drive .* (vdd - vth * 0.001)) * 1000;

%|=========================================================================
%|  Calculations for conventional CMOS at standard voltage levels.
%|vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

max_leakage_frac_conv = 0.1;     % Leakage power at most 10% of total.
%max_leakage_frac_conv = 0.01;     % Leakage power at most 10% of total.

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y]  Gate density cap for leakage power = 10% 
    %|      (Mgates/mm^2):

dens_cap_conv = (power_dens' * max_leakage_frac_conv) ./ Pleak_cell_std / 1e8;

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y]  Maximum gate density consistent with capping
    %|      leakage to 10% of power density constraint (Mgates/mm^2):

gate_dens_aug_conven = ones(NP, 1) * gate_dens;
gate_dens_aug_conven(gate_dens_aug_conven > dens_cap_conv) = ...
	dens_cap_conv(gate_dens_aug_conven > dens_cap_conv);

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y]  Leakage Power Density (W/cm^2):

P_leak_conv = Pleak_cell_std .* (gate_dens_aug_conven * 1e8);

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y]  Power density available for switching (W/cm^2):

power_dens_subleak_conven = max(power_dens' - P_leak_conv, ...
	power_dens' * (1 - max_leakage_frac_conv));
    % Don't really need the max; first term is enough.

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y]  Maximum conventional frequency @ maximum gate 
    %|      density given power density constraint (MHz):

max_freq_conven = power_dens_subleak_conven ./ ...
	(gate_dens_aug_conven * 1e8 .* 1/2 .* energy_per_sw * 1e-15) / 1e6;

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y]  Slowdown factor for conventional switching
    %|    relative to the "ideal" max frequency (w/o power constraint):

y_conven = (ideal_speed(ones(1,NP),:) * 1e9) ./ (max_freq_conven * 1e6);

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y]  Maximum throughput density for conventional
    %|      switching at standard voltages (P bit-ops/s/cm^2):

conven_thruput = max_freq_conven * 1e6 .* gate_dens_aug_conven * 1e8 / 1e15;


%|=========================================================================
%|  Calculations for adiabatic CMOS at standard voltage levels.
%|vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

max_leakage_frac_adia = 0.5;     % Leakage power at most 50% of total.
%max_leakage_frac_adia = 0.1;     % Leakage power at most 10% of total.

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y]  Gate density cap for leakage power = 50% 
    %|      (Mgates/mm^2):

dens_cap_adia = power_dens' * max_leakage_frac_adia ./ Pleak_cell_std / 1e8;

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y]  Maximum gate density consistent with capping
    %|      leakage to 50% of power density constraint (Mgates/mm^2):

gate_dens_aug_adia = ones(NP, 1) * gate_dens;
gate_dens_aug_adia(gate_dens_aug_adia > dens_cap_adia) = ...
	dens_cap_adia(gate_dens_aug_adia > dens_cap_adia);

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y]  Leakage Power Density (W/cm^2):

P_leak_adia = Pleak_cell_std .* (gate_dens_aug_adia * 1e8);

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y]  Power density available for switching (W/cm^2):

power_dens_subleak_adia = max(power_dens' - P_leak_adia, power_dens' * ...
	(1 - max_leakage_frac_adia));
    % Don't really need the max; first term is enough.

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y]  Maximum frequency for (4-tick/cycle) adiabatic 
    %|      switching at max density within power constraint (MHz):

adia_max_freq = 0.5 * sqrt(power_dens_subleak_adia ./ ((gate_dens_aug_adia ...
	* 1e8) .* energy_per_sw * 1e-15 .* rc_delay * 1e-12)) / 1e6;

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y]  Slowdown factor for adiabatic clocking versus
    %|      a reference maximum frequency based on RC delay.

y = (ideal_speed(ones(1,NP),:) * 1e9) ./ (adia_max_freq * 1e6);

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y]  Speedup factor for adiabatic clocking versus
    %|      conventional.

x = adia_max_freq ./ max_freq_conven;

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y]  Maximum throughput density for adiabatic
    %|      switching at standard voltages (P ops/s/cm^2):

adia_thruput = adia_max_freq * 1e6 .* gate_dens_aug_adia * 1e8 / 1e15;


%|=========================================================================
%|  Calculations for conventional CMOS at optimized voltage levels.
%|vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y]  Optimal operating frequency given power density 
    %|      constrint for Conventional Voltage-scaled case. (MHz):

vsc_freq_final = zeros(NP,NY);

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y]  Optimal gate density for given power density
    %|      constraint constrained by leakage power consuming at 
    %|      most 10% of allotted power for the Conventional, Voltage-
    %|      scaled case. (Mgates/mm^2):

vsc_dens_final = zeros(NP,NY);

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y]  Maximum throughput density for the Conventional,
    %|      Voltage-scaled case. (Pops/cm^2):

sc_max_thruput = zeros(NP,NY);

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y] Optimal operating voltage (Vdd = v_sc) to achieve
    %|      maximum throughput density for given power density constraint,
    %|      for Conventional Voltage-scaled case. (V):
    
v_sc_final = zeros(NP,NY);

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  Non-linear optimization of throughput is achieved by sweeping v_sc
    %|      and identifying the maximum achieved throughput. 
    
for yr = 1:NY % 
    
        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [cols=V] Operating voltage test points. (V):

    %v_sc = 0.0001*vth(yr):0.001:1.5*vdd(yr);
    v_sc = 0.001*ss(yr):0.001:1.5*vdd(yr);
    NV=length(v_sc);

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [cols=V] Leakage power per cell. (W):
   
    Pleak_cell = v_sc .* v_sc / (Roff_cell(yr) * 1e9);

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [cols=V] Device on-conductance is calculated as a piecewise
        %|      function of v_sc.  G1_dev is calculated as a decaying 
        %|      exponential sub-threshold accoring to ss and ioff.  G2_dev 
        %|      is calculated as the approximate saturation conductance 
        %|      for v_sc = Vgs > Vth. (S):
    
    G1_dev = ioff_dev(yr) * 1e-12 ./ vdd(yr) * 10 .^ (v_sc/(ss(yr) * 0.001));
    G2_dev = beta(yr) * .001 .* (v_sc - vth(yr) * .001);

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [cols=V] Device on-conductance is calculated as a piecewise
        %|      function of v_sc.  Gtot = {G1 for v_sc = Vgs < Vsw; G2 for
        %|      v_sc = Vgs > Vsw}, where Vsw is the intersection of G1 and
        %|      G2 close to v_sc = Vth. (S):
    
    Gtot_dev = G2_dev;
    Gtot_dev((G2_dev < G1_dev) & (v_sc < (vth(yr) * 0.0011))) = ...
		G1_dev((G2_dev < G1_dev) & (v_sc < (vth(yr) * 0.0011)));

    % Save the Gtot data for later plotting
    if(yr == 1)
        Gtot_2022 = Gtot_dev;
    end
    if (yr == 2)
        Gtot_2025 = Gtot_dev;
    end
    if (yr == 3)
        Gtot_2028 = Gtot_dev;
    end
    if (yr == 4)
        Gtot_2031 = Gtot_dev;
    end
    if (yr == 5)
        Gtot_2034 = Gtot_dev;
    end
    if (yr == 6)
        Gtot_2037 = Gtot_dev;
    end
    
        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [cols=V] Cell on-conductance as a function of v_sc. (S):
    
    Gtot_cell = Gtot_dev * d2c_mul(yr);
    
        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [rows=P,cols=V] Gate density assuming leakage power
        %|      consumes exactly 10% of allotted power. (Mgates/mm^2):
        
    max_dens = max_leakage_frac_conv * power_dens' ./ Pleak_cell / 1e8;

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [rows=P,cols=V] Maximum gate density assuming leakage power
        %|      consumes at most 10% of allotted power. (Mgates/mm^2):
    
    gate_dens_aug_vsc = ones(NP,NV) * gate_dens(yr);
    gate_dens_aug_vsc(gate_dens_aug_vsc > max_dens) = ...
		max_dens(gate_dens_aug_vsc > max_dens);

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [rows=P,cols=V] Recalculated leakage power density with 
        %|      updated gate density. (W/cm^2):
    
    P_leak_Vsc = v_sc .* v_sc ./ (Roff_cell(yr) .* 1e9) .* ...
		(gate_dens_aug_vsc * 1e8);
    
        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [cols=V] RC delay according to calculated cell drive and load 
        %|      resistance. (ps):
    
    RC_delay_sc = (1 ./ Gtot_cell + wireload_resistance(yr)) * ...
		load_cap(yr) * 1e-15 * 1e12;

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [cols=V] Maximum operating frequency assuming a safe clock
        %|      period is at most 4RC. (MHz):
    
    max_scaled_freq = 1 / 4 ./ (RC_delay_sc * 1e-12) / 1e6;

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [cols=V] Energy per switching event in voltage scaled regime.
        %|      (fJ/switch):
    
    energy_per_sw_sc = energy_per_sw(yr) .* v_sc .* v_sc / vdd(yr) / vdd(yr);

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [rows=P,cols=V] Maximum frequency is capped by power density 
        %|      constraint and RC delay after voltage scaling. (MHz):
    
    max_freq_sc = min((power_dens' - P_leak_Vsc) ./ ...
		((gate_dens_aug_vsc * 1e8) .* (energy_per_sw_sc / 2) * 1e-15) ...
		/ 1e6, max_scaled_freq(ones(1,NP),:));

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [rows=P,cols=V] Maxmimum throughput for power density contraint
        %|      and scaled voltage. (P bit-ops/s/cm^2):
    
    sc_thruput = (max_freq_sc * 1e6) .* (gate_dens_aug_vsc * 1e8) / 1e15;

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [rows=P,cols=Y] Maxmimum throughput for power density 
        %|      contraint with any v_sc. (P bit-ops/s/cm^2):
    
    [sc_max_thruput(:,yr), final_index] = max(sc_thruput');
    
        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [rows=P,cols=Y] Optimal v_sc granting the maxmimum throughput.
        %|      (V):
    
    v_sc_final(:,yr) = v_sc(final_index);

    for p=1:NP
        
            %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %|  [rows=P,cols=Y]  Optimal operating frequency given power 
            %|      density constrint. (MHz):
        
        vsc_freq_final(p,yr) = max_freq_sc(p,final_index(p));
        
            %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            %|  [rows=P,cols=Y]  Optimal gate density for given power 
            %|      density constraint. (Mgates/mm^2):
        
        vsc_dens_final(p,yr) = gate_dens_aug_vsc(p,final_index(p));
    end

end

%|=========================================================================
%|  Calculations for adiabatic CMOS at optimized voltage levels.
%|vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  Minimum factor of clock period increase compared with RC time 
    %|      constant for proper adiabatic charging. Based on heuristics: 
    %|      (T_adia >> 2*RC). (unitless):

adiabatic_slowdown = 4; 

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y] Maximum throughput density given power density 
    %|      constraint for Adiabatic, Voltage-scaled case.
    %|      (P bit-ops/s/cm^2):

adia_sc_max_thruput = zeros(NP, NY);

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y] Optimal scaled voltage achieving maximum throughput
    %|      density given power density constraint for Adiabatic, Voltage-
    %|      scaled case. (V):

adia_v_sc = zeros(NP, NY);

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y] Optimal operating frequency achieving maximum 
    %|      throughput density given power density constraint for 
    %|      Adiabatic, Voltage-scaled case. (MHz):

adia_freq = zeros(NP, NY);

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  [rows=P,cols=Y] Optimal gate density for given power density
    %|      constraint constrained by leakage power consuming at 
    %|      most 50% of allotted power for the Adiabatic, Voltage-
    %|      scaled case. (Mgates/mm^2)

adia_dens = zeros(NP, NY);

    %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    %|  Non-linear optimization of throughput is achieved by sweeping v_sc
    %|      and identifying the maximum achieved throughput. 

for yr = 1:NY
    for v_sc = 0.001*ss(yr):0.0001:1.5*vdd(yr)

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  CV^2 Energy per switching event in voltage scaled regime.
        %|      (fJ/switch):
    
        energy_per_sw_sc = energy_per_sw(yr) .* v_sc .* v_sc ...
			./ vdd(yr) ./ vdd(yr);

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  Leakage power per cell. (W):
        
        Pleak_cell = v_sc * v_sc / (Roff_cell(yr) * 1e9);

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [rows=P] Gate density assuming leakage power consumes exactly 
        %|      50% of allotted power. (Mgates/mm^2):
        
        max_dens = max_leakage_frac_adia * power_dens' / Pleak_cell / 1e8;

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [rows=P] Maximum gate density assuming leakage power consumes 
        %|      at most 50% of allotted power. (Mgates/mm^2):
        
        gate_dens_aug_adia_vsc = ones(NP, 1) * gate_dens(yr);
        gate_dens_aug_adia_vsc(gate_dens_aug_adia_vsc > max_dens) = ...
			max_dens(gate_dens_aug_adia_vsc > max_dens);

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [rows=P] Leakage power calculated from density per cell 
        %|      with updated gate density. (W/cm^2):
        
        P_leak_Vsc = Pleak_cell * (gate_dens_aug_adia_vsc * 1e8);

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [rows=P] Switching power density available subtracting leakage
        %|      power density. (W/cm^2):
        
        Psw_dens = power_dens' - P_leak_Vsc;

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  Device on conductance is calculated as a piecewise
        %|      function of v_sc.  G1_dev is calculated as a decaying 
        %|      exponential sub-threshold accoring to ss and ioff.  G2_dev 
        %|      is calculated as the approximate saturation conductance 
        %|      for v_sc = Vgs > Vth. (S):
        
        G1_dev = ioff_dev(yr) * 1e-12 ./ vdd(yr) * 10 .^ (v_sc/(ss(yr) * 0.001));
        G2_dev = beta(yr) * .001 .* (v_sc - vth(yr) * .001);

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  Device on conductance is calculated as a piecewise
        %|      function of v_sc.  Gtot = {G1 for v_sc = Vgs < Vsw; G2 for
        %|      v_sc = Vgs > Vsw}, where Vsw is the intersection of G1 and
        %|      G2 close to v_sc = Vth. (S):
        
        Gtot_dev = G2_dev;
        Gtot_dev((G2_dev < G1_dev) & (v_sc < (vth(yr) * 0.0011))) = ...
			G1_dev((G2_dev < G1_dev) & (v_sc < (vth(yr) * 0.0011)));

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  Cell on conductance as a function of v_sc. (S):
        
        Gtot_cell = Gtot_dev * d2c_mul(yr);

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  RC delay according to calculated cell drive and load 
        %|      resistance. (ps):
        
        RC_delay_sc = (1 ./ Gtot_cell + wireload_resistance(yr)) * ...
			load_cap(yr) * 1e-15 * 1e12;
        
        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  Maximum operating frequency assuming a minimum clock period 
        %|      of 4RC is necessary for propoer adiabatic charging. (MHz):
        
        max_freq_given_RC = 1e-6/(RC_delay_sc * adiabatic_slowdown * 1e-12);

        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [rows=P] Adiabatic optimal frequency given power available for 
        %|      switching & cell density. (MHz):
        
        adia_opt_freq = 0.5 * sqrt(Psw_dens ./ ((gate_dens_aug_adia_vsc * ...
			1e8) * (energy_per_sw_sc * 1e-15) * (rc_delay * 1e-12))) / 1e6;
        
        %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [rows=P] Cap the adiabatic frequency based on the scaled RC
        %|      delay.

        for power_index = 1:NP
            if (adia_opt_freq(power_index) > max_freq_given_RC)
                adia_opt_freq(power_index) = max_freq_given_RC;
            end
        end

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [rows=P] Maxmimum throughput for power density contraint
        %|      and scaled voltage. (P bit-ops/s/cm^2):
        
        adia_sc_thruput = (adia_opt_freq * 1e6) .* (gate_dens_aug_adia_vsc ...
			* 1e8) / 1e15;

        %|~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        %|  [rows=P,cols=Y] Update best solution to identify maximum 
        %|      throughput density for each power level and year.
        
        for power_index = 1:NP
            if (adia_sc_thruput(power_index) > ...
				adia_sc_max_thruput(power_index,yr))
				
                adia_sc_max_thruput(power_index,yr) = ...
					adia_sc_thruput(power_index);
                adia_v_sc(power_index, yr) = v_sc;
                adia_freq(power_index, yr) = adia_opt_freq(power_index);
                adia_dens(power_index, yr) = gate_dens_aug_adia_vsc(power_index);
            end
        end
    end
end

%|=========================================================================
%|  Figure plotting.
%|vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Figure 1: Quadrillions of Switching Operations per Second per Square 
%     Centimeter (Pops/cm^2) under various power density constraints 
%     (W/cm^2) for Conventional, Adiabatic, Voltage-scaled, and Adiabatic 
%     voltage-scaled methodologies. 

figure;

ms = 5;  % marker size

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Conventional (non-adiabatic, non-voltage-scaled):  SQUARE MARKERS.

loglog(power_dens, conven_thruput(:,1), 'rs-','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','red');
hold on;
grid on;
loglog(power_dens, conven_thruput(:,2), 'ys-','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','yellow');
loglog(power_dens, conven_thruput(:,3), 'gs-','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','green');
loglog(power_dens, conven_thruput(:,4), 'cs-','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','cyan');
loglog(power_dens, conven_thruput(:,5), 'bs-','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','blue');
loglog(power_dens, conven_thruput(:,6), 'ms-','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','magenta');

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Adiabatic (non-voltage-scaled):                    CIRCLE MARKERS.

loglog(power_dens, adia_thruput(:,1),'ro-','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','red');
loglog(power_dens, adia_thruput(:,2),'yo-','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','yellow');
loglog(power_dens, adia_thruput(:,3),'go-','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','green');
loglog(power_dens, adia_thruput(:,4),'co-','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','cyan');
loglog(power_dens, adia_thruput(:,5),'bo-','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','blue');
loglog(power_dens, adia_thruput(:,6),'mo-','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','magenta');

xlabel('Power Density Constraint, Watts/cm^2');
ylabel('Quadrillions of Switching Operations per Second per Square Centimeter (Pops/cm^2)');

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Voltage-scaled:                                     DIAMOND MARKERS.

loglog(power_dens, sc_max_thruput(:,1), 'rd--','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','red');
loglog(power_dens, sc_max_thruput(:,2), 'yd--','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','yellow');
loglog(power_dens, sc_max_thruput(:,3), 'gd--','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','green');
loglog(power_dens, sc_max_thruput(:,4), 'cd--','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','cyan');
loglog(power_dens, sc_max_thruput(:,5), 'bd--','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','blue');
loglog(power_dens, sc_max_thruput(:,6), 'md--','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','magenta');

    %~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    % Adiabatic voltage-scaled:                           STAR MARKERS.

loglog(power_dens, adia_sc_max_thruput(:,1),'rp--','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','red');
loglog(power_dens, adia_sc_max_thruput(:,2),'yp--','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','yellow');
loglog(power_dens, adia_sc_max_thruput(:,3),'gp--','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','green');
loglog(power_dens, adia_sc_max_thruput(:,4),'cp--','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','cyan');
loglog(power_dens, adia_sc_max_thruput(:,5),'bp--','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','blue');
loglog(power_dens, adia_sc_max_thruput(:,6),'mp--','LineWidth',1.5, ...
	'MarkerSize',ms,'MarkerFaceColor','magenta');

title('Throughput Density vs. Power Density');

set(gcf,'position',[10,10,400,700]);

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Figure 2: Optimal gate density (Mgates/mm^2) for various power density 
%     constraints (W/cm^2) constrained by leakage power consuming at most 
%     10% of allotted power for the Conventional, Standard-voltage case.

figure;
loglog(power_dens, gate_dens_aug_conven(:,1), 'r');
hold on;
grid on;
loglog(power_dens, gate_dens_aug_conven(:,2), 'y');
loglog(power_dens, gate_dens_aug_conven(:,3), 'g');
loglog(power_dens, gate_dens_aug_conven(:,4), 'c');
loglog(power_dens, gate_dens_aug_conven(:,5), 'b');
loglog(power_dens, gate_dens_aug_conven(:,6), 'm');
xlabel('Power Density Constraint, Watts/cm^2');
ylabel('Optimal Gate Density, Standard-Voltage Conventional (Mgates/mm^2)')
title('Optimal Gate Density vs. Power Density', ...
    'Standard-Voltage Conventional Switching');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Figure 3: Optimal operating frequency (MHz) for various power density 
%     constraints (W/cm^2) for the Conventional, Standard-voltage case.

figure;
%semilogx(power_dens, vsc_freq_final(:,1),'r');
loglog(power_dens, max_freq_conven(:,1),'r');
hold on;
grid on;
plot(power_dens, max_freq_conven(:,2),'y');
plot(power_dens, max_freq_conven(:,3),'g');
plot(power_dens, max_freq_conven(:,4),'c');
plot(power_dens, max_freq_conven(:,5),'b');
plot(power_dens, max_freq_conven(:,6),'m');
xlabel('Power Density Constraint, Watts/cm^2');
ylabel('Optimal Frequency, Standard-Voltage Conventional (MHz)')
title('Optimal Switching Frequency vs. Power Density', ...
    'Standard-Voltage Conventional Switching');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Figure 4: Optimal gate density (Mgates/mm^2) for various power density 
%     constraints (W/cm^2) constrained by leakage power consuming at most 
%     50% of allotted power for the Adiabatic, Standard-voltage case.

figure;
loglog(power_dens, gate_dens_aug_adia(:,1), 'r');
hold on;
grid on;
loglog(power_dens, gate_dens_aug_adia(:,2), 'y');
loglog(power_dens, gate_dens_aug_adia(:,3), 'g');
loglog(power_dens, gate_dens_aug_adia(:,4), 'c');
loglog(power_dens, gate_dens_aug_adia(:,5), 'b');
loglog(power_dens, gate_dens_aug_adia(:,6), 'm');
xlabel('Power Density Constraint, Watts/cm^2');
ylabel('Optimal Gate Density, Standard-Voltage Adiabatic (Mgates/mm^2)')
title('Optimal Gate Density vs. Power Density', ...
    'Standard-Voltage Adiabatic Switching');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Figure 5: Optimal operating frequency (MHz) for various power density
%     constraints (W/cm^2) for the Adiabatic, Standard-voltage case.

figure;
%semilogx(power_dens, vsc_freq_final(:,1),'r');
loglog(power_dens, adia_max_freq(:,1),'r');
hold on;
grid on;
plot(power_dens, adia_max_freq(:,2),'y');
plot(power_dens, adia_max_freq(:,3),'g');
plot(power_dens, adia_max_freq(:,4),'c');
plot(power_dens, adia_max_freq(:,5),'b');
plot(power_dens, adia_max_freq(:,6),'m');
xlabel('Power Density Constraint, Watts/cm^2');
ylabel('Optimal Frequency, Standard-Voltage Adiabatic (MHz)')
title('Optimal Switching Frequency vs. Power Density', ...
    'Standard-Voltage Adiabatic Switching');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Figure 6: Optimal scaled operating voltage (V) for various power density 
%     constraints (W/cm^2) to maximize throughput for the Conventional,
%     Voltage-scaled case.

figure;
semilogx(power_dens, v_sc_final(:,1),'r');
hold on;
grid on;
yline(vdd(1),'r--');
yline(vth(1)/1000, 'r--')
plot(power_dens, v_sc_final(:,2),'y');
yline(vdd(2),'y--');
yline(vth(2)/1000, 'y--')
plot(power_dens, v_sc_final(:,3),'g');
yline(vdd(3),'g--');
yline(vth(3)/1000, 'g--')
plot(power_dens, v_sc_final(:,4),'c');
yline(vdd(4),'c--');
yline(vth(4)/1000, 'c--')
plot(power_dens, v_sc_final(:,5),'b');
yline(vdd(5),'b--');
yline(vth(5)/1000, 'b--')
plot(power_dens, v_sc_final(:,6),'m');
yline(vdd(6),'m--');
yline(vth(6)/1000, 'm--')
ylim([0 vdd(2)*1.1]);
xlabel('Power Density Constraint, Watts/cm^2');
ylabel('Optimal Scaled Voltage, Conventional (V)')
title('Optimal Operating Voltage vs. Power Density', ...
    'for Conventional Switching');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Figure 7: Optimal gate density (Mgates/mm^2) for various power density 
%     constraints (W/cm^2) constrained by leakage power consuming at most 
%     10% of allotted power for the Conventional, Voltage-scaled case.

figure;
loglog(power_dens, vsc_dens_final(:,1), 'r');
hold on;
grid on;
loglog(power_dens, vsc_dens_final(:,2), 'y');
loglog(power_dens, vsc_dens_final(:,3), 'g');
loglog(power_dens, vsc_dens_final(:,4), 'c');
loglog(power_dens, vsc_dens_final(:,5), 'b');
loglog(power_dens, vsc_dens_final(:,6), 'm');
xlabel('Power Density Constraint, Watts/cm^2');
ylabel('Optimal Gate Density, Voltage-Scaled Conventional (Mgates/mm^2)')
title('Optimal Gate Density vs. Power Density', ...
    'Low-Voltage Conventional Switching');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Figure 8: Optimal operating frequency (MHz) for various power density
%     constraints (W/cm^2) for the Conventional, Voltage-scaled case.

figure;
%semilogx(power_dens, vsc_freq_final(:,1),'r');
loglog(power_dens, vsc_freq_final(:,1),'r');
hold on;
grid on;
plot(power_dens, vsc_freq_final(:,2),'y');
plot(power_dens, vsc_freq_final(:,3),'g');
plot(power_dens, vsc_freq_final(:,4),'c');
plot(power_dens, vsc_freq_final(:,5),'b');
plot(power_dens, vsc_freq_final(:,6),'m');
xlabel('Power Density Constraint, Watts/cm^2');
ylabel('Optimal Frequency, Voltage-Scaled Conventional (MHz)')
title('Optimal Switching Frequency vs. Power Density', ...
    'Low-Voltage Conventional Switching');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Figure 9: Optimal scaled operating voltage (V) for various power density 
%     constraints (W/cm^2) to maximize throughput for the Adiabatic,
%     Voltage-scaled case.

figure;
semilogx(power_dens, adia_v_sc(:,1),'r');
hold on;
grid on;
yline(vdd(1),'r--');
yline(vth(1)/1000, 'r--')
plot(power_dens, adia_v_sc(:,2),'y');
yline(vdd(2),'y--');
yline(vth(2)/1000, 'y--')
plot(power_dens, adia_v_sc(:,3),'g');
yline(vdd(3),'g--');
yline(vth(3)/1000, 'g--')
plot(power_dens, adia_v_sc(:,4),'c');
yline(vdd(4),'c--');
yline(vth(4)/1000, 'c--')
plot(power_dens, adia_v_sc(:,5),'b');
yline(vdd(5),'b--');
yline(vth(5)/1000, 'b--')
plot(power_dens, adia_v_sc(:,6),'m');
yline(vdd(6),'m--');
yline(vth(6)/1000, 'm--')
ylim([0 vdd(2)*1.1]);
xlabel('Power Density Constraint, Watts/cm^2');
ylabel('Optimal Scaled Voltage, Adiabatic (V)')
title('Optimal Operating Voltage vs. Power Density', ...
    'for Adiabatic Switching');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Figure 10: Optimal gate density (Mgates/mm^2) for various power density 
%     constraints (W/cm^2) constrained by leakage power consuming at most 
%     50% of allotted power for the Adiabatic, Voltage-scaled case.

figure;
%semilogx(power_dens, adia_freq(:,1),'r');
loglog(power_dens, adia_dens(:,1),'r');
hold on;
grid on;
plot(power_dens, adia_dens(:,2),'y');
plot(power_dens, adia_dens(:,3),'g');
plot(power_dens, adia_dens(:,4),'c');
plot(power_dens, adia_dens(:,5),'b');
plot(power_dens, adia_dens(:,6),'m');
xlabel('Power Density Constraint, Watts/cm^2');
ylabel('Optimal Gate Density, Voltage-Scaled Adiabatic (Mgates/mm^2)')
title('Optimal Gate Density vs. Power Density', ...
    'Low-Voltage Adiabatic Switching');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Figure 11: Optimal operating frequency (MHz) for various power density
%     constraints (W/cm^2) for the Adiabatic, Voltage-scaled case.

figure;
%semilogx(power_dens, adia_freq(:,1),'r');
loglog(power_dens, adia_freq(:,1),'r');
hold on;
grid on;
plot(power_dens, adia_freq(:,2),'y');
plot(power_dens, adia_freq(:,3),'g');
plot(power_dens, adia_freq(:,4),'c');
plot(power_dens, adia_freq(:,5),'b');
plot(power_dens, adia_freq(:,6),'m');
xlabel('Power Density Constraint, Watts/cm^2');
ylabel('Optimal Frequency, Voltage-Scaled Adiabatic (MHz)')
title('Optimal Switching Frequency vs. Power Density', ...
    'Low-Voltage Adiabatic Switching');

%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Figure 12: Device on-conductance (S) for various scaled voltages, Vdd =
%     vsc (V).

figure;
semilogy(Gtot_2022, 'r');
hold on;
grid on;
semilogy(Gtot_2025, 'y');
semilogy(Gtot_2028, 'g');
semilogy(Gtot_2031, 'c');
semilogy(Gtot_2034, 'b');
semilogy(Gtot_2037, 'm');
xlabel('Index in vsc array');
ylabel('On-conductance for Vdd=vsc (S)');
title('Device On-Conductance','Across Range of Considered Voltages');
