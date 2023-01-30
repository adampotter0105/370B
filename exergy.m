function X = exergy_mass(gas)
    [X_tm, X_c] = exergy_calc(gas);
    X = X_tm + X_c;
end

function [X_tm, X_c] = exergy_calc(gas)
    % Dead State
    T0 = 25+273.15;
    P0 = oneatm;
    R = 8.314; % universal gas constant
    dead = GRI30;  
    set(dead, "T", T0, "P", P0, 'X', 'N2:0.757223,O2:0.202157,H2O:0.031208,Ar:0.009015,CO2:0.000397')
    % Thermo-Mechanical Equillibrium
    P_in = pressure(gas);
    T_in = temperature(gas);
    H_in = enthalpy_mass(gas);
    tm_eq = gas; 
    set(tm_eq, "T", T0, "P", P0);
    gibbs_tm = gibbs_mass(tm_eq);
    U_tm = intEnergy_mass(tm_eq);
    rho_tm = 1/density(tm_eq);
    S_tm = entropy_mass(tm_eq);
    try
        set(gas, "T", T_in, "P", P_in); % Return gas to normal
    catch
        setState_TH(gas,[T_in, H_in])
    end
    
    % Thermo-Mechanical Exergy
    X_tm = (intEnergy_mass(gas)-U_tm) + P0*(1/density(gas) - rho_tm) - T0*(entropy_mass(gas)-S_tm); % J/kg
    
    

    if size(moleFractions(gas),1)>1 % Object is likely GRI30
        mu = chemPotentials(dead);
    
        % Decomposing Input to Environmental 
        % [ C, H, N, O, Ar]
        mm_all = molecularWeights(gas);
        mm = [mm_all(speciesIndex(gas,'C')), mm_all(speciesIndex(gas,'H')), mm_all(speciesIndex(gas,'N')), mm_all(speciesIndex(gas,'O')), mm_all(speciesIndex(gas,'Ar'))]; % molar masses (g/mol)
        mass_frac = [ elementalMassFraction(gas,'C'),  elementalMassFraction(gas,'H'),  elementalMassFraction(gas,'N'), elementalMassFraction(gas,'O'),  elementalMassFraction(gas,'Ar')];
        mole_frac = (mass_frac./mm)./(sum(mass_frac./mm)); % normalize to 1mol of elements
        avg_atoms = meanMolecularWeight(gas)/sum(mm.*mole_frac); % average atoms per molecule
        comb_frac = mole_frac*avg_atoms; % atoms of each element from combusting a mole of input gas
    
        % Assign Inputs and Outputs
        n_out_co2 = comb_frac(1);
        n_out_h20 = comb_frac(2)/2;
        n_in_o2 = max(n_out_h20/2 + n_out_co2 - (comb_frac(4)/2), 0);
        n_out_o2 = max((comb_frac(4)/2) - n_out_h20/2 - n_out_co2 , 0);
        n_in_n2 = 3.76*n_in_o2;
        n_out_n2 = n_in_n2 + comb_frac(3)/2;
        n_in_ar = n_in_o2*(0.009015/0.202157);
        n_out_ar = n_in_ar + comb_frac(5);
        
        % Sum Chemical Potentials
        mu_in = n_in_o2*(mu(speciesIndex(dead, 'O2')) + R*T0*log(moleFraction(dead,'O2'))) ... 
            + n_in_n2*(mu(speciesIndex(gas, 'N2')) + R*T0*log(moleFraction(dead,'N2'))) ...
            + n_in_ar*(mu(speciesIndex(gas, 'Ar')) + + R*T0*log(moleFraction(dead,'Ar')));
        mu_out = n_out_co2*(mu(speciesIndex(dead, 'CO2')) + R*T0*log(moleFraction(dead,'CO2'))) ...
            + n_out_h20*(mu(speciesIndex(gas, 'H2O')) + R*T0*log(moleFraction(dead,'H2O'))) ...
            + n_out_n2*(mu(speciesIndex(gas,'N2')) + R*T0*log(moleFraction(dead,'N2'))) ...
            + n_out_ar*(mu(speciesIndex(gas, 'Ar')) + R*T0*log(moleFraction(dead,'Ar'))) ...
            + n_out_o2*(mu(speciesIndex(gas, 'O2')) + R*T0*log(moleFraction(dead,'O2')));
        % These chem potentials are [J/kmol of gas in]
        X_c = gibbs_tm - (mu_out - mu_in)/meanMolecularWeight(gas); % J/Kg
    else
        X_c = 0;
    end
    
end

function X = flowExergy_mass(gas)
    P_in = pressure(gas);
    P0 = oneatm;
    X_in = exergy_mass(gas);
    X_ex = ((1/density(gas))*(P_in - P0)); % [J/kg]
    % Ignoring potential and kinetic energy
    X = X_in + X_ex;
end
