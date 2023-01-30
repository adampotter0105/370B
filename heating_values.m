function hhv = HHV_mass(gas)
    T0 = 298.15;
    P0 = oneatm;
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
    
    % Compare enthalpies
    gas0 = GRI30;
    set(gas0, "T", T0, "P", P0, "X", moleFractions(gas))
    n2 = GRI30; set(n2, "T", T0, "P", P0, "X", 'N2:1');
    o2 = GRI30; set(o2, "T", T0, "P", P0, "X", 'O2:1');
    h0 = enthalpy_mole(gas0) + n_in_n2*enthalpy_mole(n2) + n_in_o2*enthalpy_mole(o2); % [J/kmol]

    h20 = Water; set(h20, "T", T0, "P", P0);
    co2 = GRI30; set(co2, "T", T0, "P", P0, "X", 'CO2:1');
    h1 = n_out_n2*enthalpy_mole(n2) + n_out_h20*enthalpy_mole(h20) + n_out_co2*enthalpy_mole(co2) + n_out_o2*enthalpy_mole(o2); % [J/kmol]
    
    hhv = (h0-h1)/meanMolecularWeight(gas); % [J/kg]
end

function lhv = LHV_mass(gas)
    T0 = 298.15;
    P0 = oneatm;
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
    
    % Compare enthalpies
    gas0 = GRI30;
    set(gas0, "T", T0, "P", P0, "X", moleFractions(gas))
    n2 = GRI30; set(n2, "T", T0, "P", P0, "X", 'N2:1');
    o2 = GRI30; set(o2, "T", T0, "P", P0, "X", 'O2:1');
    h0 = enthalpy_mole(gas0) + n_in_n2*enthalpy_mole(n2) + n_in_o2*enthalpy_mole(o2); % [J/kmol]

    h20 = Water; set(h20, "T", T0, "P", P0);
    co2 = GRI30; set(co2, "T", T0, "P", P0, "X", 'CO2:1');
    h1 = n_out_n2*enthalpy_mole(n2) + n_out_h20*enthalpy_mole(h20) + n_out_co2*enthalpy_mole(co2) + n_out_o2*enthalpy_mole(o2); % [J/kmol]
    
    h2o_temp = Water; set(h2o_temp, "T", T0, "P", P0); 
    setState_satVapor(h2o_temp); h1_h2o = enthalpy_mole(h2o_temp); 
    setState_satLiquid(h2o_temp); h0_h2o = enthalpy_mole(h2o_temp);
    delta_h2o = h1_h2o - h0_h2o;

    lhv = (h0-h1 - delta_h2o*n_out_h20)/meanMolecularWeight(gas); % [J/kg]]
end