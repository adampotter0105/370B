function [total_w, H, X] = turb_work(eta, Pin, Pout, h0)
    w = Water;
    set(w, "H", h0, "P", Pin) % Initial Conditions of water
    n = 100; % discretization points for pressure
    h = (Pout-Pin)/n; % integration step
    total_w = 0; % total work output of turbine
    for p = Pin+h:h:Pout
        h1 = enthalpy_mass(w);
        s0 = entropy_mass(w);
        setState_SP(w, [s0, p]) % isentropic expansion
        h2 = enthalpy_mass(w);
        dw_s = h2-h1; % isentropic work
    
        dw = dw_s*eta; % actual work
        setState_HP(w, [h1+dw, p]) % actual new state of steam 
        total_w = total_w + dw; % sum all steps
    end
    X = vaporFraction(w);
    H = enthalpy_mass(w);
end

function [total_w, H] = pump_work(eta, Pin, Pout, h0)
    w = Water;
    set(w, "H", h0, "P", Pin) % Initial Conditions of water
    n = 100; % discretization points for pressure
    h = (Pout-Pin)/n; % integration step
    total_w = 0; % total work output of turbine
    for p = Pin+h:h:Pout
        h1 = enthalpy_mass(w);
        s0 = entropy_mass(w);
        setState_SP(w, [s0, p]) % isentropic expansion
        h2 = enthalpy_mass(w);
        dw_s = h2-h1; % isentropic work
    
        dw = dw_s/eta; % actual work
        setState_HP(w, [h1+dw, p]) % actual new state of steam 
        total_w = total_w + dw; % sum all steps
    end
    H = enthalpy_mass(w);
end