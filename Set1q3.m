% Set 1 Problem 3

clc;

% givens
gas = HFC134a();
x6 = 1;
T6 = 5+273.15;
x3 = 0;
T3 = 50+273.15;
Ecomp = 0.85;
Eex = 0.8;
% mflow = ;

% calculate pressures for all states and entropies for known states
set(gas,'T',T6,'Vapor',x6);
P6 = pressure(gas);
h6 = enthalpy_mass(gas);
P1 = 0.9*P6;
P5 = P6/.9;

set(gas,'T',T3,'Vapor',x3);
P3 = pressure(gas);
h3 = enthalpy_mass(gas);
P2 = P3/0.9;
P4 = 0.9*P3;

% Calculate heat transferred through exchanger
set(gas,'T',T3,'P',P1);
ha = enthalpy_mass(gas);
set(gas,'T',T6,'Vapor',x6); %use instead of P6 because it is saturated state
hb = enthalpy_mass(gas);
qmax = ha-hb; % J/K
qex = Eex*qmax; % J/K

% calculate h1 and h4 from qex
h1 = qex+h6;
h4 = h3-qex;
% set h5=h4 through throttling valve
h5 = h4;

% calculate state 2 using compressor work
[win, h2] = comp_work(Ecomp, P1, P2, h1);

% calculate COP and capacity
qevap = h6-h5;
COP = qevap/win;
capacity = qevap/1e3; % kJ/kg

pressures = [P1 P2 P3 P4 P5 P6 P1];
enthalpies = [h1 h2 h3 h4 h5 h6 h1];
semilogy(enthalpies/1e3, pressures/1e6, '--ok', 'LineWidth',1.7);
hold on;

% Create Vapor Dome
% Set up a fluid object to work with in Cantera.
fluid = HFC134a();

% Get the critical point props.
Tc = critTemperature(fluid);
Pc = critPressure(fluid);
set(fluid,'T',Tc,'P',Pc);
sc = entropy_mass(fluid);
hc = enthalpy_mass(fluid);
vc = 1/density(fluid);

% Set the triple point props.  Use a value epsilon above Tt to get the
% bottom edge of the vapor dome.
Tt = 170;
Pt = satPressure(fluid,Tt);
setState_Tsat(fluid,[Tt 0]);
sft = entropy_mass(fluid);
hft = enthalpy_mass(fluid);
vft = 1/density(fluid);
setState_Tsat(fluid,[Tt 1]);
sgt = entropy_mass(fluid);
hgt = enthalpy_mass(fluid);
vgt = 1/density(fluid);

% Set the limits for data curves.
Tmin = Tt+1;
Tmax = 450;
Pmin = 1.01*Pt;
Pmax = 500e5;

% Make a vapor dome.
dT = (Tc-Tt)/100;
i = 1;
%setState_Tsat(fluid,[Tmin 0]);
for T=Tt:dT:Tc-dT    % Stop short of the critical point.
    Tsatline(i) = T;
    setState_Tsat(fluid,[T 0]);
    sliqline(i) = entropy_mass(fluid);
    hliqline(i) = enthalpy_mass(fluid);
    Pliqline(i) = pressure(fluid);
    setState_Tsat(fluid,[T 1]);
    svapline(i) = entropy_mass(fluid);
    hvapline(i) = enthalpy_mass(fluid);
    Pvapline(i) = pressure(fluid);
    i = i+1;
end
Tsatline(i) = Tc;   % Add the critical point now.
sliqline(i) = sc;
hliqline(i) = hc;
Pliqline(i) = Pc;
svapline(i) = sc;
hvapline(i) = hc;
Pvapline(i) = Pc;

% Start a set of isotherms.
Tlist = [Tmin 200 250 300 350 400 450];
for j=1:1:3    % Do the part of list below the critical point. 
    fluid = Solution('liquidvapor.xml','hfc134a');
    T = Tlist(j);

    % Do the compressed liquid side.
    setState_Tsat(fluid,[T 0]);
    Psat = pressure(fluid);
    logdP = (log(Pmax) - log(Psat))/50;
    i = 1;
    for logP=log(Pmax):-logdP:log(Psat)+logdP  % Stop short of saturation.
        P = exp(logP);
        Ttempline(j,i) = T;
        set(fluid,'T',T,'P',P);
        stempline(j,i) = entropy_mass(fluid);
        htempline(j,i) = enthalpy_mass(fluid);
        Ptempline(j,i) = pressure(fluid);
        i = i+1;
    end
    Ttempline(j,i) = T;   % Add the saturation point now.
    setState_Tsat(fluid,[T 0]);
    stempline(j,i) = entropy_mass(fluid);
    htempline(j,i) = enthalpy_mass(fluid);
    Ptempline(j,i) = pressure(fluid);
    i = i+1;

    % Now go across the dome.
    dq = 1/50;
    for q=0+dq:dq:1-dq  % Stop short of saturation.
        Ttempline(j,i) = T;
        setState_Psat(fluid,[Psat q]);
        stempline(j,i) = entropy_mass(fluid);
        htempline(j,i) = enthalpy_mass(fluid);
        Ptempline(j,i) = pressure(fluid);
        i = i+1;
    end
    Ttempline(j,i) = T; % Add the saturation point now.
    setState_Psat(fluid,[Psat 1]);
    stempline(j,i) = entropy_mass(fluid);
    htempline(j,i) = enthalpy_mass(fluid);
    Ptempline(j,i) = pressure(fluid);
    i = i+1;

    % Do the vapor side.
    logdP = (log(Psat) - log(Pmin))/50;
    for logP=log(Psat)-logdP:-logdP:log(Pmin)  % Stop short of saturation.
        P = exp(logP);
        Ttempline(j,i) = T;
        set(fluid,'T',T,'P',P);
        stempline(j,i) = entropy_mass(fluid);
        htempline(j,i) = enthalpy_mass(fluid);
        Ptempline(j,i) = pressure(fluid);
        i = i+1;
    end
end

% Add isotherms above the critical temperature.
for j=j+1:1:length(Tlist)
    fluid = Solution('liquidvapor.xml','hfc134a');
    T = Tlist(j);

    logdP = (log(Pmax) - log(Pmin))/150;
    i = 1;
    for logP=log(Pmax):-logdP:log(Pmin)  % Stop short of saturation.
        P = exp(logP);
        Ttempline(j,i) = T;
        set(fluid,'T',T,'P',P);
        stempline(j,i) = entropy_mass(fluid);
        htempline(j,i) = enthalpy_mass(fluid);
        Ptempline(j,i) = pressure(fluid);
        i = i+1;
    end
end

% Start a set of isentropes.
slist     = [sft   500   1000   1500];
Phighlist = [Pmax  Pmax  Pmax   10e6];
lPlow = log(Pmin);
lPhigh = log(10e6);
for j=1:1:length(slist)
    s = slist(j);
    lPhigh = log(Phighlist(j));
    ldP = (lPhigh - lPlow)/150;
    i = 1;
    for lP=lPlow:ldP:lPhigh  % Stop short of saturation.
        P = exp(lP);
        Pentrline(j,i) = P;
        setState_SP(fluid,[s,P]);
        hentrline(j,i) = enthalpy_mass(fluid);
        i = i+1;
    end
end


figure(1)
i = 1;
semilogy(htempline(i,:)/1e3,Ptempline(i,:)/1e6,'r', 'LineWidth',1.7);
semilogy(hentrline(i,:)/1e3,Pentrline(i,:)/1e6,'Color','b', 'LineWidth',1.7)
semilogy(hc/1e3,Pc/1e6,'kd', 'LineWidth',1.7);
%semilogy([hft hgt]/1e3,[Pt Pt]/1e6,'ko--');
for i=1:1:length(Tlist)
    semilogy(htempline(i,:)/1e3,Ptempline(i,:)/1e6,'r', 'LineWidth',1.7)
end
semilogy(hliqline/1e3,Pliqline/1e6,'k', 'LineWidth',1.7)
semilogy(hvapline/1e3,Pvapline/1e6,'k', 'LineWidth',1.7)
%semilogy([hft hgt]/1e3,[Pt Pt]/1e6,'ko--')

xlabel('Specific Enthalpy (kJ/kg)')
ylabel('Pressure (MPa)')
for i=1:1:length(Tlist)
    text(htempline(i,150)/1e3,Ptempline(i,150)/1e6,...
        num2str(Tlist(i)),'Color','r')
end

i=1;
semilogy(hentrline(i,:)/1e3,Pentrline(i,:)/1e6,'Color','b')
for i=1:1:length(slist)
    semilogy(hentrline(i,:)/1e3,Pentrline(i,:)/1e6,'Color','b', 'LineWidth',1.7)
end
for i=1:1:length(slist)
    text(hentrline(i,150)/1000,Pentrline(i,150)/1e6,...
        num2str(slist(i)/1e3),'Color','b')
end

% plot cycle
semilogy(enthalpies/1e3, pressures/1e6, '--ok', 'LineWidth',1.7);
hold on;
for i = 1:6
    text(enthalpies(i)/1e3+7,pressures(i)/1e6-.1,int2str(i), "FontWeight", "bold");
end
ylim([5*1e-5, 1e2]);
text(-190, .0002, strcat("COP: ", num2str(COP)), "FontWeight", "bold");
text(-190, .0001, strcat("Capacity: ", num2str(capacity), " kJ/kg"), "FontWeight", "bold");
legend('Vapor Compression Cycle','Isotherms (K)', 'Isentropes (kJ/kg-K)');
title('HFC134a Vapor-Compression Refrigeration');
hold off;
%plotfixer;

% gives the total work and the enthalpy of the compressor
function [total_w, H] = comp_work(eta, Pin, Pout, h0)
    w = HFC134a;
    set(w, "H", h0, "P", Pin) % Initial Conditions of fluid
    n = 101; % discretization points for pressure
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

