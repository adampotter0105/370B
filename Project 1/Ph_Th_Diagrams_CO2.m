% Make T-s and P-h diagrams for carbon dioxide.
% C.F. Edwards, 12/21/10

clear all
format compact
fprintf('\n********************************************************\n')

% Set up a fluid object to work with in Cantera.
fluid = Solution('liquidvapor.cti','carbondioxide');

% Get the critical point props.
Tc = critTemperature(fluid)
Pc = critPressure(fluid)
set(fluid,'T',Tc,'P',Pc);
sc = entropy_mass(fluid);
hc = enthalpy_mass(fluid);

% Set the triple point props.  Use a value epsilon above Tt to get the
% bottom edge of the vapor dome.
Tt = 216.6
Pt = satPressure(fluid,Tt)
setState_Tsat(fluid,[Tt 0]);
sft = entropy_mass(fluid);
hft = enthalpy_mass(fluid);
vft = 1/density(fluid);
setState_Tsat(fluid,[Tt 1]);
sgt = entropy_mass(fluid);
hgt = enthalpy_mass(fluid);
vgt = 1/density(fluid);

% Set the limits for data curves.
Tmin = Tt;
Tmax = 500;
Pmin = Pt;
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

% Start a set of isobaric curves.
Plist = [Pt/1e5 10 20 50];  % Pressure in bar
for j=1:1:length(Plist)
    fluid = Solution('liquidvapor.xml','carbondioxide');
    P = Plist(j)*1e5

    % Do the compressed liquid side.
    setState_Psat(fluid,[P 0]);
    Tsat = temperature(fluid)
    dT = (Tsat - Tmin)/50;
    for i=1:1:50  % Stop short of saturation.
        T = Tmin + (i-1)*dT;
        Tpresline(j,i) = T;
        set(fluid,'T',T,'P',P);
        hpresline(j,i) = enthalpy_mass(fluid);
    end
    i = i+1;
    setState_Psat(fluid,[P 0]);
    Tpresline(j,i) = Tsat;   % Add the saturation point now.
    hpresline(j,i) = enthalpy_mass(fluid);
    
    % Now go across the dome.
    i = i+1;
    dq = 1/50;
    for q=0+dq:dq:1-dq  % Stop short of saturation.
        Tpresline(j,i) = Tsat;
        setState_Psat(fluid,[P q]);
        hpresline(j,i) = enthalpy_mass(fluid);
        i = i+1;
    end
    Tpresline(j,i) = Tsat;   % Add the saturation point now.
    setState_Psat(fluid,[P 1]);
    hpresline(j,i) = enthalpy_mass(fluid);
    
    % Do the vapor side.
    i = i+1;
    dT = (Tmax - Tsat)/50;
    for ii=1:1:50
        T = Tsat + ii*dT;  % Start just above saturation.
        Tpresline(j,i) = T;
        set(fluid,'T',T,'P',P);
        hpresline(j,i) = enthalpy_mass(fluid);
        i = i+1;
    end
end

% Add isobars above the critical pressure.
Plist = [Plist 100 200 500];  % Add to the list
dT = (Tmax - Tmin)/150;
for j=j+1:1:length(Plist);
    P = Plist(j)*1e5;
    i = 1;
    for T=Tmin:dT:Tmax  % Stop short of saturation.
        Tpresline(j,i) = T;
        set(fluid,'T',T,'P',P);
        hpresline(j,i) = enthalpy_mass(fluid);
        i = i+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start a set of isotherms.
Tlist = [Tmin 250 275 300 350 400 500];
for j=1:1:4    % Do the part of list below the critical point. 
    fluid = Solution('liquidvapor.xml','carbondioxide');
    T = Tlist(j)

    % Do the compressed liquid side.
    setState_Tsat(fluid,[T 0]);
    Psat = pressure(fluid);
    logdP = (log(Pmax) - log(Psat))/50;
    for i=1:1:50    % Stop short of saturation.
        logP = log(Pmax) - (i-1)*logdP
        P = exp(logP);
        Ttempline(j,i) = T;
        set(fluid,'T',T,'P',P);
        htempline(j,i) = enthalpy_mass(fluid);
        Ptempline(j,i) = pressure(fluid);
    end
    i = i+1;
    Ttempline(j,i) = T;   % Add the saturation point now.
    setState_Tsat(fluid,[T 0]);
    htempline(j,i) = enthalpy_mass(fluid);
    Ptempline(j,i) = pressure(fluid);
    
    % Now go across the dome.
    i = i+1;
    dq = 1/50;
    for q=0+dq:dq:1-dq  % Stop short of saturation.
        Ttempline(j,i) = T;
        setState_Psat(fluid,[Psat q]);
        htempline(j,i) = enthalpy_mass(fluid);
        Ptempline(j,i) = pressure(fluid);
        i = i+1;
    end
    Ttempline(j,i) = T; % Add the saturation point now.
    setState_Psat(fluid,[Psat 1]);
    htempline(j,i) = enthalpy_mass(fluid);
    Ptempline(j,i) = pressure(fluid);
    
    % Do the vapor side.
    i = i+1;
    logdP = (log(Psat) - log(Pmin))/50;
    for ii=1:1:50
        logP = log(Psat) - (ii)*logdP
        P = exp(logP);
        Ttempline(j,i) = T;
        set(fluid,'T',T,'P',P);
        htempline(j,i) = enthalpy_mass(fluid);
        Ptempline(j,i) = pressure(fluid);
        i = i+1;
    end
end

% Add isotherms above the critical temperature.
for j=j+1:1:length(Tlist)
    fluid = Solution('liquidvapor.xml','carbondioxide');
    T = Tlist(j)

    logdP = (log(Pmax) - log(Pmin))/150;
    i = 1;
    for logP=log(Pmax):-logdP:log(Pmin)  % Stop short of saturation.
        P = exp(logP);
        Ttempline(j,i) = T;
        set(fluid,'T',T,'P',P);
        htempline(j,i) = enthalpy_mass(fluid);
        Ptempline(j,i) = pressure(fluid);
        i = i+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start a set of isochores.
set(fluid,'T',Tmin,'P',Pmin);
vmin = 1/density(fluid);
vmax = vgt;
set(fluid,'T',Tmax,'V',vmin);
%Pmax = pressure(fluid);
vlist = [0.00125 .0025 .005 0.01 .02 .05 0.1];
for j=1:1:length(vlist)
    v = vlist(j)
    dT = (Tmax - Tmin)/150;
    i = 1;
    for T=Tmin:dT:Tmax  % Stop short of saturation.
        Tvoluline(j,i) = T;
        set(fluid,'T',T,'V',v);
        hvoluline(j,i) = enthalpy_mass(fluid);
        Pvoluline(j,i) = pressure(fluid);
        i = i+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start a set of isentropes.
slist     = [sft/1e3 2.75 3 3.25 3.5 3.75 4 4.25 4.5]*1e3;
lPlow  = log(Pmin);
lPhigh = log(Pmax);
for j=1:1:length(slist)
    s = slist(j)
    ldP = (lPhigh - lPlow)/150;
    i = 1;
    for lP=lPlow:ldP:lPhigh  % Stop short of saturation.
        P = exp(lP);
        try
            setState_SP(fluid,[s,P]);
            Tentrline(j,i) = temperature(fluid);
            hentrline(j,i) = enthalpy_mass(fluid);
            Pentrline(j,i) = P;
        catch
            disp('Trouble finding isentrope...')
            fluid = importPhase('liquidvapor.xml','carbondioxide');
            Tentrline(j,i) = Tentrline(j,i-1);
            hentrline(j,i) = hentrline(j,i-1);
            Pentrline(j,i) = Pentrline(j,i-1);
        end
        i = i+1;
    end
end

% Make the plots.
figure(1)
clf
i = 1;
semilogy(htempline(i,:)/1e3,Ptempline(i,:)/1e6,'r')
hold on
semilogy(hvoluline(i,:)/1e3,Pvoluline(i,:)/1e6,'Color',[0 .6 0])
semilogy(hentrline(i,:)/1e3,Pentrline(i,:)/1e6,'Color','b')
semilogy(hc/1e3,Pc/1e6,'kd')
semilogy([hft hgt]/1e3,[Pt Pt]/1e6,'ko--')
for i=1:1:length(Tlist)
    semilogy(htempline(i,:)/1e3,Ptempline(i,:)/1e6,'r')
end
for i=1:1:length(vlist)
    semilogy(hvoluline(i,:)/1e3,Pvoluline(i,:)/1e6,'Color',[0 .6 0])
end
for i=1:1:length(slist)
    semilogy(hentrline(i,:)/1e3,Pentrline(i,:)/1e6,'Color','b')
end
semilogy(hliqline/1e3,Pliqline/1e6,'k')
semilogy(hvapline/1e3,Pvapline/1e6,'k')
semilogy([hft hgt]/1e3,[Pt Pt]/1e6,'ko--')
hold off
xlabel('Specific Enthalpy (kJ/kg)')
ylabel('Pressure (MPa)')
for i=1:1:length(Tlist)
    text(htempline(i,150)/1e3,Ptempline(i,150)/1e6,...
        num2str(Tlist(i)),'Color','r')
end
for i=1:1:length(vlist)
    text(hvoluline(i,150)/1e3,Pvoluline(i,150)/1e6,...
        num2str(vlist(i)),'Color',[0 .6 0])
end
for i=1:1:length(slist)
    text(hentrline(i,150)/1000,Pentrline(i,150)/1e6,...
        num2str(slist(i)/1e3),'Color','b')
end
legend('Isotherm (K)','Isochore (m^3/kg)','Isentrope (kJ/kg-K)',...
    'Critical Point','Triple Line')
title('Carbon Dioxide')
scale = axis;
axis([scale(1) -8600 .4 50]);
plotfixer

figure(2)
clf
hold on
i = 1;
plot(hpresline(i,:)/1000,Tpresline(i,:),'b')
plot(hvoluline(i,:)/1000,Tvoluline(i,:),'-','Color',[0 .6 0])
plot(hentrline(i,:)/1000,Tentrline(i,:),'r')
plot(hc/1000,Tc,'kd')
plot([hft/1000 hgt/1000],[Tt Tt],'ko--')
for i=1:1:length(Plist)
    plot(hpresline(i,:)/1000,Tpresline(i,:),'b')
end
for i=1:1:length(vlist)
    plot(hvoluline(i,:)/1000,Tvoluline(i,:),'-','Color',[0 .6 0])
end
for i=1:1:length(slist)
    plot(hentrline(i,:)/1000,Tentrline(i,:),'r')
end
plot(hliqline/1000,Tsatline,'k')
plot(hvapline/1000,Tsatline,'k')
plot([hft/1000 hgt/1000],[Tt Tt],'ko--')
hold off
xlabel('Specific Enthalpy (kJ/kg-K)')
ylabel('Temperature (K)')
for i=1:1:length(Plist)
    text(hpresline(i,150/2)/1000,Tpresline(i,150/2),num2str(Plist(i)),'Color','b')
end
for i=1:1:length(vlist)
    text(hvoluline(i,150/3)/1000,Tvoluline(i,150/3),num2str(vlist(i)),'Color',[0 .6 0])
end
for i=1:1:length(slist)
    text(hentrline(i,1)/1000,Tentrline(i,1),num2str(slist(i)/1e3),'Color','r')
end
legend('Isobar (bar)','Isochore (m^3/kg)','Isentrope (kJ/kg-K)','Critical Point','Triple Line')
title('Carbon Dioxide')
scale = axis;
axis([scale(1) -8700 200 500])
plotfixer
