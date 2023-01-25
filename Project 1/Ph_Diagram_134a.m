% Make an internal energy-entropy diagram for fluid (T-s too).
% C.F. Edwards, 12/20/09

clear all
format compact
fprintf('\n********************************************************\n')

% Set up a fluid object to work with in Cantera.
fluid = Solution('liquidvapor.cti','hfc134a');

% Get the critical point props.
Tc = critTemperature(fluid)
Pc = critPressure(fluid)
set(fluid,'T',Tc,'P',Pc);
sc = entropy_mass(fluid);
hc = enthalpy_mass(fluid);
vc = 1/density(fluid);

% Set the triple point props.  Use a value epsilon above Tt to get the
% bottom edge of the vapor dome.
Tt = 170
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

% Start a set of isobaric curves.
Plist = [Pmin/1e5 .01 .1 1 10];  % Pressure in bar
for j=1:1:length(Plist)
    fluid = Solution('liquidvapor.xml','hfc134a');
    P = Plist(j)*1e5

    % Do the compressed liquid side.
    setState_Psat(fluid,[P 0]);
    Tsat = temperature(fluid)
    dT = (Tsat - Tmin)/50;
    for i=1:1:50  % Stop short of saturation.
        T = Tmin + (i-1)*dT;
        Tpresline(j,i) = T;
        set(fluid,'T',T,'P',P);
        spresline(j,i) = entropy_mass(fluid);
    end
    i = i+1;
    setState_Psat(fluid,[P 0]);
    Tpresline(j,i) = Tsat;   % Add the saturation point now.
    spresline(j,i) = entropy_mass(fluid);
    
    % Now go across the dome.
    i = i+1;
    dq = 1/50;
    for q=0+dq:dq:1-dq  % Stop short of saturation.
        Tpresline(j,i) = Tsat;
        setState_Psat(fluid,[P q]);
        spresline(j,i) = entropy_mass(fluid);
        i = i+1;
    end
    Tpresline(j,i) = Tsat;   % Add the saturation point now.
    setState_Psat(fluid,[P 1]);
    spresline(j,i) = entropy_mass(fluid);
    
    % Do the vapor side.
    i = i+1;
    dT = (Tmax - Tsat)/50;
    for ii=1:1:50
        T = Tsat + ii*dT;  % Start just above saturation.
        Tpresline(j,i) = T;
        set(fluid,'T',T,'P',P);
        spresline(j,i) = entropy_mass(fluid);
        i = i+1;
    end
end

% Add an isobar above the critical pressure.
P = 100;            % In bar
Plist = [Plist P];  % Add to the list
P = P*100000;       % In Pascals
j = length(Plist);
dT = (Tmax - Tmin)/150;
i = 1;
for T=Tmin:dT:Tmax  % Stop short of saturation.
    Tpresline(j,i) = T;
    set(fluid,'T',T,'P',P);
    spresline(j,i) = entropy_mass(fluid);
    i = i+1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start a set of isotherms.
Tlist = [Tmin 200 250 300 350 400 450];
for j=1:1:3    % Do the part of list below the critical point. 
    fluid = Solution('liquidvapor.xml','hfc134a');
    T = Tlist(j)

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
    T = Tlist(j)

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start a set of isochores.
set(fluid,'T',Tmin,'P',Pmin);
vmin = 1/density(fluid);
vmax = vgt;
set(fluid,'T',Tmax,'V',vmin);
%Pmax = pressure(fluid);
vlist = [vft .0007 .0008 .0009 .001 vc .01 .1 1 10 vgt];
for j=1:1:length(vlist)
    v = vlist(j)
    dT = (Tmax - Tmin)/150;
    i = 1;
    for T=Tmin:dT:Tmax  % Stop short of saturation.
        Tvoluline(j,i) = T;
        set(fluid,'T',T,'V',v);
        svoluline(j,i) = entropy_mass(fluid);
        hvoluline(j,i) = enthalpy_mass(fluid);
        Pvoluline(j,i) = pressure(fluid);
        if(pressure(fluid) > Pmax)
            Tvoluline(j,i) = Tvoluline(j,i-1);
            svoluline(j,i) = svoluline(j,i-1);
            hvoluline(j,i) = hvoluline(j,i-1);
            Pvoluline(j,i) = Pvoluline(j,i-1);
        end
        i = i+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start a set of isenthalps.
hlist     = [-100 -50  0    50   100  150  200  250  300]*1e3;
Plowlist  = [.004 .004 .004 .004 .004 .004 .004 .004 .004]*1e5;  
Phighlist = [30   30   70   100  30   100  100  100  100]*1e5;  
for j=1:1:length(hlist)
    h = hlist(j)
    lPlow  = log(Plowlist(j));
    lPhigh = log(Phighlist(j));
    ldP = (lPhigh - lPlow)/150;
    i = 1;
    for lP=lPlow:ldP:lPhigh  % Stop short of saturation.
        P = exp(lP);
        setState_HP(fluid,[h,P]);
        Tenthline(j,i) = temperature(fluid);
        senthline(j,i) = entropy_mass(fluid);
        if(Tenthline(j,i) > Tmax)
            Tenthline(j,i) = Tenthline(j,i-1);
            senthline(j,i) = senthline(j,i-1);
        end
        i = i+1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Start a set of isentropes.
slist     = [sft   500   1000   1500];
Phighlist = [Pmax  Pmax  Pmax   10e6];
lPlow = log(Pmin);
lPhigh = log(10e6);
for j=1:1:length(slist)
    s = slist(j)
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

% Make the plots.
figure(1)
clf
hold on
i = 1;
plot(spresline(i,:)/1000,Tpresline(i,:),'b')
plot(svoluline(i,:)/1000,Tvoluline(i,:),'-','Color',[0 .6 0])
plot(senthline(i,:)/1000,Tenthline(i,:),'r')
plot(sc/1000,Tc,'kd')
plot([sft/1000 sgt/1000],[Tt Tt],'ko--')
for i=1:1:length(Plist)
    plot(spresline(i,:)/1000,Tpresline(i,:),'b')
end
for i=1:1:length(vlist)
    plot(svoluline(i,:)/1000,Tvoluline(i,:),'-','Color',[0 .6 0])
end
for i=1:1:length(hlist)
    plot(senthline(i,:)/1000,Tenthline(i,:),'r')
end
plot(sliqline/1000,Tsatline,'k')
plot(svapline/1000,Tsatline,'k')
plot([sft/1000 sgt/1000],[Tt Tt],'ko--')
hold off
xlabel('Specific Entropy (kJ/kg-K)')
ylabel('Temperature (K)')
for i=1:1:length(Plist)
    text(spresline(i,150)/1000,Tpresline(i,150),num2str(Plist(i)),'Color','b')
end
for i=1:1:length(vlist)
    text(svoluline(i,150)/1000,Tvoluline(i,150),num2str(vlist(i)),'Color',[0 .6 0])
end
for i=1:1:length(hlist)
    text(senthline(i,1)/1000,Tenthline(i,1),num2str(hlist(i)/1e3),'Color','r')
end
legend('Isobar (bar)','Isochore (m^3/kg)','Isenthalp (kJ/kg)','Critical Point','Triple Line');
title('HFC134a')
scale = axis;
axis([scale(1) scale(2) 150 500])
plotfixer

figure(2)
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
    'Critical Point','Triple Line');
title('HFC134a')
plotfixer
