clear all

%% Set up Air properties
N = 3;
Setup_Air_Props
Mair = 28.9586;                 % kg/kmol
Tmaxcondentherm = 132.6312;     % K
Pmaxcondentherm = 3.78502e6;    % Pa
rmaxcondentherm = 10.4477*Mair; % kg/m3
Tmaxcondenbar   = 132.6035;     % K
Pmaxcondenbar   = 3.7891e6;     % Pa
rmaxcondenbar   = 11.0948*Mair; % kg/m3
Tcritair        = 132.5306;     % K
Pcritair        = 3.7860e6;     % Pa
rcritair        = 11.8308*Mair; % kg/m3
% Bottom of dome conditions for air:
Tsolidair       = 59.75;        % K
Psolidair       = 5265;         % Pa
rsolidair       = 33.067*Mair;  % kg/m3
% Upper limit conditions for the model:
Tupper          = 870;          % K
Pupper          = 100e6;        % Pa
rupper          = rsolidair;
% Lower limit to stay just above solid air:
Tlower          = 60;           % K
% Molar masses
MW_O2 = 31.999; % (g/mol)
MW_N2 = 28.0134; % (g/mol)
MW_Ar = 39.948; % (g/mol)

Setup_Air_Props

% Functions for saturated liquids and vapors
% [T rl rv y] = Fast_Bubble_cP(x,P,varargin)
% [T rv rl x] = Fast_Dew_cP(y,P,varargin)
%         Tstart  = varargin{1};
%         xstart  = varargin{2};
%         rlstart = varargin{3};
%         rvstart = varargin{4};
% Functions for multi-component flash processes under the vapor dome
% [T x y rf rg] = Flash_zqP(z,q,P,varargin)
% [T q V y x rg rf] = Flash_zhP(z,h,P,varargin)
% [q V y x rg rf] = Flash_zTP(z,T,P,varargin)
%         T  = varargin{1};
%         y  = varargin{2};
%         rf = varargin{3};
%         rg = varargin{4};

% Sate 1: Air at ambient conditions
c(N2) = 0.79; c(O2) = 0.21; % Binary Air
T1 = 298.15;
P1 = oneatm;
rv1 = rv_cTP(c,T1, P1);
s1 = s_crT(c, rv1, T1);

% State 2: Air compressed isothermally to 200 bar
T2 = T1;
P2 = 200e5;
rv2 = rv_cTP(c,T2, P2);
s2 = s_crT(c, rv2, T2);

% State 3: Air Cooled isobarically to 130k
T3 = 130;
P3 = P2;
rv3 = rv_cTP(c,T3, P3);
s3 = s_crT(c, rv3, T3);
h3 = h_crT(c, rv3, T3);

% State 4: Air flashed to 1 bar
P4 = 1e5;
[T4, q4, V4, y4, x4, rg4, rf4] = Flash_zhP(c, h3, P4);

%% Now iterate over output composition until column top matches air flash

% Column Parameters
n_trays = 3;
reboil_quality = 0.6;
P_col = 1e5;

% NR Parameters
imax = 10; i = 0;
x_inc = 1e-3;
err_tol = 1e-3;

% Make a guess
x_out(N2) = 0.9;
x_out(O2) = 1 - x_out(N2);

while true
    i = i + 1;
    
    [~, vapout(1), liqin(1)] = Reboiler_xqP(x_out,reboil_quality,P_col);
    for i = 2:n_trays+1
        [vapout(i), liqin(i)] = upTray(vapout(i-1), liqin(i-1));
    end
    err = abs( liqin(end).c(N2) - x4(N2));

    % Do NR
    % High
    x_high(N2) = x_out(N2) + x_inc;
    x_high(O2) = 1 - x_high(N2);
    [Q, vapout_high(1), liqin_high(1)] = Reboiler_xqP(x_high,reboil_quality,P_col);
    for i = 2:n_trays+1
        [vapout_high(i), liqin_high(i)] = upTray(vapout_high(i-1), liqin_high(i-1));
    end
    err_high = abs( liqin_high(end).c(N2) - x4(N2));
    % Low
    x_low(N2) = x_out(N2) - x_inc;
    x_low(O2) = 1 - x_low(N2);
    [Q, vapout_low(1), liqin_low(1)] = Reboiler_xqP(x_low,reboil_quality,P_col);
    for i = 2:n_trays+1
        [vapout_low(i), liqin_low(i)] = upTray(vapout_low(i-1), liqin_low(i-1));
    end
    err_low = abs( liqin_low(end).c(N2) - x4(N2));

    % Increment x_out
    dedx = (err_high - err_low)/(2*x_inc);
    x_out(N2) = x_out(N2) - err/dedx;
    x_out(O2) = 1 - x_out(N2);

    if i>imax
        disp("Failed to Converge Column")
        break
    elseif err < err_tol
        disp("Successfully solved for outlet composition")
    end

end
