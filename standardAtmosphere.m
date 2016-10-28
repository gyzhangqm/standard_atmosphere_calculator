%% STANDARD ATMOSPHERE CALCULATOR ===================================

clc
clear
close all

%% User input:
H = 30;                         % altitude [km]
FONT = 10;                      % fontsize for graphs
ftxt = 10;                      % fontsize for zone names
sealevel = [3 76 169] ./ 255;   % sea level conditions colour
target = [241 51 47] ./ 255;    % target conditions colour
zones = [140 148 152] ./ 255;	% atmospheric zones colour
bound = [62 139 14] ./ 255;     % zone boundaries colour
xtext = 185;                    % starting position of the zone name

%% Sea level condititions:
p0 = 101325;    % pressure [Pa]
T0 = 288.15;    % temperature [K]
rho0 = 1.225;   % density [kg/m^3]

%% Constants:
H = H*1000;     % altitude [m]
g = 9.80665;    % gravity constant [m/s^2]
R = 287.057;	% ideal gas constant (dry air) [m^2/(s^2 * K)]
step = 1000;	% step for the graphs

%% Lapse rates at different altitudes:
% TROPOSPHERE .......................................... (0-10.999)km
h_ts = 0;        % [m]
a_ts = -0.0065; % [K/m]
% TROPOPAUSE ========================================== (11-19.999)km
h_tp = 11000;   % [m]
a_tp = 0;       % [K/m] (isothermal)
% STRATOSPHERE ........................................ (20-31.999)km
h_ss1 = 20000;  % [m]
a_ss1 = 0.001;  % [K/m]
% ..................................................... (32-46.999)km
h_ss2 = 32000;  % [m]
a_ss2 = 0.0028; % [K/m]
% STRATOPAUSE ========================================= (47-50.999)km
h_sp = 47000;   % [m]
a_sp = 0;       % [K/m] (isothermal)
% MESOSPHERE .......................................... (51-70.999)km
h_ms1 = 51000;  % [m]
a_ms1 = -0.0028;% [K/m]
% ......................................................... (71-85)km
h_ms2 = 71000;  % [m]
a_ms2 = -0.002; % [K/m]
% ===================================================================
h_fin = 85000;  % [m]

%% Upper boundary parameters:
% Upper boundary of troposphere: ....................................
T_1 = T0 + a_ts * (h_tp - h_ts);
p_1 = p0*(T_1/T0)^(-g/(a_ts*R));
rho_1 = rho0*(T_1/T0)^(-g/(a_ts*R) - 1);
% Graph plotting data:
Y1 = h_ts:step:h_tp;
X1 = T0 + a_ts*(Y1 - h_ts);

% Upper boundary of tropopause: .....................................
T_2 = T_1;
p_2 = p_1 * exp(-(g/(R*T_2)) * (h_ss1 - h_tp));
rho_2 = rho_1 * exp(-(g/(R*T_2)) * (h_ss1 - h_tp));
% Graph plotting data:
Y2 = h_tp:step:h_ss1;
X2 = T_1 + a_tp*(Y2 - h_tp);

% Upper boundary of stratosphere (1): ...............................
T_3 = T_2 + a_ss1 * (h_ss2 - h_ss1);
p_3 = p_2*(T_3/T_2)^(-g/(a_ss1*R));
rho_3 = rho_2*(T_3/T_2)^(-g/(a_ss1*R) - 1);
% Graph plotting data:
Y3 = h_ss1:step:h_ss2;
X3 = T_2 + a_ss1*(Y3 - h_ss1);

% Upper boundary of stratosphere (2): ...............................
T_4 = T_3 + a_ss2 * (h_sp - h_ss2);
p_4 = p_3*(T_4/T_3)^(-g/(a_ss2*R));
rho_4 = rho_3*(T_4/T_3)^(-g/(a_ss2*R) - 1);
% Graph plotting data:
Y4 = h_ss2:step:h_sp;
X4 = T_3 + a_ss2*(Y4 - h_ss2);

% Upper boundary of stratopause: ....................................
T_5 = T_4;
p_5 = p_4 * exp(-(g/(R*T_5)) * (h_ms1 - h_sp));
rho_5 = rho_4 * exp(-(g/(R*T_5)) * (h_ms1 - h_sp));
% Graph plotting data:
Y5 = h_sp:step:h_ms1;
X5 = T_4 + a_sp*(Y5 - h_sp);

% Upper boundary of mezosphere (1): .................................
T_6 = T_5 + a_ms1 * (h_ms2 - h_ms1);
p_6 = p_5*(T_6/T_5)^(-g/(a_ms1*R));
rho_6 = rho_5*(T_6/T_5)^(-g/(a_ms1*R) - 1);
% Graph plotting data:
Y6 = h_ms1:step:h_ms2;
X6 = T_5 + a_ms1*(Y6 - h_ms1);

%% Temperature, pressure and density calculation:
if H > h_fin
    disp(['Altitude specified larger than 85 km.']);
elseif H < 0
    disp(['Altitude cannot be negative!']);
else
    disp(['Calculating for altitude ', num2str(H/1000), ' km.']);

    if H >= h_ts && H < h_tp  
        disp(['You are in the troposphere.']);
        zone = 'troposphere';
        % In the troposphere:
        T_fin = T0 + a_ts * (H - h_ts);
        p_fin = p0*(T_fin/T0)^(-g/(a_ts*R));
        rho_fin = rho0*(T_fin/T0)^(-g/(a_ts*R) - 1);
        
        Yend = h_ts:step:H;
        Xend = T0 + a_ts*(Yend - h_ts);
   
    elseif H >= h_tp && H < h_ss1  
        disp(['You are in the tropopause.']);
        disp(['The temperature is constant in this zone.']);
        zone = 'tropopause';
        % In the tropopause:
        T_fin = T_1;
        p_fin = p_1 * exp(-(g/(R*T_fin)) * (H - h_tp));
        rho_fin = rho_1 * exp(-(g/(R*T_fin)) * (H - h_tp));
        
        plot(X1, Y1/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        hold on
        grid on
        line([180 300], [h_tp/1000 h_tp/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1)
        text(xtext, h_tp/1000+1.5, 'tropopause', 'Color', bound, 'FontSize', ftxt)
        
        Yend = h_tp:step:H;
        Xend = T_1 + a_tp*(Yend - h_tp);
        
    elseif H >= h_ss1 && H < h_ss2
        disp(['You are in the stratosphere (1).']);
        zone = 'stratosphere (1)';
        % In the stratosphere (1):
        T_fin = T_2 + a_ss1 * (H - h_ss1);
        p_fin = p_2*(T_fin/T_2)^(-g/(a_ss1*R));
        rho_fin = rho_2*(T_fin/T_2)^(-g/(a_ss1*R) - 1);    
        
        plot(X1, Y1/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        hold on
        grid on
        plot(X2, Y2/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        line([180 300], [h_tp/1000 h_tp/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1)
        line([180 300], [h_ss1/1000 h_ss1/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1)
        text(xtext, h_tp/1000+1.5, 'tropopause', 'Color', bound, 'FontSize', ftxt)
        text(xtext, h_ss1/1000+1.5, 'stratosphere (1)', 'Color', bound, 'FontSize', ftxt)
        
        Yend = h_ss1:step:H;
        Xend = T_2 + a_ss1*(Yend - h_ss1);
        
    elseif H >= h_ss2 && H < h_sp
        disp(['You are in the stratosphere (2).']);
        zone = 'stratosphere (2)';
        % In the stratosphere (2):
        T_fin = T_3 + a_ss2 * (H - h_ss2);
        p_fin = p_3*(T_fin/T_3)^(-g/(a_ss2*R));
        rho_fin = rho_3*(T_fin/T_3)^(-g/(a_ss2*R) - 1);     

        plot(X1, Y1/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        hold on
        grid on
        plot(X2, Y2/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        plot(X3, Y3/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        line([180 300], [h_tp/1000 h_tp/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1)
        line([180 300], [h_ss1/1000 h_ss1/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1)        
        line([180 300], [h_ss2/1000 h_ss2/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1)  
        text(xtext, h_tp/1000+1.5, 'tropopause', 'Color', bound, 'FontSize', ftxt)
        text(xtext, h_ss1/1000+1.5, 'stratosphere (1)', 'Color', bound, 'FontSize', ftxt)
        text(xtext, h_ss2/1000+1.5, 'stratosphere (2)', 'Color', bound, 'FontSize', ftxt)
        
        Yend = h_ss2:step:H;
        Xend = T_3 + a_ss2*(Yend - h_ss2);
        
    elseif H >= h_sp && H < h_ms1
        disp(['You are in the stratopause.']);
        disp(['The temperature is constant in this zone.']);
        zone = 'stratopause';
        % In the stratopause:
        T_fin = T_4;
        p_fin = p_4 * exp(-(g/(R*T_fin)) * (H - h_sp));
        rho_fin = rho_4 * exp(-(g/(R*T_fin)) * (H - h_sp));

        plot(X1, Y1/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        hold on
        grid on
        plot(X2, Y2/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        plot(X3, Y3/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        plot(X4, Y4/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        line([180 300], [h_tp/1000 h_tp/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1)
        line([180 300], [h_ss1/1000 h_ss1/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1)        
        line([180 300], [h_ss2/1000 h_ss2/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1)          
        line([180 300], [h_sp/1000 h_sp/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1)   
        text(xtext, h_tp/1000+1.5, 'tropopause', 'Color', bound, 'FontSize', ftxt)
        text(xtext, h_ss1/1000+1.5, 'stratosphere (1)', 'Color', bound, 'FontSize', ftxt)
        text(xtext, h_ss2/1000+1.5, 'stratosphere (2)', 'Color', bound, 'FontSize', ftxt)
        text(xtext, h_sp/1000+1.5, 'stratopause', 'Color', bound, 'FontSize', ftxt)
        
        Yend = h_sp:step:H;
        Xend = T_4 + a_sp*(Yend - h_sp);

    elseif H >= h_ms1 && H < h_ms2
        disp(['You are in the mezosphere (1).']);      
        zone = 'mezosphere (1)';
        % In the mezosphere (1):
        T_fin = T_5 + a_ms1 * (H - h_ms1);
        p_fin = p_5*(T_fin/T_5)^(-g/(a_ms1*R));
        rho_fin = rho_5*(T_fin/T_5)^(-g/(a_ms1*R) - 1);     

        plot(X1, Y1/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        hold on
        grid on
        plot(X2, Y2/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        plot(X3, Y3/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        plot(X4, Y4/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        plot(X5, Y5/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        line([180 300], [h_tp/1000 h_tp/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1)
        line([180 300], [h_ss1/1000 h_ss1/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1)        
        line([180 300], [h_ss2/1000 h_ss2/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1)          
        line([180 300], [h_sp/1000 h_sp/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1) 
        line([180 300], [h_ms1/1000 h_ms1/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1)
        text(xtext, h_tp/1000+1.5, 'tropopause', 'Color', bound, 'FontSize', ftxt)
        text(xtext, h_ss1/1000+1.5, 'stratosphere (1)', 'Color', bound, 'FontSize', ftxt)
        text(xtext, h_ss2/1000+1.5, 'stratosphere (2)', 'Color', bound, 'FontSize', ftxt)
        text(xtext, h_sp/1000+1.5, 'stratopause', 'Color', bound, 'FontSize', ftxt)
        text(xtext, h_ms1/1000+1.5, 'mezosphere (1)', 'Color', bound, 'FontSize', ftxt)
        
        Yend = h_ms1:step:H;
        Xend = T_5 + a_ms1*(Yend - h_ms1);

    elseif H >= h_ms2 && H <= h_fin
        disp(['You are in the mezosphere (2).']);
        zone = 'mezosphere (2)';
        % In the mezosphere (2):
        T_fin = T_6 + a_ms2 * (H - h_ms2);
        p_fin = p_6*(T_fin/T_6)^(-g/(a_ms2*R));
        rho_fin = rho_6*(T_fin/T_6)^(-g/(a_ms2*R) - 1);    
        
        plot(X1, Y1/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        hold on
        grid on
        plot(X2, Y2/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        plot(X3, Y3/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        plot(X4, Y4/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        plot(X5, Y5/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        plot(X6, Y6/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
        line([180 300], [h_tp/1000 h_tp/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1)
        line([180 300], [h_ss1/1000 h_ss1/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1)        
        line([180 300], [h_ss2/1000 h_ss2/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1)          
        line([180 300], [h_sp/1000 h_sp/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1) 
        line([180 300], [h_ms1/1000 h_ms1/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1)
        line([180 300], [h_ms2/1000 h_ms2/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1)   
        if H == 85000
            line([180 300], [h_fin/1000 h_fin/1000], 'Color', bound, 'LineStyle', '--', 'LineWidth', 1) 
        end
        text(xtext, h_tp/1000+1.5, 'tropopause', 'Color', bound, 'FontSize', ftxt)
        text(xtext, h_ss1/1000+1.5, 'stratosphere (1)', 'Color', bound, 'FontSize', ftxt)
        text(xtext, h_ss2/1000+1.5, 'stratosphere (2)', 'Color', bound, 'FontSize', ftxt)
        text(xtext, h_sp/1000+1.5, 'stratopause', 'Color', bound, 'FontSize', ftxt)
        text(xtext, h_ms1/1000+1.5, 'mezosphere (1)', 'Color', bound, 'FontSize', ftxt)
        text(xtext, h_ms2/1000+1.5, 'mezosphere (2)', 'Color', bound, 'FontSize', ftxt)

        Yend = h_ms2:step:H;
        Xend = T_6 + a_ms2*(Yend - h_ms2);
        
    end

%% Displaying the results:
disp([' ']);
disp(['Temperature: ', num2str(T_fin), ' K']);
disp(['Pressure: ', num2str(p_fin), ' Pa']);
disp(['Density: ', num2str(rho_fin), ' kg/m^3']);

% Graph for the temperature change until the target altitude:
hfig1 = figure(1);
set(hfig1, 'units','normalized');
set(hfig1, 'Position', [0.5 0.1 0.3 0.8])
plot(Xend, Yend/1000, 'Color', zones, 'LineStyle', '-', 'LineWidth', 2)
hold on
grid on
plot(T_fin, H/1000, 'Color', target, 'LineStyle', '.', 'MarkerSize', 30)
[M] = graphs(FONT);
title(['Altitude = ', num2str(H/1000), ' km, ', zone]);
xlabel('T(h) [K]');
ylabel('altitude [km]')
text(xtext, h_ts/1000+1.5, 'troposphere', 'Color', bound, 'FontSize', ftxt)
xlim([180 300]);
ylim([h_ts/1000 (H+step)/1000]);
daspect([2 1 1])

print('-dpng', '-r300', ['altitude_', num2str(H/1000), '_km.png'])

% Subplot for temperature, pressure and density at the target altitude:
hfig2 = figure(2);
set(gcf, 'units','normalized');
set(hfig2, 'Position', [0.1 0.1 0.3 0.8])
subplot(3,1,1)
line([0, 1], [T0, T0], 'Color', sealevel, 'LineWidth', 2)
hold on
line([0, 1], [T_fin, T_fin], 'Color', target, 'LineWidth', 2)
legend({['T_0 = ', num2str(T0), ' K'], ['T_1 = ', num2str(T_fin), ' K']})
ylim([0 1.01*T0])
[M] = graphs(FONT);
xlabel('T(h) [K]');
title(['Altitude = ', num2str(H/1000), ' km, ', zone]);

subplot(3,1,2)
line([0, 1], [p0/1000, p0/1000], 'Color', sealevel, 'LineWidth', 2)
hold on
line([0, 1], [p_fin/1000, p_fin/1000], 'Color', target, 'LineWidth', 2)
legend({['p_0 = ', num2str(p0/1000), ' kPa'], ['p_1 = ', num2str(p_fin/1000), ' kPa']})
ylim([0 1.01*p0/1000])
[M] = graphs(FONT);
xlabel('p(h) [kPa]');

subplot(3,1,3)
line([0, 1], [rho0, rho0], 'Color', sealevel, 'LineWidth', 2)
hold on
line([0, 1], [rho_fin, rho_fin], 'Color', target, 'LineWidth', 2)
legend({['\rho_0 = ', num2str(rho0), ' kg/m^3'], ['\rho_1 = ', num2str(rho_fin), ' kg/m^3']})
ylim([0 1.01*rho0])
[M] = graphs(FONT);
xlabel('\rho(h) [kg/m^3]');

print('-dpng', '-r300', ['atmosphere_', num2str(H/1000), '_km.png'])

end