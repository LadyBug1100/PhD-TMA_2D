clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time-modulated Array (TMA), Uniform Spacing
% 2D Patch Array 4x4 (16 elements) on YZ plane
% Taylor Coefficients (reshaped to 2D)
%
% Beam steering of the first harmonic m=1 to (th0_1, ph0_1). Utilizing 't1n' - switch
% start parameter)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Load Electric field of a single element

Ns = 1800;

THETA = ((0:Ns)*(180/Ns));            % Theta angles [deg]
N_Theta = length(THETA);            % Number of Theta points

Phi = ((0:Ns)*(360/Ns))-180;              % Phi angles [deg]
N_Phi = length(Phi);                % Number of Phi points

[phi_mesh, theta_mesh] = meshgrid(Phi, THETA);

path = '\\file\Usersa$\ave68\Home\My Documents\Research_Anastasia\____________PhD\________AntennaArray\__2__Non_Uniform_Amplitude_Uniform_Spacing\';
load(strcat(path,'E_cut.mat'));

%% 1D Geometry Array Inputs

Ny = 4; % number of elements in y-axis
Nz = 4; % number of elements in z-axis

f = 2.5e9; % Frequency [Hz]

lambda = 3e8/f;    % Wavelength [m]
k = 2*pi/lambda;   % Wavenumber [rad/m]

dy = lambda;    % Distance between the elements [m] in y-axis
dz = lambda;    % Distance between the elements [m] in z-axis

%% Steering angle

th0 = 90; % Scanning theta angle [deg]
ph0 = 0; % Scanning phi angle [deg]

% m = 1 (first harmonic)
th0_1 = 80; % Steering theta angle [deg] of first harmonic (m=1)
ph0_1 = 10; % Steering phi angle [deg] of first harmonic (m=1)

%% Spacing distances

dz_vec = [];
for n = 1 : Nz
    dz_vec = [dz_vec (n-1)*dz]; % uniform spacing array
end

dy_vec = [];
for n = 1 : Ny
    dy_vec = [dy_vec (n-1)*dy]; % uniform spacing array
end

%% Time-modulated array Inputs

Tp = 1; % normalized modulation (switching) period
M = 3; % number of harmonics to plot (m = 0, 1, ... M)

%% Switch on-time (duration) for all elements (should be normalized from 0 to 1)

% Amplitude Coefficients - Taylor
load(strcat(path,'Coeff.mat')); % upload array coefficients generated in "Taylor_Coeff_v2.m" ("save_ON" variable in "Taylor_Coeff_2.m" is 1)
tau_n = Coeff;     % non uniform amplitude Taylor array (current in [V])

tau_n = reshape(tau_n,[Nz,Ny])';
tau_n = [tau_n(1:Nz/2,:); flip(tau_n((Nz/2)+1:end,:),2)];


%% Switch start time t1n for all elements

% t1n= zeros(Nz,Ny);
% % for nz=1:Nz
% %     for ny=1:Ny
% %         t1n(nz, ny) = 0;
% %     end
% % end

t1n = zeros(Nz,Ny);
t1n_check = zeros(Nz,Ny);
for nz=1:Nz
    for ny=1:Ny

        temp =  (k/(2*pi)) * ( dz_vec(nz)*cosd(th0_1) + dy_vec(ny)*sind(th0_1).*sind(ph0_1) ); 

        temp_check = mod(temp, 1);
    
        t1n(nz,ny) = temp;  % start time is zero (for simplicity)

        t1n_check(nz,ny) = temp_check;
    
    end
end
t1n
t1n = t1n_check
diff(t1n) % between rows
diff(t1n,1,2) % between columns

% t1n_check = zeros(Nz,Ny);
% for nz=1:Nz
%     for ny=1:Ny
% 
%         t1n_check(nz,ny) = mod(t1n(nz,ny), 1);
% 
%     end
% end
% 
% t1n_check
% diff(t1n_check) % between rows
% diff(t1n_check,1,2) % between columns

%% Phase Distribution for Scanning

for nz=1:Nz
    for ny=1:Ny
        Data_Phase(nz, ny) = exp(-j*k*( dz_vec(nz)*cosd(th0) + dy_vec(ny)*sind(th0).*sind(ph0) ) );
    end
end
 
Phase = angle(Data_Phase)*180/pi;  % Phase distribution [deg]

%% Linear Array Factor for M harmonics

amn = zeros(Nz, Ny, M+1);
AF = cell(1, M+1);
AF_dB = cell(1, M+1);
for m = 1:M+1

    AF{m} = zeros(N_Theta,N_Phi);
    
    for nz=1:Nz
        for ny=1:Ny

            amn(nz,ny,m) = (tau_n(nz,ny)/Tp) * sinc((m-1) * (tau_n(nz,ny)/Tp)) * exp(-j*pi*(m-1)*( (2*t1n(nz,ny) + tau_n(nz,ny))/Tp)); % Coefficients
            temp = exp(j*Phase(nz,ny)*pi/180) * amn(nz,ny,m) * exp(j*k*( dz_vec(nz)*cosd(theta_mesh) + dy_vec(ny)*sind(theta_mesh).*sind(phi_mesh) ) ); 
            AF{m} = AF{m} + temp;
    
            AF_norm = abs(AF{m})/max(abs(AF{1}(:)));
            AF_dB{m} = 20*log10(AF_norm); % Array Factor in [dB] 

        end
    end

end


%% Find Azimuth and Elevation Cuts

% [theta_CUT, ~ ] = find(theta_mesh == th0);
% theta_CUT = theta_CUT(1);
% 
% [~, phi_CUT] = find(phi_mesh == ph0);
% phi_CUT = phi_CUT(1);

[theta_CUT, ~ ] = find(theta_mesh == th0_1);
theta_CUT = theta_CUT(1);

[~, phi_CUT] = find(phi_mesh == ph0_1);
phi_CUT = phi_CUT(1);

AF_Azimuth_Cut = cell(1, M+1);
AF_Elevation_Cut = cell(1, M+1);
for m = 1:M+1

    AF_Azimuth_Cut{m} = AF_dB{m}(theta_CUT,:);
    AF_Elevation_Cut{m} = AF_dB{m}(:,phi_CUT);


end

%% Plot Directivity Pattern

txt =  ['2D Patch Array YZ, Nelem = ' num2str(Nz) 'x' num2str(Ny) ', d = ' num2str(dy/lambda) '\lambda, ({\theta}_1, {\phi}_1) = (' num2str(th0_1) '\circ, '  num2str(ph0_1) '\circ)'];

% figure()
% bar(1:Nz,tau_n)
% title(['Switching sequence of ' num2str(Nz) ' element array'])
% ylabel('Normalized ON-time')
% xlabel('Element Number')
% xticks(1:1:Nz)
% ax = gca;
% ax.FontSize=16;
% grid on


% Array Factor, Elevation cut 
Legend_txt = cell(1, M+1);
figure()
for m = 1:M+1

    Legend_txt{m} = strcat('m= ', num2str((m-1)));

    plot(THETA, AF_Elevation_Cut{m}, 'linewidth', 2); hold on
    title({['Normalized Radiation Pattern, Array Factor, Elevation Cut (\phi=' num2str(ph0_1) '\circ)']})
    subtitle(txt)
    xlabel('\theta [deg]')
    ylabel('AF [dB]')
    xlim([0 180])
    xticks(0:10:180)
    ylim([-50 0])
    ax = gca;
    ax.FontSize=18;
    grid on
end
legend(Legend_txt)


% Array Factor, Azimuth cut 
figure();
for m = 1:M+1
    plot(Phi, AF_Azimuth_Cut{m}, 'linewidth', 2); hold on
    title({['Normalized Radiation Pattern, Array Factor, Azimuth Cut (\theta=' num2str(th0_1) '\circ)']})
    subtitle(txt)
    xlabel('\phi [deg]')
    ylabel('AF [dB]')
    xlim([-180 180])
    xticks(-180:20:180)
    ylim([-50 0])
    % yticks(-30:5:10)
    ax = gca;
    ax.FontSize=18;
    grid on
end
legend(Legend_txt)

% figure()
% plot(THETA, AF_Elevation_Cut{2}, 'linewidth', 2, 'Color', "#D95319"); hold on
% title({['Normalized Radiation Pattern, Array Factor m=1, Elevation Cut (\phi=' num2str(ph0_1) '\circ)']})
% subtitle(txt)
% xlabel('\theta [deg]')
% ylabel('AF [dB]')
% xlim([0 180])
% xticks(0:10:180)
% ylim([-50 0])
% ax = gca;
% ax.FontSize=18;
% grid on
% 
% figure()
% plot(Phi, AF_Azimuth_Cut{2}, 'linewidth', 2, 'Color', "#D95319"); hold on
% title({['Normalized Radiation Pattern, Array Factor m=1, Azimuth Cut (\theta=' num2str(th0_1) '\circ)']})
% subtitle(txt)
% xlabel('\phi [deg]')
% ylabel('AF [dB]')
% xlim([-180 180])
% xticks(-180:20:180)
% ylim([-50 0])
% ax = gca;
% ax.FontSize=18;
% grid on



