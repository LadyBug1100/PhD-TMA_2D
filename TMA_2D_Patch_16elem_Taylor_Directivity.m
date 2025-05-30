clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time-modulated Array (TMA), Uniform Spacing
% 2D Patch Array 4x4 (16 elements) on YZ plane
% Taylor Coefficients (reshaped to 2D)
%
% Calculates Radiation Pattern (Directivity) for central frequency and harmonics 
% and Excitation Coefficients (for central frequency m=0 and harmonics m=1:3) 
% for ANSYS Simulation
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

% dy = lambda;    % Distance between the elements [m] in y-axis
% dz = lambda;    % Distance between the elements [m] in z-axis

dy = 3*lambda;    % Distance between the elements [m] in y-axis
dz = 3*lambda;    % Distance between the elements [m] in z-axis

%% Steering angle

th0 = 90; % Scanning theta angle [deg]
ph0 = 0; % Scanning phi angle [deg]

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

t1n= zeros(Nz,Ny);
% for nz=1:Nz
%     for ny=1:Ny
%         t1n(nz, ny) = 0;
%     end
% end

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
for m = 1:M+1

    AF{m} = zeros(N_Theta,N_Phi);
    
    for nz=1:Nz
        for ny=1:Ny

            amn(nz,ny,m) = (tau_n(nz,ny)/Tp) * sinc((m-1) * (tau_n(nz,ny)/Tp)) * exp(-j*pi*(m-1)*( (2*t1n(nz,ny) + tau_n(nz,ny))/Tp)); % Coefficients
            temp = exp(j*Phase(nz,ny)*pi/180) * amn(nz,ny,m) * exp(j*k*( dz_vec(nz)*cosd(theta_mesh) + dy_vec(ny)*sind(theta_mesh).*sind(phi_mesh) ) ); 
            AF{m} = AF{m} + temp;
    

        end
    end

end

fprintf('Coefficients for ANSYS\n')
amn_ANSYS = amn; % z-rows, y-columns
amn_ANSYS_MAG = zeros(Nz, Ny, M+1);
amn_ANSYS_PHASE = zeros(Nz, Ny, M+1);
amn_MAG = zeros(Nz, Ny, M+1);
amn_PHASE = zeros(Nz, Ny, M+1);
Pin = zeros(1, M + 1);
for m = 1:M+1
    
    fprintf('Harmonics m = %d\n',m-1)
    amn_ANSYS_MAG(:,:,m) = abs(amn_ANSYS(:,:,m).^2); % non uniform amplitude array (power in [W]) for ANSYS
    amn_ANSYS_PHASE(:,:,m) = angle(amn_ANSYS(:,:,m))*(180/pi);

    amn_MAG(:,:,m) = abs(amn(:,:,m)); 
    amn_PHASE(:,:,m) = angle(amn(:,:,m))*(180/pi);

    fprintf('MAGNITUDE\n')
    amn_ANSYS_MAG(:,:,m)
    fprintf('Phase\n')
    amn_ANSYS_PHASE(:,:,m)

    Pin(m) = sum(sum(amn_ANSYS_MAG(:,:,m)));
end

Pin

% max_AF = max(max(AF{1}))

%% Calculation of Normalized Radiation Intensity

U_Cut = (1/(2*120*pi)) * E_Cut.^2;   

max_U = max(max(U_Cut))

%% Calculation of Normalized Radiation Pattern of the Array


F_Cut = cell(1, M+1);
max_F = zeros(1, M + 1);
for m = 1:M+1

    F_Cut{m} = U_Cut.*(abs(AF{m}).^2);

    max_F(m) = max(max(F_Cut{m}));

end

max_F

%% Calculation of Total Radiated Power

Prad = zeros(1, M + 1);
for m = 1:M+1
    temp = F_Cut{m}.*sind(theta_mesh);
    Prad(m) = sum(temp(:))*(pi/N_Theta)*(2*pi/N_Phi);
end
Prad
Prad_total = sum(Prad)

%% Calculation of Directivity of the Array

Dir = cell(1, M+1);
Dir_dBi = cell(1, M+1);
for m = 1:M+1

    Dir{m} = 4*pi*F_Cut{m}/Prad_total;  % Directivity
    % Dir{m} = 4*pi*F_Cut{m}/Prad(m);  % Directivity


    % Dir_norm = abs(Dir{m})/max(abs(Dir{m}(:))); % Normalized Directivity

    % Dir_dBi{m} = 10*log10(Dir_norm); % Normalized Directivity in [dBi]
    Dir_dBi{m} = 10*log10(Dir{m}); % Directivity in [dBi]

end

%% Find Azimuth and Elevation Cuts

[theta_CUT, ~ ] = find(theta_mesh == th0);
theta_CUT = theta_CUT(1);

[~, phi_CUT] = find(phi_mesh == ph0);
phi_CUT = phi_CUT(1);

Azimuth_Cut = cell(1, M+1);
Elevation_Cut = cell(1, M+1);
for m = 1:M+1

    Azimuth_Cut{m} = Dir_dBi{m}(theta_CUT,:);
    Elevation_Cut{m} = Dir_dBi{m}(:,phi_CUT);

end

%% Directivity - Steered angle

Dir_STEERANGLE = zeros(1, M + 1);
for m = 1:M+1
    Dir_STEERANGLE(m) = Dir_dBi{m}(theta_CUT,phi_CUT); % at (th0, ph0) angle
end

Dir_STEERANGLE

%% Plot Directivity Pattern

txt =  ['2D Patch Array YZ, Nelem = ' num2str(Nz) 'x' num2str(Ny) ', d = ' num2str(dy/lambda) '\lambda, {\theta}_0 = ' num2str(th0) '\circ, {\phi}_0 = '  num2str(ph0) '\circ'];


% Directivity, Elevation cut 
Legend_txt = cell(1, M+1);
figure()
for m = 1:M+1

    Legend_txt{m} = strcat('m= ', num2str((m-1)));

    plot(THETA, Elevation_Cut{m}, 'linewidth', 2); hold on
    title('Radiation Pattern, Directivity, Elevation Cut (\phi=0\circ)')
    subtitle(txt)
    xlabel('\theta [deg]')
    ylabel('Dir [dB]')
    xlim([0 180])
    xticks(0:10:180)
    ylim([-30 20])
    ax = gca;
    ax.FontSize=18;
    grid on
end
legend(Legend_txt)

% Directivity, Azimuth cut 
figure();
for m = 1:M+1
    plot(Phi, Azimuth_Cut{m}, 'linewidth', 2); hold on
    title('Radiation Pattern, Directivity, Azimuth Cut (\theta=90\circ)')
    subtitle(txt)
    xlabel('\phi [deg]')
    ylabel('Dir [dB]')
    xlim([-180 180])
    xticks(-180:20:180)
    ylim([-30 20])
    ax = gca;
    ax.FontSize=18;
    grid on
end
legend(Legend_txt)




