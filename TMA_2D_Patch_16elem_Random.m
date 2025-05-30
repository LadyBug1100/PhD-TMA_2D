clc
clear all
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time-modulated Array (TMA), Uniform Spacing
% 2D Patch Array 4x4 (16 elements) on YZ plane
%
% Random coefficients (durations tau_n), starting time is zero for now
%
% Non-uniform Random spacing (uniform distribution with random offsets)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% rng(100);

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

d0y = 3*lambda;    % Distance between the elements [m] in y-axis
d0z = 3*lambda;    % Distance between the elements [m] in z-axis

Box_y = 2*lambda; % Offset boundary [m] in y-direction
Box_z = 2*lambda; % Offset boundary [m] in z-direction

%% Generate Random offsets

% Y offsets
lb_y = -Box_y/2;   % lower Y boundary
ub_y =  Box_y/2;   % upper Y boundary

delta_y_vec = [];
for n = 1:Ny*Nz
    temp = lb_y + (ub_y - lb_y).*rand(1);               
    delta_y_vec = [delta_y_vec temp];
    % delta_y_vec = [delta_y_vec 0];
end

% Z offsets
lb_z = -Box_z/2;   % lower Z boundary
ub_z =  Box_z/2;   % upper Z boundary

delta_z_vec = [];
for n = 1:Ny*Nz
    temp = lb_z + (ub_z - lb_z).*rand(1);               
    delta_z_vec = [delta_z_vec temp];
    % delta_z_vec = [delta_z_vec 0];
end

delta_vec = [delta_z_vec delta_y_vec];


%% Generated dy and dz vectors

dz_vec = [];
dy_vec = [];
n=1;
for nz=1:Nz
    for ny = 1:Ny

        dz = (nz-1)*d0z + delta_vec(n);
        dz_vec = [dz_vec dz];

        dy = (ny-1)*d0y + delta_vec(Nz*Ny + n);
        dy_vec = [dy_vec dy];

        n=n+1;
    end
end
% dz_vec
% dy_vec

%% Steering angle

th0 = 90; % Scanning theta angle [deg]
ph0 = 0; % Scanning phi angle [deg]

%% Time-modulated array Inputs

Tp = 1; % normalized modulation (switching) period
M = 2; % number of harmonics to plot (m = 0, 1, ... M)

%% Switch on-time (duration) for all elements (should be normalized from 0 to 1)

% tau_n= ones(1,Ny*Nz); % duration time is 1 for all elements (ON all the time)

% load(strcat(path,'Coeff.mat')); % upload array coefficients generated in "Taylor_Coeff_v2.m" ("save_ON" variable in "Taylor_Coeff_2.m" is 1)
% tau_n = Coeff;     % non uniform amplitude Taylor array (current in [V])


% Generate Random time durations between [0 Tp]
tau_n= [];
for n = 1:Ny*Nz
    temp = Tp*rand(1);               
    tau_n = [tau_n temp];
end


%% Switch start time t1n for all elements

% t1n= zeros(1,Ny*Nz); % starting time is 0 for all elements

% Generate Random start time between [0 Tp]
t1n= [];
for n = 1:Ny*Nz
    temp = Tp*rand(1);               
    t1n = [t1n temp];
end

time_vec = [tau_n t1n];


%% Variable vectorconsists of 1) Z offsets 2) Y offsets 3) TAU time durations 4) T1n starting time

variable_vec = [delta_vec time_vec];


%% Phase Distribution for Scanning

n = 1;
for nz=1:Nz
    for ny=1:Ny
        
        Data_Phase(nz, ny) = exp(-j*k*( ( (nz-1)*d0z + variable_vec(n) )*cosd(th0) + (  (ny-1)*d0y + variable_vec(Nz*Ny + n) )*sind(th0).*sind(ph0) ) );
        n=n+1;

    end
end
 
Phase = angle(Data_Phase)*180/pi;  % Phase distribution [deg]


%% Array Factor for M harmonics

tStart = tic;

amn = zeros(Nz, Ny, M+1);
AF = zeros(N_Theta, N_Phi, M+1);
AF_temp = zeros(N_Theta,N_Phi);
for m = 1:M+1

    AF_temp = zeros(N_Theta,N_Phi);
    n=1;

    for nz=1:Nz
        for ny=1:Ny

            amn(nz,ny,m) = (variable_vec(2*Nz*Ny + n)/Tp) * sinc((m-1) * (variable_vec(2*Nz*Ny + n)/Tp)) * exp(-j*pi*(m-1)*( (2*variable_vec(3*Nz*Ny + n) + variable_vec(2*Nz*Ny + n))/Tp)); % Coefficients

            temp = exp(j*Phase(nz, ny)*pi/180) * amn(nz,ny,m) * exp(j*k*( ( (nz-1)*d0z + variable_vec(n) )*cosd(theta_mesh) + (  (ny-1)*d0y + variable_vec(Nz*Ny + n) )*sind(theta_mesh).*sind(phi_mesh) ) ); 
           
            AF_temp = AF_temp + temp;

            n=n+1;

        end
    end

    AF(:,:,m)=AF_temp;

end


tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));



% fprintf('Coefficients for ANSYS\n')
% amn_ANSYS = amn; % z-rows, y-columns
% amn_ANSYS_MAG = zeros(Nz, Ny, M+1);
% amn_ANSYS_PHASE = zeros(Nz, Ny, M+1);
% amn_MAG = zeros(Nz, Ny, M+1);
% amn_PHASE = zeros(Nz, Ny, M+1);
% for m = 1:M+1
% 
%     fprintf('Harmonics m = %d\n',m-1)
%     amn_ANSYS_MAG(:,:,m) = abs(amn_ANSYS(:,:,m).^2); % non uniform amplitude array (power in [W]) for ANSYS
%     amn_ANSYS_PHASE(:,:,m) = angle(amn_ANSYS(:,:,m))*(180/pi);
% 
%     amn_MAG(:,:,m) = abs(amn(:,:,m)); 
%     amn_PHASE(:,:,m) = angle(amn(:,:,m))*(180/pi);
% 
%     fprintf('MAGNITUDE\n')
%     amn_ANSYS_MAG(:,:,m)
%     fprintf('Phase\n')
%     amn_ANSYS_PHASE(:,:,m)
% end

%% Calculation of Total Electric Field for M harmonics

E_total = zeros(N_Theta, N_Phi, M+1);
E_total_dB = zeros(N_Theta, N_Phi, M+1);
for m = 1:M+1

    E_total(:,:,m) = E_Cut.*AF(:,:,m);

    E_total_norm = abs(E_total(:,:,m));

    E_total_dB(:,:,m) = 20*log10(E_total_norm); % Total electric field in [dB] 

end


%% Find Azimuth and Elevation Cuts


[theta_CUT, ~ ] = find(theta_mesh == th0);
theta_CUT = theta_CUT(1);

[~, phi_CUT] = find(phi_mesh == ph0);
phi_CUT = phi_CUT(1);


E_total_Azimuth_Cut = zeros(N_Phi, M+1);
E_total_Elevation_Cut = zeros(N_Theta, M+1);
for m = 1:M+1


    E_total_Azimuth_Cut(:,m) = E_total_dB(theta_CUT,:,m);
    E_total_Elevation_Cut(:,m) = E_total_dB(:,phi_CUT,m);

end



%% Plot Directivity Pattern

txt =  ['2D Patch Array YZ, Nelem = ' num2str(Nz) 'x' num2str(Ny) ', dy = dz = ' num2str(d0y/lambda) '\lambda, {\theta}_0 = ' num2str(th0) '\circ, {\phi}_0 = '  num2str(ph0) '\circ'];

% Array Distribution - Elements Positions
figure('Position',[500 250 600 550]);
scatter( dy_vec, dz_vec, 90, 'r','filled'); hold on

rectangle('Position',[-Box_y/2 , -Box_z/2, Box_y, Box_z])
rectangle('Position',[(-Box_y/2 + d0y) , -Box_z/2, Box_y, Box_z])
rectangle('Position',[(-Box_y/2 + 2*d0y) , -Box_z/2, Box_y, Box_z])
rectangle('Position',[(-Box_y/2 + 3*d0y) , -Box_z/2, Box_y, Box_z])

rectangle('Position',[-Box_y/2 , (-Box_z/2 + d0z), Box_y, Box_z])
rectangle('Position',[(-Box_y/2 + d0y) , (-Box_z/2 + d0z), Box_y, Box_z])
rectangle('Position',[(-Box_y/2 + 2*d0y) , (-Box_z/2 + d0z), Box_y, Box_z])
rectangle('Position',[(-Box_y/2 + 3*d0y) , (-Box_z/2 + d0z), Box_y, Box_z])

rectangle('Position',[-Box_y/2 , (-Box_z/2 + 2*d0z), Box_y, Box_z])
rectangle('Position',[(-Box_y/2 + d0y) , (-Box_z/2 + 2*d0z), Box_y, Box_z])
rectangle('Position',[(-Box_y/2 + 2*d0y) , (-Box_z/2 + 2*d0z), Box_y, Box_z])
rectangle('Position',[(-Box_y/2 + 3*d0y) , (-Box_z/2 + 2*d0z), Box_y, Box_z])

rectangle('Position',[-Box_y/2 , (-Box_z/2 + 3*d0z), Box_y, Box_z])
rectangle('Position',[(-Box_y/2 + d0y) , (-Box_z/2 + 3*d0z), Box_y, Box_z])
rectangle('Position',[(-Box_y/2 + 2*d0y) , (-Box_z/2 + 3*d0z), Box_y, Box_z])
rectangle('Position',[(-Box_y/2 + 3*d0y) , (-Box_z/2 + 3*d0z), Box_y, Box_z])

title(['2D Array Distribution (\lambda=' num2str(lambda) 'm)'])
xlabel('y[m]')
ylabel('z[m]')
xlim([-Box_y/2 Box_y/2 + 3*d0y])
ylim([-Box_z/2 Box_z/2 + 3*d0z])
xticks(-Box_y/2:Box_y/2:Box_y/2 + 3*d0y)
yticks(-Box_z/2:Box_z/2:Box_z/2 + 3*d0z)
ax = gca;
ax.FontSize=14;
grid on

% Switching sequences
figure('Position',[500 250 600 550])
for n=1:Ny*Nz

    rectangle('Position', [ n-0.5 t1n(n) 1 tau_n(n)], 'FaceColor', 'b'); hold on

    if (t1n(n) + tau_n(n)) > Tp

        already_ON = Tp - t1n(n);
        remained_time = tau_n(n) - already_ON;

        rectangle('Position', [ n-0.5 0 1 remained_time], 'FaceColor', 'b'); hold on

    end

end
title(['Switching sequence of ' num2str(Ny*Nz) ' element array'])
ylabel('Normalized ON-time')
xlabel('Element Number')
ylim([ 0 Tp])
xlim([1-0.5 Ny*Nz+0.5])
xticks(1:1:Ny*Nz)
yticks(0:0.1:Tp)
ax = gca;
ax.FontSize=16;
grid on


% Total Electric field, Elevation cut 
Legend_txt = cell(1, M+1);
figure()
for m = 1:M+1

    Legend_txt{m} = strcat('m= ', num2str((m-1)));

    plot(THETA, E_total_Elevation_Cut(:,m), 'linewidth', 2); hold on
    title('Radiation Pattern, Total Electric Field, Elevation Cut (\phi=0\circ)')
    subtitle(txt)
    xlabel('\theta [deg]')
    ylabel('E total [dB]')
    xlim([0 180])
    xticks(0:10:180)
    ylim([-20 50])
    ax = gca;
    ax.FontSize=18;
    grid on
end
legend(Legend_txt)

% Total Electric field, Azimuth cut 
figure();
for m = 1:M+1
    plot(Phi, E_total_Azimuth_Cut(:,m), 'linewidth', 2); hold on
    title('Radiation Pattern, Total Electric Field, Azimuth Cut (\theta=90\circ)')
    subtitle(txt)
    xlabel('\phi [deg]')
    ylabel('E total [dB]')
    xlim([-180 180])
    xticks(-180:20:180)
    ylim([-20 50])
    ax = gca;
    ax.FontSize=18;
    grid on
end
legend(Legend_txt)

% return

% % Array Factor, Elevation cut 
% % Legend_txt = cell(1, M+1);
% figure()
% for m = 1:M+1
% 
%     Legend_txt{m} = strcat('m= ', num2str((m-1)));
% 
%     plot(THETA, AF_Elevation_Cut{m}, 'linewidth', 2); hold on
%     title('Normalized Radiation Pattern, Array Factor, Elevation Cut (\phi=0\circ)')
%     subtitle(txt)
%     xlabel('\theta [deg]')
%     ylabel('AF [dB]')
%     xlim([0 180])
%     xticks(0:10:180)
%     ylim([-50 0])
%     ax = gca;
%     ax.FontSize=18;
%     grid on
% end
% legend(Legend_txt)
% 
% 
% % Array Factor, Azimuth cut 
% figure();
% for m = 1:M+1
%     plot(Phi, AF_Azimuth_Cut{m}, 'linewidth', 2); hold on
%     title('Normalized Radiation Pattern, Array Factor, Azimuth Cut (\theta=90\circ)')
%     subtitle(txt)
%     xlabel('\phi [deg]')
%     ylabel('AF [dB]')
%     xlim([-180 180])
%     xticks(-180:20:180)
%     ylim([-50 0])
%     % yticks(-30:5:10)
%     ax = gca;
%     ax.FontSize=18;
%     grid on
% end
% legend(Legend_txt)
% 
% 

