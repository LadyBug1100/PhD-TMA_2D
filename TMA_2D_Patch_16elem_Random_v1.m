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

dz_vec_2D = reshape(dz_vec, [Ny,Nz]);
dy_vec_2D = reshape(dy_vec, [Ny,Nz]);

%% Calculate minimum element spacing - 2D
min2D_col_t = zeros(Ny-1,Nz);
for nz=1:Nz
    for ny = 1:Ny-1

        min2D_col_t(ny,nz) = sqrt( (dz_vec_2D(ny,nz) - dz_vec_2D(ny+1,nz))^2 + (dy_vec_2D(ny,nz) - dy_vec_2D(ny+1,nz))^2 );

    end
end
min2D_col = min(min2D_col_t(:));

min2D_row_t = zeros(Ny,Nz-1);
for nz=1:Nz-1
    for ny = 1:Ny


        min2D_row_t(ny,nz) = sqrt( (dz_vec_2D(ny,nz) - dz_vec_2D(ny,nz+1))^2 + (dy_vec_2D(ny,nz) - dy_vec_2D(ny,nz+1))^2 );

    end
end
min2D_row = min(min2D_row_t(:));

min2D = min([min2D_col min2D_row]);
min2D_lambda = min2D/lambda;

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

t1n= zeros(1,Ny*Nz); % starting time is 0 for all elements

% % Generate Random start time between [0 Tp]
% t1n= [];
% for n = 1:Ny*Nz
%     temp = Tp*rand(1);               
%     t1n = [t1n temp];
% end

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

%% ANSYS Coefficients and input power Pin

% fprintf('Coefficients for ANSYS\n') % z-rows, y-columns
amn_ANSYS_MAG = zeros(Nz, Ny, M+1);
amn_ANSYS_PHASE = zeros(Nz, Ny, M+1);
Pin = zeros(1, M + 1);
for m = 1:M+1
    
%     fprintf('Harmonics m = %d\n',m-1)
    amn_ANSYS_MAG(:,:,m) = abs(amn(:,:,m).^2); % non uniform amplitude array (power in [W]) for ANSYS
    amn_ANSYS_PHASE(:,:,m) = angle(amn(:,:,m))*(180/pi);

%     fprintf('MAGNITUDE\n')
%     amn_ANSYS_MAG(:,:,m)
%     fprintf('Phase\n')
%     amn_ANSYS_PHASE(:,:,m)

    Pin(m) = sum(sum(amn_ANSYS_MAG(:,:,m)));
end

% Pin
Pin_total = sum(Pin);


tEnd = toc(tStart);
fprintf('%d minutes and %f seconds\n', floor(tEnd/60), rem(tEnd,60));

%% Calculation of Normalized Radiation Intensity

U_Cut = (1/(2*120*pi)) * E_Cut.^2;   

% max_U = max(max(U_Cut))

%% Calculation of Normalized Radiation Pattern of the Array

F_Cut = zeros(N_Theta, N_Phi, M+1);
for m = 1:M+1

    F_Cut(:,:,m) = U_Cut.*(abs(AF(:,:,m)).^2);

end

%% Calculation of Total Radiated Power

Prad = zeros(1, M + 1);
for m = 1:M+1
    temp = F_Cut(:,:,m).*sind(theta_mesh);
    Prad(m) = sum(temp(:))*(pi/N_Theta)*(2*pi/N_Phi);
end
% Prad
Prad_total = sum(Prad);

Prad_perc = zeros(1, M + 1); % radiated power distribution (percentage per frequency)
for m = 1:M+1
    Prad_perc(m) = Prad(m)/Prad_total *100;
end

%% Calculation of Directivity of the Array

Dir = zeros(N_Theta, N_Phi, M+1);
Dir_dBi = zeros(N_Theta, N_Phi, M+1);
for m = 1:M+1

    Dir(:,:,m) = 4*pi*F_Cut(:,:,m)/Prad_total;  % Directivity
    % Dir(:,:,m) = 4*pi*F_Cut(:,:,m)/Prad(m);  % Directivity

    temp1 = Dir(:,:,1);
    Dir_norm = abs(Dir(:,:,m))/max(abs(temp1(:)));

    Dir_dBi(:,:,m) = 10*log10(Dir_norm); % Directivity in [dBi]

end


%% Find Azimuth and Elevation Cuts

[theta_CUT, ~ ] = find(theta_mesh == th0);
theta_CUT = theta_CUT(1);

[~, phi_CUT] = find(phi_mesh == ph0);
phi_CUT = phi_CUT(1);


Dir_Azimuth_Cut_dB = zeros(N_Phi, M+1);
Dir_Elevation_Cut_dB = zeros(N_Theta, M+1);
Dir_Azimuth_Cut_lin = zeros(N_Phi, M+1);
Dir_Elevation_Cut_lin = zeros(N_Theta, M+1);
for m = 1:M+1

    Dir_Azimuth_Cut_dB(:,m) = Dir_dBi(theta_CUT,:,m);
    Dir_Elevation_Cut_dB(:,m) = Dir_dBi(:,phi_CUT,m);

    Dir_Azimuth_Cut_lin(:,m) = Dir(theta_CUT,:,m);
    Dir_Elevation_Cut_lin(:,m) = Dir(:,phi_CUT,m);

end

%% Calculate Parameters 

%% Maximum Value at Steered angle for m=0 (centre frequency)


Dir_MAX = 10*log10(Dir(theta_CUT,phi_CUT, 1)); % at (th0, ph0) angle


%% Calculate HPBW (half-power beamwidth) and Sidelobe Level (SLL) for m=0 (centre frequency) and Sideband levels (SBL) for m=1...M (harmonics)

SBL = NaN(1, M);
for m = 1:M+1

    if m ==1 % centre frequency

        % Azimuth
        [peaks_az, locs_az, w_az] = findpeaks(Dir_Azimuth_Cut_lin(:,m),'SortStr','descend' ,'WidthReference','halfheight');
        HPBW_az = w_az(locs_az==phi_CUT)*(Phi(2)-Phi(1));
        locs_az(locs_az==phi_CUT)=[];
        SLL_az = Dir_Azimuth_Cut_dB(locs_az(1), m);
        
        % Elevation
        [peaks_el, locs_el, w_el] = findpeaks(Dir_Elevation_Cut_lin(:,m),'SortStr','descend','WidthReference','halfheight');
        HPBW_el = w_el(locs_el==theta_CUT)*(THETA(2)-THETA(1));
        locs_el(locs_el==theta_CUT)=[];
        SLL_el = Dir_Elevation_Cut_dB(locs_el(1), m);
        
        SLL = max([ SLL_az SLL_el]); % Sidelobe Level 
        HPBW = max([ HPBW_az HPBW_el]); % Beamwidth

    else % harmonics

        % Azimuth
        [peaks_az, locs_az] = findpeaks(Dir_Azimuth_Cut_lin(:,m),'SortStr','descend');
        SBL_az = Dir_Azimuth_Cut_dB(locs_az(1), m);
       
        % Elevation
        [peaks_el, locs_el, w_el] = findpeaks(Dir_Elevation_Cut_lin(:,m),'SortStr','descend');
        SBL_el = Dir_Elevation_Cut_dB(locs_el(1), m);

        SBL(m) = max([ SBL_az SBL_el]);

    end
end
% SBL

% Prad_total/Pin_total *100

%% Summarize Results

fprintf('\nArray Results:\n')
fprintf('\n')
for m = 1:M+1
    fprintf([strcat('Pin_',num2str(m-1)) ': %.2f [W] '], Pin(m))
end
fprintf('\nPin TOTAL: %.2f [W]\n', Pin_total)
for m = 1:M+1
    fprintf([strcat('Prad_',num2str(m-1)) ': %.2f [W] (%.2f %%) '], Prad(m), Prad_perc(m))
end
fprintf('\nPrad TOTAL: %.2f [W]\n', Prad_total)
fprintf('\ndmin: %.2f [%c]\n', min2D_lambda, 955)
fprintf('HPBW: %.2f [deg]\n', HPBW)
fprintf('MAX Dir: %.2f [dB]\n', Dir_MAX)
fprintf('SLL_0: %.2f [dB]\n', SLL)
for m = 2:M+1
    fprintf([strcat('SBL_',num2str(m-1)) ': %.2f [dB] '], SBL(m))
end
fprintf('\n')


%% Plot Directivity Pattern

txt =  ['2D Patch Array YZ, N = ' num2str(Nz) 'x' num2str(Ny) ', d0 = ' num2str(d0y/lambda) '\lambda, Box_{yz} = ' num2str(Box_y/lambda) '\lambda, ({\theta}_0, {\phi}_0) = (' num2str(th0) '\circ, '  num2str(ph0) '\circ)'];

%% Array Distribution - Elements Positions
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

%% Switching sequences
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

%% Radiated Power Distribution

figure()
bar(0:M, Prad_perc,'g')
title('Radiated Power')
ylabel('Power [%]')
xlabel('Harmonic, m')
ylim([0 100])
yticks(0:10:100)
ax = gca;
ax.FontSize=16;
grid on


%% Directivity, Elevation cut 
Legend_txt = cell(1, M+1);
figure('Position',[200 100 800 450])
for m = 1:M+1

    Legend_txt{m} = strcat('m= ', num2str((m-1)));

    plot(THETA, Dir_Elevation_Cut_dB(:,m), 'linewidth', 2); hold on
    title('Radiation Pattern, Directivity, Elevation Cut (\phi=0\circ)')
    subtitle(txt)
    xlabel('\theta [deg]')
    ylabel('Dir [dB]')
    xlim([0 180])
    xticks(0:10:180)
    ylim([-50 0])
    ax = gca;
    ax.FontSize=18;
    grid on
end
xline([THETA(theta_CUT) - HPBW_el/2 THETA(theta_CUT) + HPBW_el/2],'black','LineStyle','--','LineWidth',1); hold on
yline(-3,'-','BW','FontSize',12, 'Color', 'r','LineStyle','--','LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment','center','LineWidth',1); hold on
yline(SLL,'-','SLL_0','FontSize',12, 'Color', 'black','LineStyle','--','LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment','left','LineWidth',1); hold on 
for m = 2:M+1
    yline(SBL(m),'-',strcat('SBL_',num2str(m-1)),'FontSize',12, 'Color', 'black','LineStyle','--','LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment','left','LineWidth',1);
end
legend(Legend_txt)

%% Directivity, Azimuth cut 

figure('Position',[1000 100 800 450]);
for m = 1:M+1
    plot(Phi, Dir_Azimuth_Cut_dB(:,m), 'linewidth', 2); hold on
    title('Radiation Pattern, Directivity, Azimuth Cut (\theta=90\circ)')
    subtitle(txt)
    xlabel('\phi [deg]')
    ylabel('Dir [dB]')
    xlim([-180 180])
    xticks(-180:20:180)
    ylim([-50 0])
    ax = gca;
    ax.FontSize=18;
    grid on
end
xline([Phi(phi_CUT) - HPBW_az/2 Phi(phi_CUT) + HPBW_az/2],'black','LineStyle','--','LineWidth',1); hold on
yline(-3,'-','BW','FontSize',12, 'Color', 'r','LineStyle','--','LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment','center','LineWidth',1); hold on
yline(SLL,'-','SLL_0','FontSize',12, 'Color', 'black','LineStyle','--','LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment','left','LineWidth',1); hold on
for m = 2:M+1
    yline(SBL(m),'-',strcat('SBL_',num2str(m-1)),'FontSize',12, 'Color', 'black','LineStyle','--','LabelVerticalAlignment','bottom', 'LabelHorizontalAlignment','left','LineWidth',1);
end
legend(Legend_txt)



