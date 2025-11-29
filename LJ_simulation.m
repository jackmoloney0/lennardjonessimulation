% ----------------------------------------------------
% Script Name:
% Author: Jack Moloney
% Date Created 20-Nov-2025
% Purpose: Lennard-Jones Potential Simulation
% ----------------------------------------------------

% Load & Extract the Lennard-Jones parameters from lj_param.dat
A = load('lj_param.dat');
eps_Ar    = A (1,1);
sigma_Ar  = A (1,2);
eps_CH    = A (2,1);
sigma_CH  = A (2,2);

eptot = [eps_Ar, eps_CH];
stot = [sigma_Ar, sigma_CH];

% Range of r Value in Angstroms
r = linspace (2, 10, 1000); 

% Calculate the Lennard-Jones potential for Argon & Methane
E_Ar = 4*eps_Ar* ((sigma_Ar./r).^12- (sigma_Ar./r).^6);
E_CH = 4*eps_CH* ((sigma_CH./r).^12- (sigma_CH./r).^6);

% Plot the results
figure;
plot (r, E_Ar, 'b', 'LineWidth', 3.5);
hold on;
plot (r, E_CH, 'r', 'LineWidth', 3.5);
% Axis Limits
ylim ([-200 200]); xlim ([3 10]); 

% Add labels & formatting 
xlabel ('\bf r_{ij} [Å]','FontSize', 16.5, 'Fontname', 'Cambria Math');
ylabel (' \bf E_{ij} [kJ mol^{-1}]','FontSize', 16.5, 'Fontname', 'Cambria Math');
title ('\bf Lennard-Jones Model: Ar–Ar vs CH_4–CH_4','FontSize', 17, 'Fontname', 'Cambria Math');
legend ('Ar', 'CH_4','FontSize', 14, 'Fontname', 'Cambria Math'); 
grid on;

% Request Input for Simulation
N = input ('\nEnter the Number of Particles (Positive Integer):\n\n');
L = input ('\nEnter length of Cube Side (Positive):\n\n ');

% Create Empty Matrix for Particle Vector Position
xyz = zeros (N, 3);
% Initialize Counter for Particles Inserted
icount = 0;

% Loop to Insert Particles
for i = 1:N
    xyz (i,1) = L*rand (1)-0.5*L; % X
    xyz (i,2) = L*rand (1)-0.5*L; % Y
    xyz (i,3) = L*rand (1)-0.5*L; % Z
    icount = icount + 1;
end
fprintf ('\nIntial Positions of Particle: %fi %fj %fk \n', xyz);
fprintf ('\nInserted %+1.0f Particles into the Cube Plane.\n\n', icount);

% Randomly assign N particles as Argon(1) / Methane(2)
species = randi([1 2], N, 1);
% Initialize Total Energy
etot = 0;   

% Calculate Total Lennard-Jones Energy 
for particle1 = 1:N-1
    for particle2 = particle1+1:N 
        dx = xyz(particle1, 1) - xyz(particle2, 1); % Subtract X Components
        dy = xyz(particle1, 2) - xyz(particle2, 2); % Subtract Y Components
        dz = xyz(particle1, 3) - xyz(particle2, 3); % Subtract Z Components
        rij = sqrt(dx^2 + dy^2 + dz^2); % Find Straight Line Distance
        species_particle1 = species(particle1); % Assign Argon(1) / Methane(2)
        species_particle2 = species(particle2); % Assign Argon(1) / Methane(2)
        % If same, use same Argon/Methane Epsilon & Sigma Values
        if species_particle1 ==  species_particle2 
            epsilon_ij = eptot(species_particle1);
            sigma_ij = stot(species_particle1);
        else
            % If species different, i.e. Ar-CH4, Use Lorentz-Berthelot
            % to Find Combined Epsilon & Sigma Values
            epsilon_ij = sqrt(eptot(species_particle1)* eptot(species_particle2));
            sigma_ij = (stot(species_particle1) + stot(species_particle2))/2;
        end
        % Calculate Energy
        eij = 4 * epsilon_ij* ((sigma_ij / rij)^12 - (sigma_ij / rij)^6);
        % Keep Total Updated
        etot = etot + eij;
    end
end

fprintf('\nTotal Initial Energy of the system: %f kJ/mol\n\n', etot);

% Input from User N & L 
nmove = input ('\nEnter number of moves: \n\n');
delta = input ('\nEnter maximum displacement: \n\n');

% Loop for each Displacement Move
for move = 1:nmove
    selec_particle = randi(N); % Radnom Integer [0,N]    
    displacement = delta*(rand(1,3) - 0.5); % Random Displacement
    % Applying Displacement to each Vector Component
    xyz(selec_particle, :) = xyz(selec_particle, :) + displacement; 
   % Pasted Previous Energy Function
    etot = 0;
   for particle1 = 1:N-1
    for particle2 = particle1+1:N 
        dx = xyz(particle1, 1) - xyz(particle2, 1); 
        dy = xyz(particle1, 2) - xyz(particle2, 2); 
        dz = xyz(particle1, 3) - xyz(particle2, 3); 
        rij = sqrt(dx^2 + dy^2 + dz^2); 
        species_particle1 = species(particle1);
        species_particle2 = species(particle2);
        if species_particle1 ==  species_particle2 
            epsilon_ij = eptot(species_particle1);
            sigma_ij = stot(species_particle1);
        else
            epsilon_ij = sqrt(eptot(species_particle1)* eptot(species_particle2));
            sigma_ij = (stot(species_particle1) + stot(species_particle2))/2;
        end
        eij = 4 * epsilon_ij* ((sigma_ij / rij)^12 - (sigma_ij / rij)^6);
        etot = etot + eij;
    end
   end 
   % Display Move No. & Energy up to Move
    fprintf('\nMove %d: Energy Total = %f kJ/mol\n\n', move, etot);
end

% Create a Writing File to Input Results
fid = fopen('final_positions2.txt','w');
fprintf(fid, 'Species\t                x\ty\tz\n\n');
species_names = {'Argon', 'Methane'}; % Assign Argon to 1, Methane to 2
for i = 1:N
    % Loop to Print Species Type & Co-Ordinates for N Particles
    fprintf(fid, '%s\t%f\t%f\t%f\n', species_names{species(i)}, xyz(i,1), xyz(i,2), xyz(i,3));
end
% Total Overall Lennard-Jones Energy
fprintf(fid, '\nTotal Energy: %f kJ/mol\n\n', etot);
fclose(fid);

% Formatting & Structuring 3D Plot
figure;
plot3(xyz(:,1), xyz(:,2), xyz(:,3), 'o','MarkerSize', 8, 'MarkerFaceColor','b');
xlabel('\bf x [Å]', 'FontSize', 16, 'FontName','Cambria Math');
ylabel('\bf y [Å]', 'FontSize', 16, 'FontName','Cambria Math');
zlabel('\bf z [Å]', 'FontSize', 16, 'FontName','Cambria Math');
grid on; axis equal; box on;
% Limiting Axis to Constraints
xlim([-L/2, L/2]); ylim([-L/2, L/2]); zlim([-L/2, L/2]);
title('\bf Final Particle Positions','FontSize', 18, 'FontName','Cambria Math');




