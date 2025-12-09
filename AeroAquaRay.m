%% Aero-Aqua Ray: Simulation
clear; clc; close all;

%% 1. System Parameters
m_robot = 2.0;          
g = 9.81;               
rho_water = 1000;
rho_air = 1.225;

% Areas
S_wing = 0.6;           % Planform Area (for Lift)
S_frontal = 0.06;       % Frontal Area (for Water Drag)

% Coefficients
Cd_water = 0.3;         
Kt_swim = 150.0;        
Kt_burst = 1000.0;      

% Glider Polar
Cl_slope = 2*pi;        
AoA_trim = deg2rad(6);  
Cd0_air = 0.02;         
k_ind = 0.04;           

% CPG Parameters
alpha = 5.0;
omega_swim = 2*pi*2.0;
omega_burst = 2*pi*4.0; 
mu_swim = 1.0;
mu_glide = -1.0;

%% 2. Simulation Setup
T_total = 9;           
dt = 0.002;             
time = 0:dt:T_total;
N = length(time);

state = zeros(6, N);
state(:, 1) = [0; -4.0; 2.0; 0; 0.1; 0]; 

% --- Enhanced Data Logging ---
medium_log = zeros(1, N);
pitch_log = zeros(1, N); 
gamma_log = zeros(1, N); % Flight Path Angle
aoa_log = zeros(1, N);   % Angle of Attack
mode_log = zeros(1, N); 

% Force Logging: [Thrust; Drag; Lift; Net_X; Net_Z]
force_log = zeros(5, N); 

current_pitch = 0; 

%% 3. Main Loop
for k = 1:N-1
    x = state(1, k);  z = state(2, k);
    vx = state(3, k); vz = state(4, k);
    v_sq = vx^2 + vz^2;
    v_mag = sqrt(v_sq);
    gamma = atan2(vz, vx);
    
    % --- Control Logic ---
    if z < 0
        % === WATER PHASE ===
        in_water = true;
        mode = 1;
        rho = rho_water;
        
        if z > -2.0
            % [Breach]
            target_pitch = deg2rad(60); 
            target_mu = mu_swim * 2.0;
            target_omega = omega_burst;
            thrust_gain = Kt_burst;
        else
            % [Cruise]
            target_pitch = deg2rad(15);
            target_mu = mu_swim;
            target_omega = omega_swim;
            thrust_gain = Kt_swim;
        end
        aoa = 0; % AoA undefined/irrelevant in simple water model
        
    else
        % === AIR PHASE ===
        in_water = false;
        rho = rho_air;
        target_mu = mu_glide; 
        thrust_gain = 0; 
        target_omega = 0;
        
        if vz > 0.5 
            mode = 2; % Ascent
            target_pitch = gamma; 
        elseif vz > -1.0 && vz <= 0.5
            mode = 3; % Pushover
            target_pitch = deg2rad(-2); 
        else
            mode = 4; % Glide
            if v_mag < 4.0
                target_pitch = deg2rad(-15);
            else
                target_pitch = gamma + AoA_trim;
                target_pitch = min(target_pitch, deg2rad(20));
            end
        end
        aoa = current_pitch - gamma;
    end
    
    % Update Logs
    medium_log(k) = in_water;
    mode_log(k) = mode;
    gamma_log(k) = gamma;
    aoa_log(k) = aoa;
    
    % Smoothing
    if in_water, smooth = 0.1; else, smooth = 0.05; end
    current_pitch = (1 - smooth) * current_pitch + smooth * target_pitch;
    pitch_log(k) = current_pitch;
    
    % CPG Step
    cpg_x = state(5, k); cpg_y = state(6, k);
    r2 = cpg_x^2 + cpg_y^2;
    state(5, k+1) = cpg_x + (alpha*(target_mu - r2)*cpg_x - target_omega*cpg_y)*dt;
    state(6, k+1) = cpg_y + (alpha*(target_mu - r2)*cpg_y + target_omega*cpg_x)*dt;
    wing_action = state(5, k+1);
    
    % Force Calculation
    if in_water
        T = thrust_gain * (wing_action^2);
        D = 0.5 * rho * v_sq * S_frontal * Cd_water; % S_frontal
        L = 0; % Negligible hydrodynamic lift in this model
        B = m_robot * g * 1.05; 
        
        Fx = T*cos(current_pitch) - D*(vx/v_mag);
        Fz = T*sin(current_pitch) + B - m_robot*g - D*(vz/v_mag);
    else
        % Aerodynamics
        if abs(aoa) < deg2rad(20)
            Cl = Cl_slope * aoa;
        else
            Cl = Cl_slope * aoa * 0.5; % Stall
        end
        Cd = Cd0_air + k_ind * Cl^2;
        
        T = 0;
        L = 0.5 * rho * v_sq * S_wing * Cl; % S_wing
        D = 0.5 * rho * v_sq * S_wing * Cd;
        B = 0;
        
        Fx = -D*cos(gamma) - L*sin(gamma);
        Fz = L*cos(gamma) - D*sin(gamma) - m_robot*g;
    end
    
    force_log(:, k) = [T; D; L; Fx; Fz];
    
    % Integration
    state(3, k+1) = vx + (Fx/m_robot)*dt;
    state(4, k+1) = vz + (Fz/m_robot)*dt;
    state(1, k+1) = x + state(3, k+1)*dt;
    state(2, k+1) = z + state(4, k+1)*dt;
end
% Finalize logs
pitch_log(N) = pitch_log(N-1); mode_log(N) = mode_log(N-1);
gamma_log(N) = gamma_log(N-1); aoa_log(N) = aoa_log(N-1);

%% 4. Advanced Visualization (Publication Ready)

% Color Palette (Scientific)
col_water = [0.2 0.6 0.9];
col_air   = [0.9 0.9 0.95];
col_mode  = [0 0.447 0.741; 0.85 0.325 0.098; 0.929 0.694 0.125; 0 0 0];

% Indices
idx_exit = find(medium_log == 0, 1);
t_exit = time(idx_exit);

% --- Figure 1: The Trajectory (Macro View) ---
figure(1); set(gcf, 'Color', 'w', 'Position', [100, 100, 900, 500]);
hold on;
% Environment
fill([-5, max(state(1,:))+5, max(state(1,:))+5, -5], [-10, -10, 0, 0], ...
     col_water, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
fill([-5, max(state(1,:))+5, max(state(1,:))+5, -5], [0, 0, 10, 10], ...
     col_air, 'EdgeColor', 'none', 'FaceAlpha', 0.2);
yline(0, 'b-', 'LineWidth', 1.5); 

% Plot by segment
scatter(state(1,:), state(2,:), 15, mode_log, 'filled');
colormap(col_mode);
c = colorbar; c.Ticks = [1.4 2.1 3 3.8]; 
c.TickLabels = {'Swim', 'Ballistic Ascent', 'Pushover', 'Glide'};

% Annotations
text(state(1,idx_exit), 0.5, '\leftarrow Breach', 'FontSize', 10);
[z_max, idx_max] = max(state(2,:));
text(state(1,idx_max), z_max+0.5, sprintf('Apex: %.1fm', z_max), 'HorizontalAlignment', 'center');

%title('Fig 1: Multimodal Trans-medium Trajectory');
xlabel('Distance (m)'); ylabel('Altitude (m)');
axis equal; grid on; ylim([-5, z_max+2]); xlim([-2, state(1,end)+2]);
set(gca, 'FontSize', 10, 'LineWidth', 1.2);
exportgraphics(gcf, 'trajectory.pdf', 'ContentType', 'vector');

% --- Figure 2: Kinematics (Velocity) ---
figure(2); set(gcf, 'Color', 'w', 'Position', [150, 150, 600, 600]);
subplot(2,1,1);
v_total = sqrt(state(3,:).^2 + state(4,:).^2);
plot(time, v_total, 'k-', 'LineWidth', 2); hold on;
xline(t_exit, 'b--', 'Surface Breach');
ylabel('Total Speed (m/s)'); grid on;
%title('Fig 2A: Velocity Profile');
legend('Speed', 'Breach Event');

subplot(2,1,2);
plot(time, state(3,:), 'r', 'LineWidth', 1.5); hold on;
plot(time, state(4,:), 'b', 'LineWidth', 1.5);
xline(t_exit, 'k--');
yline(0, 'k-', 'Alpha', 0.2);
ylabel('Components (m/s)'); xlabel('Time (s)'); grid on;
legend('Vx (Horizontal)', 'Vz (Vertical)');
%title('Fig 2B: Velocity Decomposition');
exportgraphics(gcf, 'velocity.pdf', 'ContentType', 'vector');

% --- Figure 3: Aerodynamics (Pitch & AoA) ---
figure(3); set(gcf, 'Color', 'w', 'Position', [200, 200, 600, 600]);
subplot(2,1,1);
plot(time, rad2deg(pitch_log), 'k', 'LineWidth', 2); hold on;
plot(time, rad2deg(gamma_log), 'g--', 'LineWidth', 1.5);
xline(t_exit, 'b--');
ylabel('Angle (deg)'); 
%title('Fig 3A: Attitude vs Flight Path');
legend('Pitch (\theta)', 'Flight Path (\gamma)'); grid on;
ylim([-30, 80]);

subplot(2,1,2);
% Filter AoA plot to only show Air phase clearly
plot(time(idx_exit:end), rad2deg(aoa_log(idx_exit:end)), 'm', 'LineWidth', 2); hold on;
yline(6, 'k--', 'Trim AoA (Target)');
yline(15, 'r--', 'Stall Limit');
ylabel('Angle of Attack (deg)'); xlabel('Time (s)');
%title('Fig 3B: Angle of Attack (Air Phase Only)');
grid on; ylim([-5, 20]);
text(time(end)-2, 7, 'Stable Glide Region');
exportgraphics(gcf, 'attitude_angle.pdf', 'ContentType', 'vector');

% --- Figure 4: Forces (Dynamics) ---
figure(4); set(gcf, 'Color', 'w', 'Position', [250, 250, 600, 500]);
% Plot Thrust, Drag, Lift
plot(time, force_log(1,:), 'r', 'LineWidth', 1.5); hold on;
plot(time, force_log(2,:), 'b', 'LineWidth', 1.5);
plot(time, force_log(3,:), 'g', 'LineWidth', 1.5);
xline(t_exit, 'k--', 'Breach');

ylabel('Force (N)'); xlabel('Time (s)');
%title('Fig 4: Dynamic Force Distribution');
legend('Thrust (Propulsion)', 'Drag (Resistance)', 'Lift (Aerodynamic)');
grid on; 
xlim([time(idx_exit)-2, time(end)]); % Zoom in around breach
ylim([-10, 200]); % Crop the huge underwater burst spike to show air details
text(t_exit+0.2, 100, 'Drag Drop ->', 'Color', 'b');
exportgraphics(gcf, 'force.pdf', 'ContentType', 'vector');

% --- Figure 5: CPG Phase Portrait ---
figure(5); set(gcf, 'Color', 'w', 'Position', [300, 300, 500, 500]);
plot(state(5,:), state(6,:), 'k'); hold on;
% Mark start and end
plot(state(5,1), state(6,1), 'go', 'MarkerFaceColor', 'g');
plot(state(5,end), state(6,end), 'rs', 'MarkerFaceColor', 'r');
%title('Fig 5: CPG Phase Portrait (Limit Cycle -> Fixed Point)');
xlabel('Neuron State u_1'); ylabel('Neuron State u_2');
axis equal; grid on;
annotation('textarrow', [0.7 0.5], [0.7 0.5], 'String', 'Convergence');
exportgraphics(gcf, 'cpg.pdf', 'ContentType', 'vector');
