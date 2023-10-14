% Written by Vinit Doke - October 2023

% WAKE UP CYCLING (PULSE_MAX - READ - PULSE_MIN READ - PULSE_MAX - READ -----)

%TODO : Implement Function for Constant Amplitude Pulse Train Testing
%TODO : Implement Function for Varying Amplitude Testing

tic

clear all; close all; clc; % clear all variables, close all figures, clear command window
% loadlibrary('wgfmu', 'wgfmu.h'); 

% UNITS
ms = 1e-3;
us = 1e-6;
ns = 1e-9;

% DEVICE CONSTANTS
GATE_LENGTH = 500; % nm
GATE_WIDTH = 80; % nm

% INSTRUMENT SETTINGS
I_RANGE_WRITE = 6002;
I_RANGE_READ  = 6003;

%% TIMING CONSTANTS :==
% PULSE TIMING PARAMETERS
t_rise_pulse  = 1*us;
t_fall_pulse  = 1*us;
t_width_pulse = 1*us;
t_delay_pulse = 10*us;

% READ TIMING PARAMETERS
t_read       = 100*us;
t_rise_read  = 1*us;
t_fall_read  = 1*us;
t_delay_read = 10*us;


% WHAT IS THIS?
It0 = 1e-6*(GATE_WIDTH/GATE_LENGTH);

% MEASUREMENT PARAMS: $ MAJOR CHANGES OCCUR HERE ---------------
REPETITIONS = 1;
TYPE = 2;
% TYPE = 2;
% TYPE = 3;

if TYPE == 1
    % TYPE 1 : RESET + SET Sequence Then Read
    V_MIN   = -4.5; % Volts
    V_MAX   =  4.5; % Volts
    N_PAIRS = 10;

    V_AMPLITUDES = zeros(1, 2*N_PAIRS);
    V_AMPLITUDES(1:2:end) = V_MIN;
    V_AMPLITUDES(2:2:end) = V_MAX;

    disp(V_AMPLITUDES);

elseif TYPE == 2
    % TYPE 2 : RESET + N Pulses of V_N Amplitude Then Read
    RESET_AMPLITUDE = -4.5; % Volts
    N_PULSES = 10;
    PULSE_AMPLITUDE = 2.5; % Volts

    V_AMPLITUDES = zeros(1, N_PULSES + 1);
    V_AMPLITUDES(1) = RESET_AMPLITUDE;
    V_AMPLITUDES(2:end) = PULSE_AMPLITUDE;

    disp(V_AMPLITUDES);
end

% VOLTAGE SWEEP RANGE
VG_READ_MIN = -0.5; % Volts
VG_READ_MAX =  2.0; % Volts

% DRAIN VOLTAGE DURING READ
VD_READ = 0.1; % Volts


% CREATE VOLTAGE AND TIME SCHEME FOR THE PULSES
v_switch = PulseTrain(V_AMPLITUDES);
t_sw     = [0 t_delay_pulse t_rise_pulse/10 t_width_pulse t_fall_pulse/10 t_delay_pulse];


tw = [];
vw = [];

for i = 1:REPETITIONS
    tw = [tw t_sw]; %FIX
    vw = [vw v_switch]; %FIX
end

vw = vw'; %CHECK IF REQUIRED
tw = tw'; %CHECK IF REQUIRED
t_finalw = cumsum(tw);


sample_interval_write = 1*us;
sample_points_write   = round(t_finalw(end)/sample_interval_write);


% CREATE VOLTAGE AND TIME SCHEME FOR THE READ SWEEP

% FOR GATE
vg_read1 = VG_READ_MIN*ones(1, length(V_AMPLITUDES));
vg_read2 = VG_READ_MAX*ones(1, length(V_AMPLITUDES));

for i=1:length(vg_read1)
    vg_read(i, :) = ([0 0 vg_read1(i) vg_read2(i) 0 0]);
end

tg_read = [0 t_delay_read t_rise_read t_read t_fall_read t_delay_read];

no_of_points = round((t_read + 2*(t_rise_pulse + t_delay_read))/us);
sample_points_read = 1000; % can bring this line to top as constant param


% FOR DRAIN
vd_read  = ([0 0 VD_READ VD_READ 0 0]);

tg  = tg_read;
vg  = vg_read;
vd  = vd_read;

vd = vd';
vg = vg';
tg = tg';

t_final = cumsum(tg);

sample_interval_read = t_final(end)/sample_points_read;

if sample_interval_read <10*ns
    sample_interval_read=10*ns;
    sample_points_read=t_final(end)/sample_interval_read;
end

pt0 = round(sample_points_read*t_final(1)/t_final(end));
pt1 = round(sample_points_read*t_final(2)/t_final(end));
pt2 = round(sample_points_read*t_final(3)/t_final(end));
pt3 = round(sample_points_read*t_final(4)/t_final(end));
pt4 = round(sample_points_read*t_final(5)/t_final(end));
pt5 = round(sample_points_read*t_final(6)/t_final(end));


m=0; % MOVE THIS

OC = length(V_AMPLITUDES);
color = distinguishable_colors(OC);

% FILENAME CREATION
filename = fileManager('wakeup', 1)


full_data_write = [];
full_data_read  = [];
VtD = [];

for i = 1:OC
    i
    Bw = [tw vw(:,i)];
    csvwrite('voltage_pattern_write.csv', Bw);

    data_out_write(:,:,i) = func_write(sample_points_write, sample_interval_write, I_RANGE_WRITE);

    V_gate(:,i) = data_out_write(:,2,i);
    I_drain(:,i) = abs(data_out_write(:,4,i));
    I_source(:,i) =abs(data_out_write(:,6,i));
    I_sub(:,i) =abs( data_out_write(:,8,i));
    T(:,i) = data_out_write(:,1,i);

    full_data_write = [full_data_write data_out_write(:,:,i)];


    % READ
    Br = [tg vg(:,i)];
    csvwrite('voltage_pattern_read_gate.csv',Br);
    Ar=[tg vd];
    csvwrite('voltage_pattern_read_drain.csv',Ar);

    data_out_read(:,:,i) = func_WGFMU(sample_points_read, sample_interval_read, I_RANGE_READ);

    Vr_gate(:,i) = data_out_read(:,2,i);
    Ir_drain(:,i) = abs(smooth(smooth(smooth(data_out_read(:,4,i),'sgolay',4))));
    Ir_source(:,i) =abs(smooth(smooth(smooth( data_out_read(:,6,i),'sgolay',4))));
    Ir_sub(:,i) = abs(smooth(smooth( data_out_read(:,8,i),'sgolay',4)));
    Tr(:,i) = data_out_read(:,1,i);

    full_data_read = [full_data_read data_out_read(:,:,i)];


    % PLOTTING

    figure(1) % WRITE VOLTAGE v/s TIME
    plot(T(:,i), V_gate(:,i), 'color', color(i,:), 'linewidth', 2);
    xlim([T(1,i) T(end, i)]);

    title('Write Voltage', 'FontSize', 16, 'FontWeight', 'b');
    xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'b');
    ylabel('Voltage (V)', 'FontSize', 16, 'FontWeight', 'b');

    ax = gca;
    set(ax, 'FontSize', 16, 'FontWeight', 'b');
    grid on;
    hold on;

    ax2 = gca;
    fig = gcf;
    set(fig, 'Color', [1,1,1]);
    grid on;
    hold on;


    figure(2) % SUBSTRATE CURRENT v/s TIME
    fig2 = semilogy(Tr(:,i), Ir_sub(:,i), 'color', color(i,:), 'linewidth', 2);
    hold on;
    xlim([Tr(1,i) Tr(end, i)]);

    title('Substrate Current', 'FontSize', 16, 'FontWeight', 'b');
    xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'b');
    ylabel('Current (A)', 'FontSize', 16, 'FontWeight', 'b');

    ax = gca;
    set(ax, 'FontSize', 16, 'FontWeight', 'b');
    hold on;
    grid on;

    fig = gcf;
    set(fig, 'Color', [1,1,1]);
    grid on;
    box on;

    figname = [filename '_substrate_current'];
    hgsave(fig2, figname);
    saveas(fig2, figname, 'png');


    figure(3) % DRAIN CURRENT v/s TIME
    fig3 = semilogy(Tr(:,i), Ir_drain(:,i), 'color', color(i,:), 'linewidth', 2);
    hold on;
    xlim([Tr(1,i) Tr(end, i)]);

    title('Drain Current', 'FontSize', 16, 'FontWeight', 'b');
    xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'b');
    ylabel('Current (A)', 'FontSize', 16, 'FontWeight', 'b');

    ax = gca;
    set(ax, 'FontSize', 16, 'FontWeight', 'b');
    hold on;
    grid on;

    fig = gcf;
    set(fig, 'Color', [1,1,1]);
    grid on;
    box on;

    figname = [filename '_drain_current'];
    hgsave(fig3, figname);
    saveas(fig3, figname, 'png');


    figure(4) % SOURCE CURRENT v/s TIME
    fig4 = semilogy(Tr(:,i), Ir_source(:,i), 'color', color(i,:), 'linewidth', 2);
    hold on;
    xlim([Tr(1,i) Tr(end, i)]);

    title('Source Current', 'FontSize', 16, 'FontWeight', 'b');
    xlabel('Time (s)', 'FontSize', 16, 'FontWeight', 'b');
    ylabel('Current (A)', 'FontSize', 16, 'FontWeight', 'b');

    ax = gca;
    set(ax, 'FontSize', 16, 'FontWeight', 'b');
    hold on;
    grid on;

    fig = gcf;
    set(fig, 'Color', [1,1,1]);
    grid on;
    box on;

    figname = [filename '_source_current'];
    hgsave(fig4, figname);
    saveas(fig4, figname, 'png');



    figure(5) % ID v/s VG
    pshift = 10;
    fig5 = semilogy(Vr_gate(pt2 + pshift : pt3 - pshift, i), Ir_source(pt2 + pshift : pt3 - pshift, i), 'color', color(i,:), 'linewidth', 2);
    hold on;

    title('ID v/s VG', 'FontSize', 16, 'FontWeight', 'b');
    xlabel('Gate Voltage (V)', 'FontSize', 16, 'FontWeight', 'b');
    ylabel('Source Current (A)', 'FontSize', 16, 'FontWeight', 'b');

    ax = gca;
    set(ax, 'FontSize', 16, 'FontWeight', 'b');
    hold on;
    grid on;

    fig = gcf;
    set(fig, 'Color', [1,1,1]);
    grid on;
    box on;

    figname = [filename ' Vgate Vs I source-t'];
    hgsave(fig5, figname);
    saveas(fig5, figname, 'png');


    % Vt Extraction
    m = m+1;
    Itd(:,i)= double(Ir_source(pt2+pshift:pt3-pshift,i)>It0);
    vtg = Vr_gate(pt2+pshift:pt3-pshift,i);

    vtd_index(m) = find(Itd(:,i),1,'first');
    VtD = [VtD vtg(vtd_index(m))];
    V_SW_S(m) = v_sw(i);

    figure(6) % Threshold Voltage v/s Switch Voltage
    fig6 = plot(V_SW_S(m), VtD(m), 'o', 'color', color(i,:), 'linewidth', 2);
    hold on;
    
    title('Threshold Voltage v/s Switch Voltage', 'FontSize', 16, 'FontWeight', 'b');
    xlabel('Switch Voltage (V)', 'FontSize', 16, 'FontWeight', 'b');
    ylabel('Threshold Voltage (V)', 'FontSize', 16, 'FontWeight', 'b');

    ax = gca;
    set(ax, 'FontSize', 16, 'FontWeight', 'b');
    hold on;
    grid on;

    fig = gcf;
    set(fig, 'Color', [1,1,1]);
    grid on;
    box on;

    figname = [filename ' Vt Vs Vsw'];
    hgsave(fig6, figname);
    saveas(fig6, figname, 'png');
end

save([filename '.mat'])

csvwrite([filename '_write.csv'], full_data_write);
csvwrite([filename '_read.csv'], full_data_read);
VtD

toc
