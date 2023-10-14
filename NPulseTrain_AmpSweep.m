% Written by Vinit Doke - October 2023
% N_PULSE_TRAIN : (ERS - N IDENTICAL PULSES OF GIVEN WIDTH AND AMPLITUDE - READ)

tic

clear all; close all; clc; % clear all variables, close all figures, clear command window
% loadlibrary('wgfmu', 'wgfmu.h'); 

% UNITS
ms = 1e-3;
us = 1e-6;
ns = 1e-9;

    % N SWEEP LINEAR : PARAMS
    N_MIN = 0;
    N_MAX = 200;
    STEPS = 10;

    EXPERIMENT_QUEUE = N_MIN : STEPS : N_MAX;
    
    % LOG SWEEP
    % N_MIN = 1; %10^1
    % N_MAX = 4; %10^4
    % STEPS = 1;
    % EXPERIMENT_QUEUE = N_MIN : STEPS : N_MAX;
    % EXPERIMENT_QUEUE = 10 .^ EXPERIMENT_QUEUE;

    % PULSE TRAIN Sweep
    V_PULSE_MIN  = -3.0; % Volts
    V_PULSE_MAX  = -2.0; % Volts
    V_PULSE_STEP = 0.5; % Volts

    V_PULSE_QUEUE = V_PULSE_MIN : V_PULSE_STEP : V_PULSE_MAX;

    t_ON_PULSE          = 1.0*us;
    t_OFF_PULSE         = 1.0*us;
    t_risefall_PULSE    = 0.01*us;
    t_delay_PULSE_TRAIN = 100*us;

    %% WAVEFORM VOLTAGE AND TIMING PARAMS:
    % ERS PULSE
    V_ERS          = 4.5; % Volts
    t_width_ERS    = 10.0*us;
    t_delay_ERS    = 1.0*us;
    t_risefall_ERS = 0.1*us;

    % VOLTAGE SWEEP READ PARAMS
    VG_READ_MIN   = -0.5; %V
    VG_READ_MAX   = 2.0; %V
    VD_READ       = 0.1; %V
    t_READ_VG     = 100*us;
    t_risefall_VG = 0.1*us;
    t_delay_READ  = 1.0*us;

    % DEVICE CONSTANTS
    GATE_LENGTH = 1000; % nm
    GATE_WIDTH = 1000; % nm
    It0 = 1e-7*(GATE_WIDTH/GATE_LENGTH); % For Calculating Threshold Voltage




% INSTRUMENT SETTINGS
I_RANGE_WRITE = 6002;
I_RANGE_READ  = 6003;

% Data Management
main_folder = fileManager('PulseTrain-AmplitudeSweep', ([N_MIN N_MAX STEPS V_PULSE_MIN V_PULSE_MAX V_PULSE_STEP t_ON_PULSE t_OFF_PULSE]));

% WRITE CONFIGURATION FILE
fid = fopen([main_folder 'config.txt'], 'w'); 
fprintf(fid, 'PARAMETER,VALUE\n');
fprintf(fid, 'N_MIN,%d\n', N_MIN);
fprintf(fid, 'N_MAX,%d\n', N_MAX);
fprintf(fid, 'STEPS,%d\n', STEPS);
fprintf(fid, 'V_PULSE_MIN,%d\n', V_PULSE_MIN);
fprintf(fid, 'V_PULSE_MAX,%d\n', V_PULSE_MAX);
fprintf(fid, 'V_PULSE_STEP,%d\n', V_PULSE_STEP);
fprintf(fid, 't_ON_PULSE,%d\n', t_ON_PULSE);
fprintf(fid, 't_OFF_PULSE,%d\n', t_OFF_PULSE);
fprintf(fid, 't_risefall_PULSE,%d\n', t_risefall_PULSE);
fprintf(fid, 't_delay_PULSE_TRAIN,%d\n', t_delay_PULSE_TRAIN);
fprintf(fid, 'VG_READ_MIN,%d\n', VG_READ_MIN);
fprintf(fid, 'VG_READ_MAX,%d\n', VG_READ_MAX);
fprintf(fid, 'VD_READ,%d\n', VD_READ);
fprintf(fid, 't_READ_VG,%d\n', t_READ_VG);
fprintf(fid, 't_risefall_VG,%d\n', t_risefall_VG);
fprintf(fid, 't_delay_READ,%d\n', t_delay_READ);
% time and date
% fprintf(fid, 'DATE,%s\n', datetime("now", 'Format', 'yyyy-MM-dd HH:mm:ss.SSS'));
fclose(fid);

% Data Containers
full_data_write = [];
full_data_read  = [];
VtD             = [];

data_out_write = [];
data_out_read = [];

% QUEUE ITERATION :
color = distinguishable_colors(length(EXPERIMENT_QUEUE));
voltage_amp_colors = distinguishable_colors(length(V_PULSE_QUEUE));

for j = 1:length(V_PULSE_QUEUE)

    % Put params in groups to pass to NPulse Time Scheme and Switch Generator
    ERS_PARAMS = {V_ERS, t_width_ERS, t_delay_ERS, t_risefall_ERS};
    PULSE_TRAIN_PARAMS = {V_PULSE_QUEUE(j), t_ON_PULSE, t_OFF_PULSE, t_risefall_PULSE, t_delay_PULSE_TRAIN};
    READ_PARAMS = {VG_READ_MIN, VG_READ_MAX, VD_READ, t_READ_VG, t_risefall_VG, t_delay_READ};


    % Data Containers
    VtD             = [];

    data_out_write = [];
    data_out_read = [];

    for i = 1:length(EXPERIMENT_QUEUE)
        N = EXPERIMENT_QUEUE(i);
        ins = NPulse(N, ERS_PARAMS, PULSE_TRAIN_PARAMS, READ_PARAMS);

        bw = ins{1}; br_vg = ins{2}; br_vd = ins{3};
        sample_interval_write = ins{4}; sample_points_write = ins{5};
        sample_interval_read = ins{6}; sample_points_read = ins{7};
        pt_packed = ins{8};

        % for self reference, not used by func write
        csvwrite([main_folder 'DATA/' 'voltage_pattern_write_' num2str(V_PULSE_QUEUE(j)) '_VPULSEAMP_' num2str(N) '_pulses.csv'], bw); 
        csvwrite([main_folder 'DATA/' 'voltage_pattern_read_gate.csv'], br_vg);
        csvwrite([main_folder 'DATA/' 'voltage_pattern_read_drain.csv'], br_vd);

        % ACTUAL FILE USED BY FUNC WRITE
        csvwrite('voltage_pattern_write.csv', bw);
        csvwrite('voltage_pattern_read_gate.csv', br_vg);
        csvwrite('voltage_pattern_read_drain.csv', br_vd);


        % send write instructions to device
        clear data_out_write;
        data_out_write = [];
        data_out_write(:,:) = func_write(sample_points_write, sample_interval_write, I_RANGE_WRITE);

        % send read instructions to device
        clear data_out_read;
        data_out_read = [];
        data_out_read(:,:) = func_WGFMU(sample_points_read, sample_interval_read,  I_RANGE_READ);


        % unpack data
        T      = data_out_write(:, 1);
        V_gate = data_out_write(:, 2);

        Tr        = data_out_read(:, 1);
        Vr_gate   = data_out_read(:, 2);
        Ir_drain  = abs(smooth(smooth(smooth(data_out_read(:, 4), 'sgolay', 4)))); %smoothing
        Ir_source = abs(smooth(smooth(smooth(data_out_read(:, 6), 'sgolay', 4)))); %smoothing 
        Ir_sub    = abs(smooth(smooth(smooth(data_out_read(:, 8), 'sgolay', 4)))); %smoothing

        % Calculate Threshold Voltage
            pshift = 10;
            Itd = double(Ir_source(pt_packed(3) + pshift : pt_packed(4) - pshift) > It0);
            vtg = Vr_gate(pt_packed(3) + pshift: pt_packed(4) - pshift);
        
            vtd_index = find(Itd, 1, 'first');
            Threshold = vtg(vtd_index);
            VtD = [VtD Threshold];

    
        PACK_DATA = {T, V_gate, Tr, Vr_gate, Ir_drain, Ir_source, Ir_sub, pt_packed, Threshold};

        % plot with data
        PlotData(PACK_DATA, N, color(i, :), voltage_amp_colors(j, :));

        % concat data
        full_data_write = [full_data_write; data_out_write];
        full_data_read  = [full_data_read; data_out_read];

        % Save Data for each N and V_PULSE_AMP
        WRITE_DATA = [T V_gate];
        READ_DATA = [Tr Vr_gate Ir_drain Ir_source Ir_sub];

        csvwrite([main_folder 'DATA/' 'WRITE_' num2str(V_PULSE_QUEUE(j)) '_VPULSEAMP_' num2str(EXPERIMENT_QUEUE(i)) '_PULSES.csv'], WRITE_DATA);
        csvwrite([main_folder 'DATA/' 'READOUT_' num2str(V_PULSE_QUEUE(j)) '_VPULSEAMP_' num2str(EXPERIMENT_QUEUE(i)) '_PULSES.csv'], READ_DATA);

    end

    ThresholdAndPulses = [EXPERIMENT_QUEUE' VtD'];

    % save data
    csvwrite([main_folder 'DATA/' num2str(V_PULSE_QUEUE(j)) '_VPULSEAMP' 'ThresholdAndPulses.csv'], ThresholdAndPulses);
    csvwrite([main_folder 'DATA/' num2str(V_PULSE_QUEUE(j)) '_VPULSEAMP_' 'full_data_write.csv'], full_data_write);
    csvwrite([main_folder 'DATA/' num2str(V_PULSE_QUEUE(j)) '_VPULSEAMP_' 'full_data_read.csv'], full_data_read);


    % save figures
    saveas(figure(1), [main_folder 'PLOTS/' num2str(V_PULSE_QUEUE(j)) '_VPULSEAMP_' 'GateVoltage_v_Time.png']);
    saveas(figure(2), [main_folder 'PLOTS/' num2str(V_PULSE_QUEUE(j)) '_VPULSEAMP_' 'DrainCurrent_v_Time.png']);
    saveas(figure(3), [main_folder 'PLOTS/' num2str(V_PULSE_QUEUE(j)) '_VPULSEAMP_' 'SourceCurrent_v_ThresholdVoltage.png']);
    saveas(figure(4), [main_folder 'PLOTS/' num2str(V_PULSE_QUEUE(j)) '_VPULSEAMP_' 'ThresholdVoltage_v_PulseNumber.png']);

    % save matlab figures
    hgsave(figure(1), [main_folder 'EXTRA/' num2str(V_PULSE_QUEUE(j)) '_VPULSEAMP_' 'GateVoltage_v_Time.fig']);
    hgsave(figure(2), [main_folder 'EXTRA/' num2str(V_PULSE_QUEUE(j)) '_VPULSEAMP_' 'DrainCurrent_v_Time.fig']);
    hgsave(figure(3), [main_folder 'EXTRA/' num2str(V_PULSE_QUEUE(j)) '_VPULSEAMP_' 'SourceCurrent_v_ThresholdVoltage.fig']);
    hgsave(figure(4), [main_folder 'EXTRA/' num2str(V_PULSE_QUEUE(j)) '_VPULSEAMP_' 'ThresholdVoltage_v_PulseNumber.fig']);

    % clear figures 1,2,3 for next V_PULSE_AMP
    clf(figure(1));
    clf(figure(2));
    clf(figure(3));

end



toc
