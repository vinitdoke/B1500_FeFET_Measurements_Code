function PlotData(DATA, N, pulse_number_color, voltage_amp_color)

    plot_width = 750; plot_height = 500;

    % unpack data
    T = DATA{1}; V_gate = DATA{2};
    Tr = DATA{3}; Vr_gate = DATA{4}; Ir_drain = DATA{5}; Ir_source = DATA{6};
    Threshold = DATA{9};


    % Gate Voltage v/s Time
    figure(1) 
    set(gcf, 'position', [0, 0, 750, 500]) % (x, y, width, height)
    plot(T, V_gate,  'linewidth', 2, 'color', pulse_number_color);
    box on
    grid on
    title(['Gate Voltage v/s Time; Latest Pulse# = ' num2str(N)]);
    xlabel('Time (s)');
    ylabel('Gate Voltage (V)');
    hold all;
    drawnow();


    % Drain Current v/s Time
    figure(2) 
    set(gcf, 'position', [plot_width, 0, plot_width, plot_height]); % (x, y, width, height)
    semilogy(Tr, Ir_drain, 'linewidth', 2, 'color', pulse_number_color) %
    box on
    grid on
    xlim([Tr(1) Tr(end)]);
    title('Drain Current v/s Time');
    xlabel('Time (s)');
    ylabel('Drain Current (A)');

    hold all;
    drawnow();


    % ID v/s VG
    figure(3) 
    set(gcf, 'position', [0, plot_height, plot_width, plot_height]);
    semilogy(Vr_gate, Ir_source, 'color', pulse_number_color, 'LineWidth', 2);
    box on
    grid on
    title('Source Current v/s Gate Voltage');
    xlabel('VG (V)');
    ylabel('ID (A)');

    hold all;
    drawnow();


    % Threshold Voltage
    figure(4) % Threshold Voltage v/s Number of Pulses
    set(gcf, 'position', [plot_width, plot_height, plot_width, plot_height]);
    plot(N, Threshold, 'o', 'color', voltage_amp_color, 'LineWidth', 2);
    box on
    grid on
    title('Threshold Voltage v/s Number of Pulses');
    xlabel('Pulse Number #');
    ylabel('Threshold Voltage (V)');

    hold all;
    drawnow();

end