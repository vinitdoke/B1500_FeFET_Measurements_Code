function instructions = NPulse(N, ERS_PARAMS, PULSE_TRAIN_PARAMS, READ_PARAMS)
    
    % FOR REFERENCE
    % ERS_PARAMS = {V_ERS, t_width_ERS, t_delay_ERS, t_risefall_ERS};
    % PULSE_TRAIN_PARAMS = {V_PULSE, t_ON_PULSE, t_OFF_PULSE, t_risefall_PULSE, t_delay_PULSE_TRAIN};
    % READ_PARAMS = {VG_READ_MIN, VG_READ_MAX, VD_READ, t_READ_VG, t_risefall_VG, t_delay_READ};

    ms = 1e-3;
    us = 1e-6;
    ns = 1e-9;


    % unpack
    V_ERS = ERS_PARAMS{1}; t_width_ERS = ERS_PARAMS{2}; t_delay_ERS = ERS_PARAMS{3}; t_risefall_ERS = ERS_PARAMS{4};
    V_PULSE = PULSE_TRAIN_PARAMS{1}; t_ON_PULSE = PULSE_TRAIN_PARAMS{2}; t_OFF_PULSE = PULSE_TRAIN_PARAMS{3}; t_risefall_PULSE = PULSE_TRAIN_PARAMS{4}; t_delay_PULSE_TRAIN = PULSE_TRAIN_PARAMS{5};
    VG_READ_MIN = READ_PARAMS{1}; VG_READ_MAX = READ_PARAMS{2}; VD_READ = READ_PARAMS{3}; t_READ_VG = READ_PARAMS{4}; t_risefall_VG = READ_PARAMS{5}; t_delay_READ = READ_PARAMS{6};

    
    % ERS
    tScheme_ERS = [0 t_delay_ERS t_risefall_ERS t_width_ERS t_risefall_ERS t_delay_ERS];
    vScheme_ERS = [0 0 V_ERS V_ERS 0 0];

    % PULSE TRAIN
    tScheme_SINGLE_PULSE = [0 t_risefall_PULSE t_ON_PULSE t_risefall_PULSE t_OFF_PULSE];
    vScheme_SINGLE_PULSE = [0 V_PULSE V_PULSE 0 0];

    tScheme_PULSE_TRAIN = [0 t_delay_PULSE_TRAIN];
    vScheme_PULSE_TRAIN = [0 0];

    %TODO: CONSIDER PREALLOCATION HERE
    for i = 1:N
        tScheme_PULSE_TRAIN = [tScheme_PULSE_TRAIN tScheme_SINGLE_PULSE];
        vScheme_PULSE_TRAIN = [vScheme_PULSE_TRAIN vScheme_SINGLE_PULSE];
    end

    % add delay at end of pulse train
    tScheme_PULSE_TRAIN = [tScheme_PULSE_TRAIN t_delay_PULSE_TRAIN];
    vScheme_PULSE_TRAIN = [vScheme_PULSE_TRAIN 0];

    % READ
    tScheme_READ = [0 t_delay_READ t_risefall_VG t_READ_VG t_risefall_VG t_delay_READ];
    vScheme_READ_VG = [0 0 VG_READ_MIN VG_READ_MAX 0 0];
    vScheme_READ_VD = [0 0 VD_READ VD_READ 0 0];

    % SAMPLE RATES
    % ERS + PULSE TRAIN
        % combine ERS and PULSE TRAIN
        tScheme_ERS_PULSE_TRAIN = [tScheme_ERS tScheme_PULSE_TRAIN];
        vScheme_ERS_PULSE_TRAIN = [vScheme_ERS vScheme_PULSE_TRAIN];

    sample_interval_write = 10*ns;
    t_final_WRITE = cumsum(tScheme_ERS_PULSE_TRAIN);
    sample_points_write = round(t_final_WRITE(end)/sample_interval_write);

    % READ
    sample_points_read = 1000; % bring this as a param to top
    t_final = cumsum(tScheme_READ);
    t_final_read = t_final(end);
    sample_interval_read = t_final_read/sample_points_read;

    if sample_interval_read < 10*ns
        sample_interval_read = 10*ns;
        sample_points_read = round(t_final_read/sample_interval_read);
    end

    pt0 = round(sample_points_read*t_final(1)/t_final(end));
    pt1 = round(sample_points_read*t_final(2)/t_final(end));
    pt2 = round(sample_points_read*t_final(3)/t_final(end));
    pt3 = round(sample_points_read*t_final(4)/t_final(end));
    pt4 = round(sample_points_read*t_final(5)/t_final(end));
    pt5 = round(sample_points_read*t_final(6)/t_final(end));

    pt_packed = [pt0 pt1 pt2 pt3 pt4 pt5];

    % GENERATE SWITCHING INSTRUCTIONS
    % ERS + PULSE TRAIN
    B_WRITE = [tScheme_ERS_PULSE_TRAIN' vScheme_ERS_PULSE_TRAIN'];

    B_READ_VG = [tScheme_READ' vScheme_READ_VG'];
    B_READ_VD = [tScheme_READ' vScheme_READ_VD'];

    instructions = {B_WRITE, B_READ_VG, B_READ_VD, sample_interval_write, sample_points_write, sample_interval_read, sample_points_read, pt_packed};
end

