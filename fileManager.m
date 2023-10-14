function path = fileManager(experiment_type, params)

    DEVICE = input('Enter DUT ID:  ', 's');
    RUN_NUMBER = input('Enter run number:  ');
    EXP_TYPE = experiment_type;

    disp(EXP_TYPE)

    if strcmp(experiment_type, 'WakeUp')
        EXP_PARAMS = input('Enter params:  ', 's');
    elseif strcmp(experiment_type, 'PulseTrain')
        % [N_MIN N_MAX STEPS V_PULSE t_ON_PULSE]
        EXP_PARAMS = ['NMIN-' num2str(params(1)) '_NMAX-' num2str(params(2)) '_STEPS-' num2str(params(3)) '_VPULSE-' num2str(params(4)) '_tONPULSE-' num2str(params(5))];
    elseif strcmp(experiment_type, 'PulseTrain-AmplitudeSweep')
        % [N_MIN N_MAX STEPS V_PULSE_MIN V_PULSE_MAX V_PULSE_STEP t_ON_PULSE t_OFF_PULSE]
        EXP_PARAMS = ['NMIN-' num2str(params(1)) '_NMAX-' num2str(params(2)) '_STEPS-' num2str(params(3)) '_VPULSEMIN-' num2str(params(4)) '_VPULSEMAX-' num2str(params(5)) '_STEPS-' num2str(params(6)) '_tONPULSE-' num2str(params(7)) '_tOFFPULSE-' num2str(params(8))];    
    end

    FOLDER_NAME = strcat(DEVICE, '_', EXP_TYPE, '_', EXP_PARAMS, '_', 'RUN_',num2str(RUN_NUMBER));

    % check if folder exists
    if exist(FOLDER_NAME, 'dir') == 7
        disp('Folder already exists');
        return;
    end

    mkdir(FOLDER_NAME);

    %subfolders - DATA, PLOTS, EXTRA

    DATA_FOLDER = strcat(FOLDER_NAME, '/DATA');
    PLOTS_FOLDER = strcat(FOLDER_NAME, '/PLOTS');
    EXTRA_FOLDER = strcat(FOLDER_NAME, '/EXTRA');

    mkdir(DATA_FOLDER);
    mkdir(PLOTS_FOLDER);
    mkdir(EXTRA_FOLDER);

    path = strcat(FOLDER_NAME, '/');

end

    


