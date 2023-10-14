function data_out = func_WGFMU(sample_points, sample_interval,i_range)
%     unloadlibrary wgfmu
    % hfile = fullfile(matlabroot,'extern','include','wgfmu.h');
   
    calllib('wgfmu', 'WGFMU_clear');
    calllib('wgfmu', 'WGFMU_setTimeout',10);
    calllib('wgfmu', 'WGFMU_openLogFile', 'log.log');

    calllib('wgfmu', 'WGFMU_createPattern','ptn_1', 0); % 0 ms, 0 V
    calllib('wgfmu', 'WGFMU_createPattern','ptn_2', 0); % 0 ms, 0 V
    calllib('wgfmu', 'WGFMU_createPattern','ptn_3', 0); % 0 ms, 0 V
    calllib('wgfmu', 'WGFMU_createPattern','ptn_4', 0); % 0 ms, 0 V
    
    %% Channel 1 pattern
    M = csvread('voltage_pattern_read_gate.csv');
    for i = 1:length(M)
        calllib('wgfmu', 'WGFMU_addVector','ptn_1', M(i,1), M(i,2)); %0.1 ms, 1 sV
    end
    
    calllib('wgfmu', 'WGFMU_setMeasureEvent', 'ptn_1', 'event_1', 0, sample_points, sample_interval, sample_interval, 12000);
    calllib('wgfmu', 'WGFMU_addSequence', 102, 'ptn_1', 1);

    %% Channel 2 pattern
    
    N = csvread('voltage_pattern_read_drain.csv');
    for i = 1:length(N)
        calllib('wgfmu', 'WGFMU_addVector','ptn_2', N(i,1), N(i,2)); %0.1 ms, 1 sV
    end
    
    calllib('wgfmu', 'WGFMU_setMeasureEvent', 'ptn_2', 'event_1', 0, sample_points, sample_interval, sample_interval, 12000);
    calllib('wgfmu', 'WGFMU_addSequence', 201, 'ptn_2', 1);
%% channel 3
 
    calllib('wgfmu', 'WGFMU_createPattern','ptn_3', 0);
    calllib('wgfmu', 'WGFMU_addVector','ptn_3', sample_points*sample_interval, 0);
    calllib('wgfmu', 'WGFMU_setMeasureEvent','ptn_3', 'event_1', 0, sample_points, sample_interval,sample_interval, 12000);
    calllib('wgfmu', 'WGFMU_addSequence', 101, 'ptn_3', 1);
        
 %% channel 4 
 
    calllib('wgfmu', 'WGFMU_createPattern','ptn_4', 0);
    calllib('wgfmu', 'WGFMU_addVector','ptn_4', sample_points*sample_interval, 0);
    calllib('wgfmu', 'WGFMU_setMeasureEvent','ptn_4', 'event_1', 0, sample_points, sample_interval,sample_interval, 12000);
    calllib('wgfmu', 'WGFMU_addSequence', 202, 'ptn_4', 1);
        
      
        
    %%
%     
% 
%     calllib('wgfmu', 'WGFMU_createPattern','ptn_3', 0);
%     calllib('wgfmu', 'WGFMU_addVector','ptn_3', sample_points*sample_interval, 0);
%     calllib('wgfmu', 'WGFMU_setMeasureEvent','ptn_3', 'event_1',  0, sample_points, sample_interval,sample_interval, 12000);
    %%
    
    %%
% 	WGFMU_setMeasureEvent('ptn_1', "event_1", 0, 499, 8e-6, 8e-6, 12000);
%/*For higher resolution use 12001 points*(1 + int(average / (5e-9))), Table 4 - 14, Page4 - 75
%if the events change averaging conditions only (not range)
%Start time of adjacent events must be >= 100ns */

% For Reset operation current range may need to be change

% WGFMU_setRangeEvent('ptn_1', "event_2", 2e-4, 6004); // 6004 for +- 1 mA,

%/* 1. If the events change Ranges Start time of adjacent events must be >= 2us
%2. The range event must be set to a term out of average defined in the measurement event.*/

    
% 
% 
%     calllib('wgfmu', 'WGFMU_createPattern','ptn_2', 0);
%     calllib('wgfmu', 'WGFMU_addVector','ptn_2', sample_points*sample_interval, 0);
%     calllib('wgfmu', 'WGFMU_setMeasureEvent','ptn_2', 'event_1',  0, sample_points, sample_interval,sample_interval, 12000);
%     
    
%     calllib('wgfmu', 'WGFMU_addSequence',102, 'ptn_2', 1);

    calllib('wgfmu', 'WGFMU_openSession', 'GPIB0::17::INSTR');
    calllib('wgfmu', 'WGFMU_initialize');

    calllib('wgfmu', 'WGFMU_setOperationMode',101, 2001);  % 2001 (FastIV Mode)
    calllib('wgfmu', 'WGFMU_setOperationMode',102, 2001);
    calllib('wgfmu', 'WGFMU_setOperationMode',201, 2001);
    calllib('wgfmu', 'WGFMU_setOperationMode',202, 2001);

%     calllib('wgfmu', 'WGFMU_setForceVoltageRange',101, 3000);%3000 for Auto(default)
%     calllib('wgfmu', 'WGFMU_setForceVoltageRange',102, 3000);
% 
%     calllib('wgfmu', 'WGFMU_setMeasureEnabled',101, 7001);  % 7001 for Enable(Default)
%     calllib('wgfmu', 'WGFMU_setMeasureEnabled',102, 7001);

    calllib('wgfmu', 'WGFMU_setMeasureMode',102, 4000);  % 4000 for Voltage measurement mode
    calllib('wgfmu', 'WGFMU_setMeasureMode',101, 4001); % 4001 for current measurement mode
    calllib('wgfmu', 'WGFMU_setMeasureMode',201, 4001)
    calllib('wgfmu', 'WGFMU_setMeasureMode',202, 4001)
    
    
    calllib('wgfmu', 'WGFMU_setMeasureVoltageRange',102, 5002);  % 5001 for + / -5V
    calllib('wgfmu', 'WGFMU_setMeasureCurrentRange',101,i_range); % 6005 for 10mA (default)
    calllib('wgfmu', 'WGFMU_setMeasureCurrentRange',201,i_range);
    calllib('wgfmu', 'WGFMU_setMeasureCurrentRange',202,6001);
    
    calllib('wgfmu', 'WGFMU_connect',101);
    calllib('wgfmu', 'WGFMU_connect',102);
    calllib('wgfmu', 'WGFMU_connect',201);
    calllib('wgfmu', 'WGFMU_connect',202);

    calllib('wgfmu', 'WGFMU_execute');
    calllib('wgfmu', 'WGFMU_waitUntilCompleted');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% data collection
    p1 = libpointer('int32Ptr', zeros(1,1));
    calllib('wgfmu', 'WGFMU_getMeasureTimeSize', 101, p1);
    p2 = libpointer('int32Ptr', zeros(1,1));
    calllib('wgfmu', 'WGFMU_getMeasureTimeSize', 102, p2);
    p3 = libpointer('int32Ptr', zeros(1,1));
    calllib('wgfmu', 'WGFMU_getMeasureTimeSize', 201, p3);
    p4 = libpointer('int32Ptr', zeros(1,1));
    calllib('wgfmu', 'WGFMU_getMeasureTimeSize', 202, p4);
    

    T1 = libpointer('doublePtr', zeros(1, p1.value));
    V1 = libpointer('doublePtr', zeros(1, p1.value));
    T2 = libpointer('doublePtr', zeros(1, p2.value));
    I2  = libpointer('doublePtr', zeros(1, p2.value));
    T3 = libpointer('doublePtr', zeros(1, p3.value));
    I3  = libpointer('doublePtr', zeros(1, p3.value));  
    T4 = libpointer('doublePtr', zeros(1, p4.value));
    I4  = libpointer('doublePtr', zeros(1, p4.value)); 
    

    calllib('wgfmu', 'WGFMU_getMeasureValues',102, 0, p1.value, T1, V1);
    calllib('wgfmu', 'WGFMU_getMeasureValues',101, 0, p2.value, T2, I2);
    calllib('wgfmu', 'WGFMU_getMeasureValues',201, 0, p3.value, T3, I3);
    calllib('wgfmu', 'WGFMU_getMeasureValues',202, 0, p4.value, T4, I4);
    
%   data_out = [T1.value' I1.value' T2.value' V2.value' ];
    data_out = [T1.value' V1.value' T2.value' I2.value' T3.value' I3.value' T4.value' I4.value'];
%     csvwrite('data_out.csv', [T1.value' V1.value' T2.value' I2.value']);

    calllib('wgfmu', 'WGFMU_exportAscii', 'data_input.csv');

    calllib('wgfmu', 'WGFMU_closeLogFile');
    calllib('wgfmu', 'WGFMU_initialize');
%    calllib('wgfmu', 'WGFMU_disconnect', 101);
    calllib('wgfmu', 'WGFMU_closeSession');

end