% endurance code without READ i.e. only SET and RESET transients
% RESET first and then SET
tic
clear all;
clc
close all;
 loadlibrary('wgfmu','wgfmu.h');
ns = 1e-9;%s
us = 1e-6;%s
t_start = 0*ns;
ms = 1e-3;
% sample_interval = ;


L=500; %Gate Length in nm
W=80; %Gate Widthth in nm
v_drain = [];
i_range_write = 6002; % 6001;
i_range_read = 6003;% edit now
rep_count = 1; % two Positive pulses and two negative pulses

% t_delay =10*us;     % longer delay is needed
t_rise = 1*us;      % check 10ns may give too much capacitive current
t_fall = t_rise;      % same as above
pw = 1*us;
delay =10*us;
% sample_points_write = 5000;

It0 = 1e-6*(W/L);
% It0 = 5e-7;
%   v_sw = [4 -0.5 -1 -1.5 -2 -2.5 -3 -3.5 -4 -4 -4 0.5 1 1.5 2 2.5 3 3.5 4];
%  z=0.25:0.25:4;
%  v_sw=[-4 -4 z];
%  v_sw=-4:0.5:4
%  v_sw=4:-1:-4;
% v_sw=0:-1:-4;
% v_sw=0:1:4;
%  v_sw = [4 -4 4 -3.5 4 -3 4 -2.5 4 -2 4 -1.5 4 -1 4 -0.5 -4 0.5 -4 1 -4 1.5 -4 2 -4 2.5 -4 3 -4 3.5 -4 4];
% v_sw=4*(2*rem(1:6,2)-1);
%  v_sw = [4.5 -4.5 4.5 -4.5 4.5 -4.5 4.5 -4.5 4.5 -4.5];
%   v_sw = [4.5];
%   v_sw = [4.5 -4.5 4.5 -4.5 4.5 -4.5 4.5 -4.5 4.5 -4.5 4.5 -4.5 4.5 -4.5 4.5 -4.5 4.5 -4.5 4.5 -4.5];
% vsw1=zeros(1,49);
% 
%  v_sw = [4.5 vsw1];

v_sw = [4.5 -4.5 4.5 -4.5 4.5 -4.5 4.5 -4.5 4.5 -4.5];
% v_sw = [-4.5 4.5 -4.5 4.5 -4.5 4.5 -4.5 4.5 -4.5 4.5 -4.5 4.5];
% v_sw = [-5 4.5 -5 4.5 -5 4.5 -5 4.5 -5 4.5 -5 4.5];

%% write voltage and time scheme
t_sw = [0 delay t_rise/10 pw t_fall/10 delay];
for h = 1:length(v_sw)
    v_switch(h,:) = ([0 0  v_sw(h) v_sw(h) 0 0 ]);
end

tw = [];
vw = [];
x = 0;

for k = 1:rep_count
        tw = [tw t_sw];
        vw = [vw v_switch];
end 

vw =vw'; % source voltage is negative of Vg
tw = tw';
t_finalw = cumsum(tw);

% sample_interval_write = t_finalw(end)/sample_points_write;

% if sample_interval_write <10*ns
%     sample_interval_write=50*ns;
        sample_interval_write=1*us;
%     sample_points_write=t_finalw(end)/sample_interval_write;
% end
sample_points_write = round(t_finalw(end)/sample_interval_write);

%% read voltage and time scheme
% gate
 vg_read1 =-0.5*ones(1,length(v_sw));
 vg_read2 =2.0*ones(1,length(v_sw));
% vg_read2 =[0.5];
% vg_read1 =[-0.5 0 -0.5 0 -0.5 0.4 -0.5 0 -0.5 0];

% % % vg_read1 =[-0.5 0.5 -0.5 0.5 -0.5 0.5 -0.5 0.5 -0.5 0.5];
% vg_read2 =[0.5 1 0.5 1 0.5 1];
% vg_read2 =[1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5 1.5];
% vg_read2 =[1 1.5 1 1.5 1 1.5 1 1.5 1 1.5];
% vg_read1 =[-0.5 0.8 -0.5 0.8 -0.5 0.8 -0.5 0.8 -0.5 0.8 -0.5 0.8 -0.5 0.8 -0.5 0.8 -0.5 1 -0.5 0.8];
% vg_read2 =[0.5 1.5 0.5 1.5 0.5 1.5 0.5 1.5 0.5 1.5 0.5 1.5 0.5 2 0.5 1.5 0.5 1.5 0.5 1.5];
% vg_read1=-0.5*ones(1,50);
% vg_read2=0.5*ones(1,50);
% vg_read1 =[-0.5 0.8 -0.5 0.8 -0.5 0.8];
% vg_read2 =[ 0.5 2 0.5 2 0.5 2 ];
% vg_read1 =[-0.5 0.8] ;
% vg_read2 =[ 0.4 2.2];
% vg_read2 = 1*ones(1,length(v_sw));
%  vg_read1=-0.5*ones(1,length(v_sw));
% vg_read2=2*ones(1,6);
tfac = 1;
t_read1 = 100*us; % changed to 1000us to 100us on 24-08-21
t_rise1 = 1*us;
t_fall1 = t_rise1;
delay1 = 10*us;
no_of_points = round((t_read1+2*t_rise+2*delay1)/us);
% sample_points_read = no_of_points*2;
sample_points_read = 1000;  % changed to 1000 from 10000 on 24-08-21

tg_read = [0 delay1 t_rise1 t_read1 t_fall1 delay1]/tfac;
for y=1:length(vg_read1)
    vg_read(y,:)=([0 0 vg_read1(y) vg_read2(y) 0 0]);
end
% vg_read = ([0 0  vg_read1 vg_read2 0 0 ]);

% drain
vd_read1=0.1;
% vd_read1 = 1;
vd_read = ([0 0  vd_read1 vd_read1 0 0 ]);

tg = [];
vd= [];
vg = [];
x = 0;
Vts = [];

        tg =  tg_read;
        vd = vd_read;
        vg = vg_read;
     

vd =vd'; % source voltage is negative of Vg
vg = vg';
tg = tg';
t_final = cumsum(tg);

% figure(11)
% plot(cumsum(tg),vg)


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

%
%     Br = [tg vg];
%     csvwrite('voltage_pattern_read_gate.csv',Br);
%     Ar=[tg vd];
%     csvwrite('voltage_pattern_read_drain.csv',Ar);
    % creating input pattern     
dvdt = [];

%  plot(t_final,vg)
m = 0;
n = 0;
%%
OC = length(v_sw);
full_data =  [];
cycles =1;
color  = distinguishable_colors(OC);
C = [];
filename = input(' filename:    ','s');   %Enter the file NAME
% y = 0;
full_data_write = [];
full_data_read = [];
% writing
Vts = [];
VtD = [];
for i = 1:OC
    i

    Bw = [tw vw(:,i)];
    csvwrite('voltage_pattern_write.csv',Bw);
    
   data_out_write(:,:,i) = func_write(sample_points_write, sample_interval_write,i_range_write);


 V_gate(:,i) = data_out_write(:,2,i);
 I_drain(:,i) = abs(data_out_write(:,4,i));
 I_source(:,i) =abs(data_out_write(:,6,i));
 I_sub(:,i) =abs( data_out_write(:,8,i));
 T(:,i) = data_out_write(:,1,i);

   full_data_write = [full_data_write data_out_write(:,:,i)];
    
%% read
if rem(i,2)==1
    Br = [tg vg(:,i)];
    csvwrite('voltage_pattern_read_gate.csv',Br);
    Ar=[tg vd];
    csvwrite('voltage_pattern_read_drain.csv',Ar);
else
     Br = [tg vg(:,i)];
    csvwrite('voltage_pattern_read_gate.csv',Br);
    Ar=[tg vd];
    csvwrite('voltage_pattern_read_drain.csv',Ar);
end



   data_out_read(:,:,i) = func_WGFMU(sample_points_read, sample_interval_read,i_range_read);
%    data_out(:,4,i) = smooth(data_out(:,4,i),.1,'rloess');



 Vr_gate(:,i) = data_out_read(:,2,i);
 Ir_drain(:,i) = abs(smooth(smooth(smooth(data_out_read(:,4,i),'sgolay',4))));
 Ir_source(:,i) =abs(smooth(smooth(smooth( data_out_read(:,6,i),'sgolay',4))));
 Ir_sub(:,i) = abs(smooth(smooth( data_out_read(:,8,i),'sgolay',4)));
 Tr(:,i) = data_out_read(:,1,i);

   full_data_read = [full_data_read data_out_read(:,:,i)];
    

%     %% plotting
       figure (7)
    fig1 = plot(T(:,i), V_gate(:,i), 'color', color(i,:), 'linewidth',2);
%     fig1 = plot(T(1000:end,i), V_gate(1000:end,i), 'color', color(i,:), 'linewidth',2);
    title([filename 'Write voltage'],'FontSize',16,'FontWeight','b');
    xlabel('Time (s)','FontSize',16,'FontWeight','b');
    ylabel ('Voltage (V)','FontSize',16,'FontWeight','b');
  xlim([T(1,i) T(end,i)]);
%     ylim([0 3.5]);
 h = gca; 
    set(h,'FontSize',16,'FontWeight','b'); 
     grid on 
    hold all; 
    h2 = gca;
    h3 = gcf;
%     set(h2,'FontSize',12,'FontWeight','b');
    set(h3,'Color',[1,1,1]);
     grid on
     box on   
    h = gca;
    set(h,'FontSize',16,'FontWeight','b');
     grid on
    hold all;
    
  
%     fig1_name = [filename,'_ResetVol'];
% fig1_name = [filename ' V-t'] ;
%     hgsave(fig1_name);
%     saveas(fig1,[fig1_name '.jpg'])
%     save([filename '.mat']);
  
    %% I source-t
%     
%     figure (2)
%     fig2 = semilogy(T(:,i), I_source(:,i), 'color', color(i,:), 'linewidth',2);hold all;
% %     fig2 = semilogy(T(1000:end,i), I_drain(1000:end,i), 'color', color(i,:), 'linewidth',2);hold all;
%     title([filename 'Write source current'],'FontSize',16,'FontWeight','b');
%     xlabel('Time (s)','FontSize',16,'FontWeight','b');
%     ylabel ('Source Current (A)','FontSize',16,'FontWeight','b');
%    xlim([T(1,i) T(end,i)]);
%     h = gca;
%     set(h,'FontSize',16,'FontWeight','b');
%     hold all;
%      grid on
% %     fig2_name = [filename,'_ResetCur'];
% fig2_name = [filename ' I-t'];
%     hgsave(fig2_name);
%     saveas(fig2,[fig2_name '.jpg'])
%     
    
%% V_gate-t
%     figure (3)
%     fig3 = plot(Tr(:,i), Vr_gate(:,i), 'color', color(i,:), 'linewidth',2);
% %     fig1 = plot(T(1000:end,i), V_gate(1000:end,i), 'color', color(i,:), 'linewidth',2);
%     title([filename 'Gate read voltage'],'FontSize',16,'FontWeight','b');
%     xlabel('Time (s)','FontSize',16,'FontWeight','b');
%     ylabel ('Voltage (V)','FontSize',16,'FontWeight','b');
%     xlim([Tr(1,i) Tr(end,i)]);
% %     ylim([0 3.5]);
%     h = gca;
%     set(h,'FontSize',16,'FontWeight','b');
%      grid on
%     hold all;
% %     fig1_name = [filename,'_ResetVol'];
% fig3_name = [filename ' Vgate_read-t'] ;
%     hgsave(fig3_name);
%     saveas(fig3,[fig3_name '.jpg'])
% %     save([filename '.mat']);
%   
    %% V_drain-t
    
     
    figure (1)
    fig1 = semilogy(Tr(:,i), Ir_sub(:,i), 'color', color(i,:), 'linewidth',2);hold all;
%     fig2 = semilogy(T(1000:end,i), I_drain(1000:end,i), 'color', color(i,:), 'linewidth',2);hold all;
    title([filename 'Substrate current'],'FontSize',16,'FontWeight','b');
    xlabel('Time (s)','FontSize',16,'FontWeight','b');
    ylabel (' Substrate Current (A)','FontSize',16,'FontWeight','b');
    xlim([Tr(1,i) Tr(end,i)]);
    h = gca;
    set(h,'FontSize',16,'FontWeight','b');
    hold all;
     grid on
      h = gca; 
    set(h,'FontSize',16,'FontWeight','b'); 
     grid on 
    hold all; 
    h2 = gca;
    h3 = gcf;
%     set(h2,'FontSize',12,'FontWeight','b');
    set(h3,'Color',[1,1,1]);
     grid on
     box on   
%     fig2_name = [filename,'_ResetCur'];
fig1_name = [filename ' Idrain_read-t'];
    hgsave(fig1_name);
    saveas(fig1,[fig1_name '.jpg'])
     
    figure (2)
    fig2 = semilogy(Tr(:,i), Ir_drain(:,i), 'color', color(i,:), 'linewidth',2);hold all;
%     fig2 = semilogy(T(1000:end,i), I_drain(1000:end,i), 'color', color(i,:), 'linewidth',2);hold all;
    title([filename ' Drain current'],'FontSize',16,'FontWeight','b');
    xlabel('Time (s)','FontSize',16,'FontWeight','b');
    ylabel ('Drain Current (A)','FontSize',16,'FontWeight','b');
%     xlim([Tr(1,i) Tr(end,i)]);
    h = gca;
    set(h,'FontSize',16,'FontWeight','b');
    hold all;
     grid on
      h = gca; 
    set(h,'FontSize',16,'FontWeight','b'); 
     grid on 
    hold all; 
    h2 = gca;
    h3 = gcf;
%     set(h2,'FontSize',12,'FontWeight','b');
    set(h3,'Color',[1,1,1]);
     grid on
     box on   
%     fig2_name = [filename,'_ResetCur'];
fig2_name = [filename ' Idrain_read-t'];
    hgsave(fig2_name);
    saveas(fig2,[fig2_name '.jpg'])
    
    figure (4)
    fig4 = semilogy(Tr(:,i), Ir_source(:,i), 'color', color(i,:), 'linewidth',2);hold all;
%     fig2 = semilogy(T(1000:end,i), I_drain(1000:end,i), 'color', color(i,:), 'linewidth',2);hold all;
    title([filename ' Source current'],'FontSize',16,'FontWeight','b');
    xlabel('Time (s)','FontSize',16,'FontWeight','b');
    ylabel ('Source Current (A)','FontSize',16,'FontWeight','b');
%     xlim([Tr(1,i) Tr(end,i)]);
    h = gca;
    set(h,'FontSize',16,'FontWeight','b');
    hold all;
     grid on
      h = gca; 
    set(h,'FontSize',16,'FontWeight','b'); 
     grid on 
    hold all; 
    h2 = gca;
    h3 = gcf;
%     set(h2,'FontSize',12,'FontWeight','b');
    set(h3,'Color',[1,1,1]);
     grid on
     box on   
%     fig2_name = [filename,'_ResetCur'];
fig4_name = [filename ' Irsouce_read-t'];
    hgsave(fig4_name);
    saveas(fig4,[fig4_name '.jpg'])
%    save([filename '.mat']);

  pshift = 10;

     figure (5)
    fig5 = semilogy(Vr_gate(pt2+pshift:pt3-pshift,i), Ir_source(pt2+pshift:pt3-pshift,i), 'color', color(i,:), 'linewidth',2);hold all;
%     fig2 = semilogy(T(1000:end,i), I_drain(1000:end,i), 'color', color(i,:), 'linewidth',2);hold all;
    title([filename 'Id-Vg'],'FontSize',16,'FontWeight','b');
    xlabel('Voltage(V)','FontSize',16,'FontWeight','b');
    ylabel ('Source current(A)','FontSize',16,'FontWeight','b');
%     xlim([Vr_gate(41,i) Vr_gate(19960,i)]);
    h = gca;
    set(h,'FontSize',16,'FontWeight','b');
    hold all;
     grid on
      h = gca; 
    set(h,'FontSize',16,'FontWeight','b'); 
     grid on 
    hold all; 
    h2 = gca;
    h3 = gcf;
%     set(h2,'FontSize',12,'FontWeight','b');
    set(h3,'Color',[1,1,1]);
     grid on
     box on   
%     fig2_name = [filename,'_ResetCur'];
fig5_name = [filename ' Vgate Vs I source-t'];
    hgsave(fig5_name);
    saveas(fig5 ,[fig5_name '.jpg'])
    
  
       %% Vt Extraction: Id=1e-6*W/L
if  mod(i,2) == 1
    m = m+1;
  Itd(:,i)= double(Ir_source(pt2+pshift:pt3-pshift,i)>It0);
  vtg = Vr_gate(pt2+pshift:pt3-pshift,i);
%   if sum(Itd(:,i)) == 0
% %     VtD = [VtD VtD(m-1)];
%     V_SW_D(m) = v_sw(i);    
%   end 
  vtd_index(m) = find(Itd(:,i),1,'first');
  VtD = [VtD vtg(vtd_index(m))];
  V_SW_S(m) = v_sw(i);
else
     m = m+1;
  Itd(:,i)= double(Ir_source(pt2+pshift:pt3-pshift,i)>It0);
  vtg = Vr_gate(pt2+pshift:pt3-pshift,i);
%   if sum(Itd(:,i)) == 0
% %     VtD = [VtD VtD(m-1)];
%     V_SW_D(m) = v_sw(i);    
%   end 
  vtd_index(m) = find(Itd(:,i),1,'first');
 VtD = [VtD vtg(vtd_index(m))];
  V_SW_S(m) = v_sw(i);
end
%   
% if  mod(i,2) == 0 
%     n = n+1;
%   Its(:,i)= double(Ir_source(:,i)>1e-6*W/L);
%   if sum(Its(:,i)) == 0
%     Vts = [Vts Vts(m-1)];
%     V_SW_S(m) = v_sw(i);    
%   end 
%   vts_index(n) = find(Its(:,i),1,'first');
%   Vts = [Vts Vr_gate(vts_index(n),i)];
%   V_SW_S(n) = v_sw(i);
% end 
    
     figure (6)
    fig6 = plot(V_SW_S(m),VtD(m),'-o', 'color', color(i,:), 'linewidth',2);hold all;
    
%     fig6 = plot(V_SW_S,Vts,'-rs', 'linewidth',2);hold all;
%     fig2 = semilogy(T(1000:end,i), I_drain(1000:end,i), 'color', color(i,:), 'linewidth',2);hold all;
    title([filename 'Vwrite vs Vt'],'FontSize',16,'FontWeight','b');
    xlabel('Switch Voltage(V)','FontSize',16,'FontWeight','b');
    ylabel ('Vt (V)','FontSize',16,'FontWeight','b');
%     xlim([0 data_out_write(end,3,i)]);
    h = gca;
    set(h,'FontSize',16,'FontWeight','b');
    hold all;
     grid on
      h = gca; 
    set(h,'FontSize',16,'FontWeight','b'); 
     grid on 
    hold all; 
    h2 = gca;
    h3 = gcf;
%     set(h2,'FontSize',12,'FontWeight','b');
    set(h3,'Color',[1,1,1]);
     grid on
     box on   
%     fig2_name = [filename,'_ResetCur'];
fig6_name = [filename ' Vsw Vs Vt'];
    hgsave(fig6_name);
    saveas(fig6,[fig6_name '.jpg'])
    
    
end
% for kk=1:length(VtD)
%     if rem(kk,2)==1
%         VtLow(kk)=VtD(kk);
%     elseif (rem(kk,2)==0)
%         VtHigh(kk)=VtD(kk);
%     end
% end
% Vt_Low=mean(VtLow);
% Vt_High=mean(VtHigh);

 save([filename '.mat']);
   
csvwrite([filename 'write.csv'],[full_data_write] );
csvwrite([filename 'read.csv'],[full_data_read]);
VtD
 toc
