function varargout = LNVmodel(ZT, kv2choice)

%
%LNVmodel (Written by P Smith, University of Bristol)
%Models activity of Drosophila LNV neuron
%
%*Call function --> LNVmodel;
%
%*Assign output --> [vrec,I1,I2,I3,I4] = LNVmodel;
%Returns membrane potential (vrec)
%        Kv1 current (I1)
%        Kv2 current (I2)
%        Kv3 current (I3)
%        Kv4 current (I4)
%
%Enter time of day as a response in the form of ZTX (zeitgeber time)
%Where ZT0 is sunrise/lights-on and ZT12 is sunset/lights-off
%ZT24 is the same as ZT0
%
%Enter form of Kv2 model as either 1, 2, or 3
%Where 1 is the normal native Kv2 current
%      2 is the native Kv2 with human wild-type Kv9
%      3 is the native Kv2 with human mutant Kv9 (c379E)


%frequency = zeros(1,24);
%k41a = zeros(24,1000001);
%k42a = zeros(24,1000001);

%for iters=1:24
%% Parameters

C=3.7; 

%call = 'What is the time (ZT)? [Input: 0-24]:';                             %User input for time of day
%ZT = input(call);

% Integration settings

DT=0.01;
t0=0;
tf=12000;                                                                   %Time is in ms
t=t0:DT:tf;
numT=length(t);

load('MODEL.mat');                                                          %Loads channel parameters
Pars1 = KV1;
Pars3 = KV32;
Pars4 = KV4;

Ena=52; Ek=-90; El=-7; Eca=132;                                             %Reversal potentials based on Nernst equation
   
g_leak = ((0.0043*cos(0.2618*ZT))+0.0906);                                  %Use user input to calculate the conductances
kv3 = ((0.07075*cos(0.2618*ZT))+0.09575);                                   %for Kv3, Kv4, and the leak conductance
kv4 = ((0.225*-cos(0.2618*ZT))+0.275);
g_k3 = kv3*Pars3(11); g_k4 = kv4*Pars4(11);
g_ca=3.86;  g_k1 = Pars1(11);                                               %Major fixed conductances
sigma=0.05; iapp = zeros(numT,1);                                           %Noise setting

%Choose a Kv2 current
%kv2call = 'What form of Shab? (1=Native, 2=With hKv9WT, 3=With hKv9MU):';
%kv2choice = input(kv2call);
if     kv2choice == 1;                                                      %Normal Kv2 current
    Pars2 = KV2; g_na=60; iapp(:) = -4.97;
elseif kv2choice == 2;                                                      %Kv2 with human wild-type Kv9
    Pars2 = KV9W; g_na=62; iapp(:) = -4.8;    
elseif kv2choice == 3;                                                      %Kv2 with mutant c379E Kv9
    Pars2 = KV9M; g_na=58.9; iapp(:) = -4.32; g_ca=g_ca*0.73; g_k3=g_k3*1.1;  
else
    disp('*Inapporpriate answer.')
    disp('*Enter Kv2 type as either: 1, 2, or 3.')
    disp('*Please try again.')
    return
end

    g_k2 = Pars2(11);

% Initial conditions

v=-55; 
mna=0.22; hna=0.02; mca=0.09; hca=0.032; 
mk1=0.13; hk1=0.31; mk2 = 0.2; hk2 = 0.46; 
mk4 = 0.45; hk4 = 0.85; mk3 = 0.52; hk3 = 0.985;

vrec=zeros(numT,1);                                                         %Handles used to retrieve variables for the entire run
vrec(1)=v;
vdiff=zeros(numT,1);
na1=zeros(numT,1);
na2=zeros(numT,1);
ca1=zeros(numT,1);
ca2=zeros(numT,1);
k11=zeros(numT,1);
k12=zeros(numT,1);
k21=zeros(numT,1);
k22=zeros(numT,1);
k31=zeros(numT,1);
k32=zeros(numT,1);
k41=zeros(numT,1);
k42=zeros(numT,1);
leak=zeros(numT,1);
I1=zeros(numT,1);
I2=zeros(numT,1);
I3=zeros(numT,1);
I4=zeros(numT,1);
%% Integrate using Euler-Maryuma method

for ix=2:numT    
    %% gating functions
         
    % sodium activation
    mna_inf = 1/(1+exp(-(v+35.2)/7.9));                                     % In the form; 1./(1+exp(-(V-vh)/k))
    tau_mna = exp(-(v+286)/160);                                            %              a*exp(-((V-b)/c))
    
    % sodium inactivation 
    hna_inf = 1/(1+exp((v+62)/5.5));                                        % In the form; 1./(1+exp((V-vh)/k));
    tau_hna = 0.51+exp(-(v+26.6)/7.1);                                      %              a*exp(-((V-b)/c));
  
    % calcium activation 
    mca_inf = 1/(1+exp(-(v+25)/7.5));
    tau_mca = 3.1;

    % calcium inactivation 
    hca_inf = 1/(1+exp((v+260)/65));
    tau_hca = exp(-(v-444)/220);
  
    % potassium1 activation
    mk1_inf = 1/((1+exp(-(v-Pars1(1))/Pars1(2))));
    tau_mk1 = Pars1(3)*exp(-((v-Pars1(4))/Pars1(5)));

    % potassium1 inactivation
    hk1_inf = 1/((1+exp((v-Pars1(6))/Pars1(7))));
    tau_hk1 = Pars1(8)*exp(-((v-Pars1(9))/Pars1(10)));
   
    % potassium2 activation
    mk2_inf = 1/((1+exp(-(v-Pars2(1))/Pars2(2))));
    tau_mk2 = Pars2(3)*exp(-((v-Pars2(4))/Pars2(5)));

    % potassium2 inactivation
    hk2_inf = 1/((1+exp((v-Pars2(6))/Pars2(7))));
    tau_hk2 = Pars2(8)*exp(-((v-Pars2(9))/Pars2(10)));
    
    % potassium3 activation
    mk3_inf = 1/((1+exp(-(v-Pars3(1))/Pars3(2))));
    tau_mk3 = Pars3(3)*exp(-((v-Pars3(4))/Pars3(5)));

    % potassium3 inactivation
    hk3_inf = 1/((1+exp((v-Pars3(6))/Pars3(7))));
    tau_hk3 = Pars3(8)*exp(-((v-Pars3(9))/Pars3(10))); 
    
    % potassium4 activation
    mk4_inf = 1/((1+exp(-(v-Pars4(1))/Pars4(2))));
    tau_mk4 = Pars4(3)*exp(-((v-Pars4(4))/Pars4(5)));

    % potassium4 inactivation
    hk4_inf = 1/((1+exp((v-Pars4(6))/Pars4(7))));
    tau_hk4 = Pars4(8)*exp(-((v-Pars4(9))/Pars4(10)));
    
    % ionic currents;
    Ina = g_na*(mna^3)*hna*(v-Ena);
    Ica = g_ca*mca*hca*(v-Eca);
    Ik1 = g_k1*(mk1^4)*hk1*(v-Ek);
    Ik2 = g_k2*(mk2^4)*(hk2^Pars2(13))*(v-Ek);
    Ik3 = g_k3*(mk3^4)*hk3*(v-Ek);
    Ik4 = g_k4*(mk4^4)*hk4*(v-Ek);
    Ileak = g_leak*(v-El);


    %% Differential equations

    dv=(iapp(ix)-Ina-Ica-Ik1-Ik2-Ik3-Ik4-Ileak)/C;
    dmna=(mna_inf-mna)/tau_mna;
    dhna=(hna_inf-hna)/tau_hna;
    dmca=(mca_inf-mca)/tau_mca;
    dhca=(hca_inf-hca)/tau_hca;
    dmk1=(mk1_inf-mk1)/tau_mk1;
    dhk1=(hk1_inf-hk1)/tau_hk1;
    dmk2=(mk2_inf-mk2)/tau_mk2;
    dhk2=(hk2_inf-hk2)/tau_hk2;
    dmk3=(mk3_inf-mk3)/tau_mk3;
    dhk3=(hk3_inf-hk3)/tau_hk3;
    dmk4=(mk4_inf-mk4)/tau_mk4;
    dhk4=(hk4_inf-hk4)/tau_hk4;
    
    W=(rand-0.5);

    v=v+dv*DT+(sigma*sqrt(DT)*W);
    mna=mna+dmna*DT;
    hna=hna+dhna*DT;
    mca=mca+dmca*DT;
    hca=hca+dhca*DT;
    mk1=mk1+dmk1*DT;
    hk1=hk1+dhk1*DT;
    mk2=mk2+dmk2*DT;
    hk2=hk2+dhk2*DT;
    mk3=mk3+dmk3*DT;
    hk3=hk3+dhk3*DT;
    mk4=mk4+dmk4*DT;
    hk4=hk4+dhk4*DT;
    
    I1(ix) = Ik1;                                                           %Handles used to retrieve channel currents for the entire run
    I2(ix) = Ik2;
    I3(ix) = Ik3;
    I4(ix) = Ik4;
    
    vdiff(ix) = dv;
    vrec(ix)=v;                                                             %Handles used to retrieve variables for the entire run
    na1(ix)=hna;
    na2(ix)=mna;
    ca1(ix)=hca;
    ca2(ix)=mca;
    k11(ix)=hk1;
    k12(ix)=mk1;
    k21(ix)=hk2;
    k22(ix)=mk2;
    k31(ix)=hk3;
    k32(ix)=mk3;
    k41(ix)=hk4;
    k42(ix)=mk4;
    leak(ix)=Ina;

end

%% Creating Figures

%pks = findpeaks(vrec,'minPeakProminence',10,'Annotate','Extents');          %Calculates AP frequency of the run
%freq = length(pks)/10;

%frequency(iters) = freq;
%k41a(iters,:) = k41(:);
%k42a(iters,:) = k42(:);

%end

%time = 0:0.00001:10; %for secs

set(0,'DefaultFigureWindowStyle','docked')                                  %Automatically docks subsequent figures
                                                                            
                                                                            %Phase-plane figures of channel activation (k*1) and inactivation (k*2)
%    figure()
%    subplot(2,2,1)
%    hold off
%    plot(vrec(5:end),k11(5:end),'Color','k','LineStyle','--')
%    hold on
%    plot(vrec(5:end),k12(5:end),'Color',[0.8 0.5 0])
%    title('Shaker')
%    ylabel('Activation')
%    subplot(2,2,2)
%    hold off
%    plot(vrec(5:end),k21(5:end),'Color','k','LineStyle','--')
%    hold on
%    plot(vrec(5:end),k22(5:end),'b')
%    title('Shab')
%    subplot(2,2,3)
%    hold off
%    plot(vrec(5:end),k31(5:end),'Color','k','LineStyle','--')
%    hold on
%    plot(vrec(5:end),k32(5:end),'r')
%    title('Shaw')
%    ylabel('Activation')
%    xlabel('Voltage (mV)')
%    subplot(2,2,4)
%    hold off
%    plot(vrec(5:end),k41(5:end),'Color','k','LineStyle','--')
%    hold on
%    plot(vrec(5:end),k42(5:end),'g')
%    title('Shal')
%    xlabel('Voltage (mV)')

%    figure()
%    plot(vrec,vdiff)
%    xlabel('Voltage (mV)')
%    ylabel('dV/dT')
    
    figure()                                                                %Plots of individual channel currents
    subplot(2,1,1)
    plot(t,vrec,'k');
%    title(['ZT',num2str(ZT),' Frequency ' num2str(freq) ' Hz'])
    ylabel('Voltage (mV)')
    xlim([200 700])
    subplot(2,1,2)
    plot(t,I1,'Color',[0.8 0.5 0])
    hold on
    plot(t,I2,'b')
    plot(t,I3,'r')
    plot(t,I4,'g')
    legend('Shaker','Shab','Shaw','Shal') 
    ylabel('Current(pA)')
    xlabel('Time (ms)')
    xlim([200 700])
    hold off
%
%    figure()                                                                %Standard voltage plot of action potentials
%    plot(t,vrec);
%    subplot(2,1,2)
%    title(['ZT',num2str(ZT),' Frequency ' num2str(freq) ' Hz'])
%    set(gca,'XTick',[0 2000 4000 6000 8000 10000])
%    set(gca,'XTickLabel',[0 2 4 6 8 10])
%    xlabel('Time (s)')
%    ylabel('Voltage (mV)')

%% Returning Outputs

    nOutputs = nargout;
    varargout = cell(1,nOutputs);
if nOutputs == 1
    varargout{1} = vrec;
elseif nOutputs == 2
    varargout{1} = vrec;
    varargout{2} = freq;
elseif nOutputs == 3
    varargout{1} = vrec;
    varargout{2} = I1;
    varargout{3} = I2;
elseif nOutputs == 4
    varargout{1} = vrec;
    varargout{2} = I1;
    varargout{3} = I2;
    varargout{4} = I3;
elseif nOutputs == 5
    varargout{1} = vrec;
    varargout{2} = I1;
    varargout{3} = I2;
    varargout{4} = I3;
    varargout{5} = I4;
end

end