ConnectivityMatrix_2D_n6_1NNonly;

%%
num = n^2;
% % set synapse parameter values
gsyn = 0.05;    %0.02; 
taus = 2;  %1.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set perturbing applied current pulse 

Ip = zeros(num,1);
% Perturbation for Figure 6A
pcells = [1:6,13:18,25:30]; %2 cluster horz stripe
Ip(pcells)=0.02;


% pulse of applied current
Tp = zeros(1,2);
Tp(1) = 1500; % ton
Tp(2) = 1800; % toff

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set time to change connectivity matrix
T_mid=1500;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solve WB network model
T0 = 0; T1 = 8000;

tspan=[T0 T1];

WBftn = @(t,y)RHSWB_phi1_perturb(t, T_mid, y, num, W, W, gsyn, taus, Ip, Tp);

% % set initial conditions
%v0 = zeros(num,1); % + 0.1*rand(num,1);
h0=zeros(num,1); % + 0.1*rand(num,1);
n0=zeros(num,1); % + 0.1*rand(num,1);
s0=zeros(num,1);

% % random initial conditions
%v0=-70*ones(num,1) + 60*rand(num,1);

% set spread initial conditions
%v0 = floor(linspace(-70,0,num))';
%v0 = v0 + 5*rand(num,1);

% n=6 3 cluster solutions
% 3 cluster diagonal stripe positive slope
% temp1=[0, -70, -40, 0, -70, -40 ];
% temp2=[-40, 0, -70, -40, 0, -70];
% temp3=[-70, -40, 0, -70, -40, 0];

% 3  cluster vertical stripe 
% Use this IC for Fig 6A
 temp1=[0, -70, -40, 0, -70, -40 ];
 temp2=[0, -70, -40, 0, -70, -40 ];
 temp3=[0, -70, -40, 0, -70, -40 ];

% 3 cluster horizontal stripe
% Use this IC for Fig 6B
% temp1=[0, 0, 0, 0, 0, 0];
% temp2=[-70, -70, -70, -70, -70, -70];
% temp3=[-40, -40, -40, -40, -40, -40];


v0 = [temp1, temp2, temp3, temp1, temp2, temp3];

% 2 cluster diagonal stripe 
% temp1=[0, -70, 0, -70, 0, -70]; %diagonal stripe tilted left 
% temp2=[-70, 0, -70, 0, -70, 0];

% 2 cluster horizontal stripe
% temp1=[0, 0, 0, 0, 0, 0];
% temp2=[-70, -70, -70, -70, -70, -70];

% v0 = [temp1, temp2, temp1, temp2, temp1, temp2];



v0 = v0';



% v0 = v0 + 2*rand(num,1);
% h0 = h0 + 0.1*rand(num,1);
% n0 = n0 + 0.1*rand(num,1);


ICs = [v0; h0; n0; s0];

[T, sol] = ode45(WBftn, tspan, ICs);
index = 1:num;
v = sol(:, index);
h = sol(:, num+index);
nn = sol(:, 2*num+index);
s = sol(:,3*num+index);

%%
close all
hh = figure(1);
plot(T, v(:,1:4),'LineWidth',2);
axis([T1-400 T1 -90 70]);

%axis([T_mid-400 T_mid+500 -90 70])
set(gca,'fontsize',25,'fontweight','bold')
xlabel('time')
ylabel('v')


%% generate voltage trace plot before and after perturbation
figure(3);
subplot(1,2,1)
plot(T, v(:,1:4),'LineWidth',2);
axis([Tp(1)-400 Tp(1) -90 70]);

%axis([T_mid-400 T_mid+500 -90 70])
set(gca,'fontsize',25,'fontweight','bold')
xlabel('time')
ylabel('v')

subplot(1,2,2)
plot(T, v(:,1:4),'LineWidth',2);
axis([Tp(2)+400 Tp(2)+800 -90 70]);

%axis([T_mid-400 T_mid+500 -90 70])
set(gca,'fontsize',25,'fontweight','bold')
xlabel('time')
ylabel('v')


%% extract spike times for raster plot
% columns: spike time, cell number, row number and even or odd cell numbers
for k=1:num
    [spkht spkind]=findpeaks(v(:,k),'minpeakheight',-10);
    spktimes=T(spkind);
    
    rownum_eo = ceil(k/n);
    if rem(k,2)==0 % check for even cell numbers
       rownum_eo = 2*rownum_eo;
    else
        rownum_eo = 2*rownum_eo-1;
    end

    if k == 1
        spiketimes=horzcat(spktimes,k*ones(length(spktimes),1), rownum_eo*ones(length(spktimes),1));
    else
        spiketimes=[spiketimes; spktimes, k*ones(length(spktimes),1), rownum_eo*ones(length(spktimes),1)];
    end
end


%% plot raster plot
% Figure for manuscript
% odd and even cells in each row are different shades of row color
% even cell numbers are darker shade
figure(2); hold on;
%rastercolor = ['k', 'b', 'c', 'g','r', 'm'];
rastercolor = zeros(12,3);
rastercolor(1,:) = [0.5 0.5 0.5]; % gray
rastercolor(3,:) = [0.1 0.6 0.9330]; % dark blue
rastercolor(4,:) = [0 0 1]; % dark blue
rastercolor(5,:) = [0 1 1]; % cyan
rastercolor(6,:) = [0 0.75 0.75]; % dark cyan
rastercolor(7,:) = [0 1 0]; % green
rastercolor(8,:) = [0 0.6 0.1]; % dark green
rastercolor(9,:) = [1 0 0]; % red
rastercolor(10,:) = [0.75 0 0]; % dark red
rastercolor(11,:) = [1 0 1]; % magenta
rastercolor(12,:) = [0.75 0 0.75]; % dark magenta

for j=1:2*n   % plotting all cells with same color
    clear inds
    inds=spiketimes(:,3)==j;
    plot(spiketimes(inds,1),spiketimes(inds,2),'s','MarkerSize',8,'MarkerFaceColor',rastercolor(j,:), 'MarkerEdgeColor',rastercolor(j,:)) 
end
% axis([T1-200 T1 0 num+1])
axis([Tp(1)-300 Tp(2)+800 0 num+1])

%axis([T_mid-400 T_mid+500 0 num+1])
set(gca,'fontsize',25,'fontweight','bold')
xlabel('time (ms)')
ylabel('No. of Cell')

% add in lines for perturbation
figure(2)
hold on
plot([1500 1500],[0 num+1],'k--');
plot([1800 1800],[0 num+1],'k--');

% raster plot at end of simulation
figure(4); hold on

for j=1:2*n   % plotting all cells with same color
    clear inds
    inds=spiketimes(:,3)==j;
    plot(spiketimes(inds,1),spiketimes(inds,2),'s','MarkerSize',8,'MarkerFaceColor',rastercolor(j,:), 'MarkerEdgeColor',rastercolor(j,:)) 
end
axis([T1-200 T1 0 num+1])
set(gca,'fontsize',25,'fontweight','bold')
xlabel('time (ms)')
ylabel('No. of Cell')
