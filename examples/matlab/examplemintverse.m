%%	An example of minimum-time VERSE applied to excitation.

FA     = 90;       % °
RFdur  = 2.56;     % ms
maxB1  = 23 *1e-6; % T
maxG   = 80 *1e-3; % T/m
maxSR  = 200;      % T/m/s
dt     = 10;       % µs
TBWP   = 10;
offset = 0;        % mm
slthk  = 5;        % mm
units  = 'T';      % 'Hz', 'T', 'G'

% Generate Sinc pulse and gradient
[B1, Gmax] = gensinc(FA, 0, RFdur, TBWP/2, TBWP/2, 'Hann', dt, units, slthk, offset); % [T, T/m]

% Create gradient option 1
rampPoints = ceil(1e6*Gmax/maxSR/dt);
G = Gmax * [(0.5:rampPoints)'/rampPoints ; ones(size(B1)) ; (rampPoints-0.5:-1:0.5)'/rampPoints];
B1 = [zeros(rampPoints,1); B1; zeros(rampPoints,1)];

% % Create gradient option 2
% % If waveforms start and end with zero, so will the VERSEd versions
% % (This will create the gradient ramps for you)
% G  = Gmax * ones(size(B1));
% B1 = [0;B1;0];
% G  = [0;G ;0];

tic
% [B1v,Gv] = mintverse(B1,G,dt*1e-6,maxB1,maxG,maxSR); % Native MATLAB
[B1v,Gv] = mintverse_c(B1,G,dt*1e-6,maxB1,maxG,maxSR); % MEXed C-code
toc

% Now, repeat but constrain the energy to 50% of the original pulse.
eb1 = sum(abs(B1).^2)*dt*1e-6;
tic
% [B1ve,Gve] = mintverse(B1,G,dt*1e-6,maxB1,maxG,maxSR,0.5*eb1); % Native MATLAB
[B1ve,Gve] = mintverse_c(B1,G,dt*1e-6,maxB1,maxG,maxSR,0.5*eb1); % MEXed C-code
toc

% Time vectors, for plotting.
t   = (0.5:length(B1 ))'  * dt*1e-6;
tv  = (0.5:length(B1v))'  * dt*1e-6;
tve = (0.5:length(B1ve))' * dt*1e-6;

% RF Energy (time integral of B1^2)
eb1   = sum(abs(B1  ).^2) * dt*1e-6;
eb1v  = sum(abs(B1v ).^2) * dt*1e-6;
eb1ve = sum(abs(B1ve).^2) * dt*1e-6;

figure('Units','normalized','Position',[0,0,1,1]);

tlim  = max([t(end),tv(end),tve(end)]);
B1lim = 1e6*max(abs([B1;B1v;B1ve]));  % μT
Glim  = 1e3*1e2*max(abs([G;Gv;Gve])); % mT/m
if Glim==0; Glim = 1; end

subplot(2,3,1); hold on
plot(t,1e6*abs(B1),'LineWidth',2);
plot(t,1e6*real(B1),'LineWidth',2);
plot(t,1e6*imag(B1),'LineWidth',2);
set(gca,'FontSize',16);
text(t(5),-0.9*B1lim,['Energy = ' num2str(eb1) ' T^2s'],'FontSize',20);
title('Initial B_1^+','FontSize',20);
xlabel('Time (s)','FontSize',20);
ylabel('B1 (µT)','FontSize',20);
xlim([0,tlim]);
ylim([-B1lim,B1lim]);
grid on; box on;

subplot(2,3,2); hold on
plot(tv,1e6*abs(B1v),'LineWidth',2);
plot(tv,1e6*real(B1v),'LineWidth',2);
plot(tv,1e6*imag(B1v),'LineWidth',2);
set(gca,'FontSize',16);
text(t(5),-0.9*B1lim,['Energy = ' num2str(eb1v) ' T^2s'],'FontSize',20);
title('VERSE B_1^+','FontSize',20);
xlabel('Time (s)','FontSize',20);
ylabel('B1 (µT)','FontSize',20);
xlim([0,tlim]);
ylim([-B1lim,B1lim]);
legend({'abs(B1)','real(B1)','imag(B1)'},'FontSize',20);
grid on; box on;

subplot(2,3,3); hold on
plot(tve,1e6*abs(B1ve),'LineWidth',2);
plot(tve,1e6*real(B1ve),'LineWidth',2);
plot(tve,1e6*imag(B1ve),'LineWidth',2);
set(gca,'FontSize',16);
text(t(5),-0.9*B1lim,['Energy = ' num2str(eb1ve) ' T^2s'],'FontSize',20);
title('Energy-Constrained VERSE B_1^+','FontSize',20);
xlabel('Time (s)','FontSize',20);
ylabel('B1 (µT)','FontSize',20);
xlim([0,tlim]);
ylim([-B1lim,B1lim]);
grid on; box on;

subplot(2,3,4); hold on
plot(t,1e3*G,'LineWidth',2);
set(gca,'FontSize',16);
title('Initial Gradient','FontSize',20);
xlabel('Time (s)','FontSize',20);
ylabel('Gradient (mT/m)','FontSize',20);
xlim([0,tlim]);
ylim([-Glim,Glim]);
grid on; box on;

subplot(2,3,5); hold on
plot(tv,1e3*Gv,'LineWidth',2);
set(gca,'FontSize',16);
title('VERSE Gradient','FontSize',20);
xlabel('Time (s)','FontSize',20);
ylabel('Gradient (mT/m)','FontSize',20);
xlim([0,tlim]);
ylim([-Glim,Glim]);
grid on; box on;

subplot(2,3,6); hold on
plot(tve,1e3*Gve,'LineWidth',2);
set(gca,'FontSize',16);
title('Energy-Constrained VERSE Gradient','FontSize',20);
xlabel('Time (s)','FontSize',20);
ylabel('Gradient (mT/m)','FontSize',20);
xlim([0,tlim]);
ylim([-Glim,Glim]);
grid on; box on;
