%%	An example of minimum-SAR VERSE applied to excitation.

FA     = 90;             % °
RFdur  = 2.56;           % ms
maxG   = 80  *1e-3*1e-2; % T/cm
maxSR  = 200 *1e-2;      % T/cm/s
dt     = 2;             % µs
TBWP   = 10;
offset = 0;              % mm
slthk  = 5;              % mm
units  = 'T';            % 'Hz', 'T', 'G' ()

[B1, Gmax] = gensinc(FA, 0, RFdur, TBWP/2, TBWP/2, 'Hann', dt, units, slthk, offset); % [T, T/cm]

% Create gradient option 1
rampPoints = 1e6*Gmax/maxSR/dt;
G = Gmax * [(0.5:rampPoints)'/rampPoints ; ones(size(B1)) ; (rampPoints-0.5:-1:0.5)'/rampPoints];
B1 = [zeros(rampPoints,1); B1; zeros(rampPoints,1)];

% % Create gradient option 2
% % If waveforms start and end with zero, so will the VERSEd versions
% % (This will create the gradient ramps for you)
% G  = Gmax * ones(size(B1));
% B1 = [0;B1;0];
% G  = [0;G ;0];

tic
% [B1v, Gv] = minsarverse_matlab_simpler(B1, G, dt*1e-6, maxG, maxSR); % Native MATLAB
% [B1v, Gv] = minsarverse_matlab(B1, G, dt*1e-6, maxG, maxSR); % Native MATLAB
[B1v, Gv] = minsarverse(B1, G, dt*1e-6, maxG, maxSR); % MEXed C-code
toc

% Change units for plotting
B1p  = B1  * 1e6;     % Convert to µT
B1vp = B1v * 1e6;     % Convert to µT
Gp   = G   * 1e2*1e3; % Convert to mT/m
Gvp  = Gv  * 1e2*1e3; % Convert to mT/m

% Time vectors, for plotting.
t  = (0.5:length(B1p ))' *dt*1e-6;
tv = (0.5:length(B1vp))' *dt*1e-6;

figure;

B1lim = max(abs([B1p;B1vp]));
Glim  = max(abs([Gp;Gvp]));
if Glim==0; Glim = 1; end

subplot(2,2,1); hold on
plot(t,abs(B1p),'LineWidth',2);
plot(t,real(B1p),'LineWidth',2);
plot(t,imag(B1p),'LineWidth',2);
grid on;
title('Starting B1 pulse vs t');
xlabel('Time (s)');
ylabel('B1 (µT)');
xlim([t(1),t(end)]);
ylim([-B1lim,B1lim]);

subplot(2,2,2); hold on
plot(tv,abs(B1vp),'LineWidth',2);
plot(tv,real(B1vp),'LineWidth',2);
plot(tv,imag(B1vp),'LineWidth',2);
grid on;
title('VERSE B1 pulse vs t');
xlabel('Time (s)');
ylabel('B1 (µT)');
xlim([t(1),t(end)]);
ylim([-B1lim,B1lim]);
legend('abs(B1)','real(B1)','imag(B1)');

subplot(2,2,3); hold on
plot(t,Gp,'LineWidth',2);
grid on;
title('Starting Gradient pulse vs t');
xlabel('Time (s)');
ylabel('Gradient (mT/m)');
xlim([t(1),t(end)]);
ylim([-Glim,Glim]);

subplot(2,2,4); hold on
plot(tv,Gvp,'LineWidth',2);
grid on;
title('VERSE Gradient pulse vs t');
xlabel('Time (s)');
ylabel('Gradient (mT/m)');
xlim([t(1),t(end)]);
ylim([-Glim,Glim]);
