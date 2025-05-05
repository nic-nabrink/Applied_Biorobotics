% Define time vector
t = linspace(0, 10, 10000); % 10000 time points from 0 to 10 seconds

%Parameters
g = 9.81; %[m/s^2]
k = 100; %[N/m]
len = 1; %[m]
h = 1.2; %[m]
m = 1; %[kg]
w = sqrt(k/m);

% Define initial conditions
y0 = h;
v0 = 0;

%% Simulate system
y = zeros(size(t));
v = zeros(size(t));
e_s = zeros(size(t));
y(1) = y0;
v(1) = v0;
% Using a auxiliary time system to each phase of the dynamics
cont_stance = 0;
cont_up = 0;
cont_down = 0;
t_aux = 0;
dt = t(2) - t(1);

for i = 2:length(t)
    % Down
    if y(i-1) > 1 && v(i-1) <= 0
        cont_up = 0;
        cont_stance = 0;
        t_aux = cont_down * dt;
        y(i) = -(g*t_aux^2)/2 + h;
        v(i) = -g*t_aux;
        cont_down = cont_down + 1;

    % Up
    elseif y(i-1) >= 1 && v(i-1) > 0
        cont_down = 0;
        cont_stance = 0;
        t_aux = cont_up * dt;
        y(i) = -(g*t_aux^2)/2 + sqrt(2*(h-len)*g)*t_aux + len;
        v(i) = -g*t_aux + sqrt(2*(h-len)*g);
        cont_up = cont_up + 1;
    
    %Stance
    else
        cont_up = 0;
        cont_down = 0;
        t_aux = cont_stance * dt;
        y(i) = (m*g/k)*cos(w*t_aux) - sqrt(2*(h-len)*g)/w*sin(w*t_aux) + len - m*g/k;
        v(i) = -(m*g/k)*w*sin(w*t_aux) - sqrt(2*(h-len)*g)*cos(w*t_aux);
        cont_stance = cont_stance + 1;

    end
    e_s(i) = (y(i) <= 1)*0.5*k*(1-y(i))^2;
end

%%

%% Energy
e_k = 0.5*m*v.^2;
e_p = m*g*y;
total = e_k + e_p + e_s;
%%


% Plot results
figure(1);
subplot(2, 1, 1);
plot(t, y);
title('Position vs Time');
xlabel('Time (s)');
ylabel('Position (m)');
grid on;

subplot(2, 1, 2);
plot(t, v, 'Color','r');
title('Velocity vs Time');
xlabel('Time (s)');
ylabel('Velocity (m/s)');
grid on;

figure(2);
plot(t, e_p);
hold on;
plot(t, e_k, 'Color', 'g');
hold on
plot(t, e_s);
hold on
plot(t, total);
title('Energy vs Time');
xlabel('Time (s)');
ylabel('Energy (J)');
legend('e_p', 'e_k', 'e_s', 'Total');
grid on;
