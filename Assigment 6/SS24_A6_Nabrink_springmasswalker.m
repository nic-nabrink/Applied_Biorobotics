clc;
clear;
close all;

%Parameters
m = 80; %[kg]
g = 9.81; %[m/s^2]
k = 10000; %[N/m]
l0 = 1; %[m]
a0 = 68 * pi/180; %[rad]

% Define initial conditions
y0 = 0.97; %[m]
dy0 = 0;
x0 = 0;
dx0 = 1.05;
fp1 = 0;
td2 = 0;

%Auxiliary fucntions
function fs = fs(x,y,fp,k,l0) %Elastic force
    fs = k*(l0-l(x,y,fp));
end

function cosa = cosa(x,y,fp) %Cos alpha
    cosa = (fp-x)/l(x,y,fp);
end

function sina = sina(x,y,fp) %sin alpha
    sina = y/l(x,y,fp);
end

function l = l(x,y,fp) %size of spring
    l = sqrt((x-fp)^2+y^2);
end

%% EoM

function odedy = odedy(y,dy,x,dx,fp1,fp2,l0,g,m,a0,k)
    if y >= l0*sin(a0) %single support
        odedy = -g + fs(x,y,fp1,k,l0)*sina(x,y,fp1)/m;
    else %double support
        odedy = -g + fs(x,y,fp1,k,l0)*sina(x,y,fp1)/m + fs(x,y,fp2,k,l0)*sina(x,y,fp2)/m;
    end
end

function odedx = odedx(y,dy,x,dx,fp1,fp2,l0,g,m,a0,k)
    if y >= l0*sin(a0) %single support
        odedx = -fs(x,y,fp1,k,l0)*cosa(x,y,fp1)/m;
    else %double support
        odedx = -fs(x,y,fp1,k,l0)*cosa(x,y,fp1)/m - fs(x,y,fp2,k,l0)*cosa(x,y,fp2)/m;
    end
end  

%%

% Initialize the time vector
h = 0.00001; %time step
t = 0:h:2.5;

% Preallocate the solution vectors
y = zeros(1, length(t));
dy = zeros(1, length(t));
x = zeros(1, length(t));
dx = zeros(1, length(t));
e_s = zeros(1, length(t));
e_k = zeros(1, length(t));
e_p = zeros(1, length(t));
e_t = zeros(1, length(t));
f_h1 = zeros(1, length(t));
f_v1 = zeros(1, length(t));
f_h2 = zeros(1, length(t));
f_v2 = zeros(1, length(t));
f_aux = zeros(1, length(t));
y(1) = y0;
dy(1) = dy0;
x(1) = x0;
dx(1) = dx0;

%Auxiliary variables
d=1; %begins in the changing from double to single support
s=0; %single support is off

%This variables are used to attribute each section of forces to legs "left
% and "right"
inverte = zeros(1, length(t));
ciclo = 0;

% Loop over each time step
for i = 1:(length(t) - 1)
    ti = t(i);
    yi = y(i);
    dyi = dy(i);
    xi = x(i);
    dxi = dx(i);

    if yi > l0
        display("Failure mode detected: Both legs on the air.");
        return;
    end

    if yi >= l0*sin(a0) && d==1 %Changing from double to single support
        fp1 = td2;
        s = 1;
        d = 0;

    elseif yi < l0*sin(a0) && s==1 %Changing from single to double support
        td2 = xi + l0*cos(a0);
        s = 0;
        d = 1;

    end 

    %Runge-Kutta 4 Method
    k1_y = h * dyi;
    k1_dy = h * odedy(yi, dyi, xi, dxi, fp1,td2,l0,g,m,a0,k);
    k1_x = h * dxi;
    k1_dx = h * odedx(yi, dyi, xi, dxi, fp1,td2,l0,g,m,a0,k);

    k2_y = h * (dyi + k1_dy/2);
    k2_dy = h * odedy(yi + k1_y/2, dyi + k1_dy/2, xi + k1_x/2, dxi + k1_dx/2, fp1,td2,l0,g,m,a0,k);
    k2_x = h * (dxi + k1_dx/2);
    k2_dx = h * odedx(yi + k1_y/2, dyi + k1_dy/2, xi + k1_x/2, dxi + k1_dx/2, fp1,td2,l0,g,m,a0,k);

    k3_y = h * (dyi + k2_dy/2);
    k3_dy = h * odedy(yi + k2_y/2, dyi + k2_dy/2, xi + k2_x/2, dxi + k2_dx/2, fp1,td2,l0,g,m,a0,k);
    k3_x = h * (dxi + k2_dx/2);
    k3_dx = h * odedx(yi + k2_y/2, dyi + k2_dy/2, xi + k2_x/2, dxi + k2_dx/2, fp1,td2,l0,g,m,a0,k);

    k4_y = h * (dyi + k3_dy);
    k4_dy = h * odedy(yi + k3_y, dyi + k3_dy, xi + k3_x, dxi + k3_dx, fp1,td2,l0,g,m,a0,k);
    k4_x = h * (dxi + k3_dx);
    k4_dx = h * odedx(yi + k3_y, dyi + k3_dy, xi + k3_x, dxi + k3_dx, fp1,td2,l0,g,m,a0,k);

    % Update the solution
    y(i + 1) = yi + (k1_y + 2*k2_y + 2*k3_y + k4_y) / 6;
    dy(i + 1) = dyi + (k1_dy + 2*k2_dy + 2*k3_dy + k4_dy) / 6;
    x(i + 1) = xi + (k1_x + 2*k2_x + 2*k3_x + k4_x) / 6;
    dx(i + 1) = dxi + (k1_dx + 2*k2_dx + 2*k3_dx + k4_dx) / 6;
    
    %Forces and Energies
    if yi > l0*sin(a0) 
        %single support
        inverte(i) = 0;
        f_v1(i) = fs(xi,yi,fp1,k,l0)*sina(xi,yi,fp1);
        f_h1(i) = fs(xi,yi,fp1,k,l0)*cosa(xi,yi,fp1);
        f_v2(i) = 0;
        f_h2(i) = 0;
        if i>1 && inverte(i)~=inverte(i-1) %inverte and ciclo are used to assign forces to the left and right legs
            if ciclo == 1
                ciclo = 0;
        	elseif ciclo == 0
                ciclo = 1;
            end
        end
        if ciclo == 1
            f_v2(i) = fs(xi,yi,fp1,k,l0)*sina(xi,yi,fp1);
            f_h2(i) = fs(xi,yi,fp1,k,l0)*cosa(xi,yi,fp1);
            f_v1(i) = 0;
            f_h1(i) = 0;
        end

        e_s(i) = 0.5*k*(l(xi,yi,fp1)-l0)^2;
        e_p(i) = m*g*yi;
        e_k(i) = 0.5*m*(dxi^2+dyi^2);
        e_t(i) = e_s(i) + e_p(i) + e_k(i);
    else
        %double support
        f_v1(i) = fs(xi,yi,fp1,k,l0)*sina(xi,yi,fp1);
        f_h1(i) = fs(xi,yi,fp1,k,l0)*cosa(xi,yi,fp1);
        f_v2(i) = fs(xi,yi,td2,k,l0)*sina(xi,yi,td2);
        f_h2(i) = fs(xi,yi,td2,k,l0)*cosa(xi,yi,td2);
        if ciclo == 1
            f_v2(i) = fs(xi,yi,fp1,k,l0)*sina(xi,yi,fp1);
            f_h2(i) = fs(xi,yi,fp1,k,l0)*cosa(xi,yi,fp1);
            f_v1(i) = fs(xi,yi,td2,k,l0)*sina(xi,yi,td2);
            f_h1(i) = fs(xi,yi,td2,k,l0)*cosa(xi,yi,td2);
        end
        e_s(i) = 0.5*k*((l(xi,yi,fp1)-l0)^2 + (l(xi,yi,td2)-l0)^2);
        e_p(i) = m*g*yi;
        e_k(i) = 0.5*m*(dxi^2+dyi^2);
        e_t(i) = e_s(i) + e_p(i) + e_k(i);
        inverte(i) = 1;
    end
end
        e_s(i+1) = e_s(i);
        e_p(i+1) = e_p(i);
        e_k(i+1) = e_k(i);
        e_t(i+1) = e_t(i);

%Plot
figure(1);
plot(x, y, 'LineWidth', 2, 'Color', 'b');
title('Position of CoM over time');
ylim([0.85 1]);
xlabel('$x_{CoM} [m]$','Interpreter','Latex', 'FontSize', 14);
ylabel('$y_{CoM} [m]$','Interpreter','Latex', 'FontSize', 14);
grid on;

figure(2)
plot(t, -f_h1, 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Left leg')
hold on
plot(t, -f_h2, 'LineWidth', 2, 'Color', 'black', 'DisplayName', 'Right leg')
grid on
hold off
xlabel('t [s]', 'FontSize', 14)
ylabel('$F_{horiz} [N]$','Interpreter','Latex', 'FontSize', 14)
title("Horizontal elastic forces over time", 'FontSize', 15)
legend('FontSize', 11)

figure(3)
plot(t, f_v1, 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Left leg')
hold on
plot(t, f_v2, 'LineWidth', 2, 'Color', 'black', 'DisplayName', 'Left leg')
hold off
xlabel('t [s]', 'FontSize', 14)
ylabel('$F_{vert} [N]$','Interpreter','Latex', 'FontSize', 14)
title("Vertical elastic forces over time")
legend('FontSize', 11)
grid on

figure(4)
plot(t,e_t , 'LineWidth', 2, 'Color', 'b')
ylim([400 1000])
xlabel('t [s]', 'FontSize', 14)
ylabel('E [J]', 'Interpreter', 'Latex', 'FontSize', 14) 
title("Total system energy over time")
grid on

figure(5)
plot(y,dy , 'LineWidth', 2, 'Color', 'b')
xlabel('y [m]', 'FontSize', 14)
ylabel('$\dot y [m/s]$', 'Interpreter', 'Latex', 'FontSize', 14) 
title("Total system energy over time")
grid on

figure(6)
meandx = mean(dx);
plot(t,dx , 'LineWidth', 2, 'Color', 'b')
hold on
yline(meandx, 'r--', 'LineWidth', 2);
xlabel('t [s]', 'FontSize', 14)
ylabel('$\dot y [m/s]$', 'Interpreter', 'Latex', 'FontSize', 14)
title("Forward velocity of the system")
text(1, meandx + 0.005, ['Mean = ', num2str(meandx)], 'Color', 'red', 'FontSize', 12);
legend('Forward Velocity', 'Mean Value')
hold off
grid on