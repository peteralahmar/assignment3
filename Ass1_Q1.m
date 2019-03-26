%Assignment 1 Peter Al-Ahmar, 100961570
%Question 1

clear;
close all;


q = 1.60217653e-19;             % electron charge
hb = 1.054571596e-34;           % Dirac constant
h = hb * 2 * pi;                % Planck constant
m = 9.10938215e-31;             % electron mass
kb = 1.3806504e-23;             % Boltzmann constant
eps = 8.854187817e-12;          % vacuum permittivity
mu = 1.2566370614e-6;           % vacuum permeability
c = 299792458;                  % speed of light
g = 9.80665;                    %metres (32.1740 ft) per s²
am = 1.66053892e-27;

iter = 1000;
T=300;

AtomSpacing = 0.5430710e-9;

nParticles = 5000;

em= 0.26*m; %effective mass

vth = sqrt(kb*T/em);

strdDev = vth/sqrt(2);

xboundary = 200e-9;
yboundary = 200e-9;
area = xboundary * yboundary;
concentration = 10e-15;

volts = 0.1;
ex = volts/xboundary;
ey = 0;

initialx = rand(1,nParticles) .*xboundary;
initialy = rand(1,nParticles) .*yboundary;

dt = 2* 10e-15;

temparray = zeros(1,1000); % temp array
temp = (1:1:1000);

%7 random electrons to be chosen
chosen = randi([1, nParticles], 1, 7);

inside_box = true;
while inside_box == true
    inside = ((initialx <= (1.15 * xboundary/2) & (initialx >= (0.85 * xboundary/2))) & ((initialy < (yboundary/3) | initialy >= (2*yboundary/3))));
    
    if (sum(inside) >0)
        initialx(inside) = rand(1,sum(inside)).* xboundary;
        initialy(inside) = rand(1,sum(inside)) .* yboundary;
    else
        inside_box = false;
    end
end

velocityx = rand(1,nParticles) .*strdDev;
velocityy = rand(1,nParticles) .*strdDev;


vrms = sqrt((velocityx .^ 2) + (velocityy .^ 2));
vrms_array =zeros(1,1000);

pScatter = 1- (exp((-1*dt) / (0.2e-12))); %probability of scatter 


current_Vector = zeros(1,iter);
time_Vector = 1:1:iter;

%Accelerations
accelerationx = zeros(1,nParticles);
accelerationy = zeros(1,nParticles);

accelerationx(1,:) = ex * (q/em);

D_Velocity = zeros(1,nParticles); % Drift Velocity
E_mobility = zeros(1,nParticles); % Electron Mobility

for i=1: iter
    
    %x-boundary using periodic boundary condition
    xbn = initialx > 200e-9;
    xbn2 = initialx < 0;
    
    initialx(xbn) = initialx(xbn) - 200e-9;
    initialx(xbn2) = initialx(xbn2) + 200e-9;
    
    
    %y-boundary reflection
    
    ybn = initialy >= 200e-9;
    ybn2 = initialy <= 0;
    
    velocityy(ybn) = - velocityy(ybn);
    velocityy(ybn2) = - velocityy(ybn2);
    
    oldx = initialx;
    oldy = initialy;
    
    velocityx = velocityx + (accelerationx .* dt);
    %update position 
    
    initialx= oldx + (velocityx .* dt);
    initialy= oldy + (velocityy .* dt);

    vrms = sqrt((velocityx .^ 2) + (velocityy .^ 2));
    
    temperature = (sqrt(2) * (mean(vrms)) ^2 *em) / kb; 
    temparray(1,i) = temperature;
    
    %electron mobility
    E_mobility = mean(vrms);
    
    %Drift Velocity
    D_Velocity = E_mobility * ex;
    
    %current calculation: charge * area * concentration * average drift velocity 
    current = concentration * area * (sum(D_Velocity)/nParticles) * q; 
    
    current_Vector(i) = current;
%     time_vector(i) = i;
    
    %scatteting with velocities
    Scattering = pScatter > rand(1,nParticles);
    
    velocityx(Scattering) = randn .* strdDev;
    velocityy(Scattering) = randn .* strdDev;
    
    figure(1);
    plot(initialx(chosen),initialy(chosen),'.');
    xlabel("x-Position");
    ylabel("y-Position");

    title(['Average Temperature: ', num2str(temperature)]);

    xlim([0,xboundary]);
    ylim([0,yboundary]);
    hold on;   
end

%plot of Temperature vs time

figure(2);
plot(temp,temparray)
title('Average Temperature')
xlabel('Time (s)')
ylabel('Temperature')
hold on

%Plot of Current vs Time
figure(3);
plot(time_Vector, current_Vector)
title('Average Current')
xlabel('Time')
ylabel('Current')
hold on

[xgr,ygr] = meshgrid(0:(xboundary/30):xboundary, 0:(yboundary/30):yboundary);
electron_M = zeros(30,30);
temperature_M = zeros(30,30);
num_Elec = 0;
Total_vel = 0;

for k = 1:30
    x_min = xgr(1,k);
    x_max = xgr(1,k+1);
    
    for m = 1:30    
        y_min = ygr(m,1);
        y_max = ygr(m+1,1);
        
        for n = 1:nParticles
            if ((initialx(n) >x_min) && (initialx(n) < x_max) && ((initialy(n) > y_min) && initialy(n) <y_max))
                num_Elec = num_Elec +1;
                electron_M(k,m) = electron_M(k,m) + 1;
                Total_vel = Total_vel + sqrt((velocityx(n) .^2) + (velocityy(n) .^2));
                if(num_Elec ~=0)
                temperature_M(k,m) = ((sqrt(2)*(Total_vel/num_Elec)^2)*em)/kb;
                end
            end
             
        end
        Total_vel = 0;
        num_Elec = 0;
        
    end
    
end
    
% Creating an Electron Map
figure(4);
surf(flipud(electron_M));
title('Electron Density Map');
    
%Creating a Temperature Map
figure(5);
surf(flipud(temperature_M));
title('Temperature Map');



