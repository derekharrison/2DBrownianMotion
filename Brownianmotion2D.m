% Simulation of two-dimenisonal brownian motion. 

clear 
clc
tic

Ns = 5;                             %Number of simulations
Nw = 50;                            %Population group number 
Dacc=0;
simc=1;
sdsum=0;
sdxsum=0;
simulatedD2sum=0;
Xf=zeros(Ns,1);
upperbound = 1;
lowerbound = -1;
stdev = 1;                          %Standard deviation, 1 for standard normal distribution
delx = (upperbound-lowerbound)/Nw;  %Population group sample width
xrange=zeros(Nw,1);
N=4000;
ti=0;
tf=1000;
dimensions = 2;             %two dimensional simulation
tau = (tf-ti)/N;            % time interval in seconds
d    = 1.0e-8;              % radius in meters
eta  = 1.0e-3;              % viscosity of water in SI units (Pascal-seconds) at 293 K
kB   = 1.38e-23;            % Boltzmann constant
T    = 293;                 % Temperature in degrees Kelvin

D    = kB * T / (3 * pi * eta * d);

writerObj = VideoWriter('C:\Users\dharrison\Desktop\Simulations','Motion JPEG AVI');
writerObj.FrameRate = 48*3;
open(writerObj);

for i=1:Ns

k = sqrt(D * dimensions * tau);
dx = k * randn(N,1);
dy = k * randn(N,1);

x = cumsum(dx);
y = cumsum(dy);
Xf(i)=x(end);

dSquaredDisplacement = (dx .^ 2) + (dy .^ 2);
squaredDisplacement = ( x .^ 2) + ( y .^ 2);

simulatedD = mean(dSquaredDisplacement)/2/dimensions/tau;
simulatedD2sum  = simulatedD2sum+squaredDisplacement(end);
standardError = std( dSquaredDisplacement ) / ( 2 * dimensions * tau * sqrt(N) );
actualError = D - simulatedD;
Dacc=Dacc+simulatedD;
meanD=Dacc/simc;
sdsum=sdsum+squaredDisplacement(end);
sdxsum=sdxsum+x(end).^2;
meansd=sdsum/simc;
meansdx=sdxsum/simc;
meansimulatedD2=simulatedD2sum/simc/2/dimensions/tf;

for j=1:N
 plot(x(1:j),y(1:j),x(j),y(j),'r.','MarkerSize',20)
 title('Particle Track of a Single Simulated Particle');
minY = -65e-5;
maxY = 65e-5;
xlabel('x-coordinate')
ylabel('y-coordinate')
set(gca,'Ylim',[minY maxY])
set(gca,'Xlim',[minY maxY])
set(gca,'nextplot','replacechildren');
set(gcf,'Renderer','zbuffer')
set(gcf,'units','normalized','outerposition',[0.2 0.10 0.56 1.6*0.56])
   frame = getframe(gcf);
   writeVideo(writerObj,frame);
simc=simc+1;
end
end
close(writerObj);

    Ngroups = zeros(Nw,1);              %Population groups array
    Ngfreq = zeros(Nw,1);               %Frequency distribution
    


    xrange(1,1) = lowerbound+delx/2;
    i=1;
    while xrange(i,1)<upperbound-delx
        xrange(i+1,1) = xrange(i,1)+delx;
        i=i+1;
    end
    
    for i=1:Ns
        for j=1:Nw
            if Xf(i)>xrange(j,1)-delx/2 && Xf(i) <= xrange(j,1)+delx/2
                Ngroups(j,1)=Ngroups(j,1)+1;
            end 
        end
    end
    
        for j=1:Nw
        Ngfreq(j,1) = Ngroups(j,1);
        end
    
       % Frequency distribution analytical
    Ngfreqana = 1/sqrt(4*pi*D*tf)*exp(-(xrange).^2./2/dimensions/D/tf)*delx*Ns;
   
        %Averageprobability density obtained from simulations up until now
   
    plot(xrange,Ngfreq,xrange,Ngfreqana)
    
Dsimulated=meanD
Dt4 = 4*D*tf
Dt2 = 2*D*tf
meansd=meansd
var=meansd
meansdx
varx=meansdx
D
meansimulatedD2
clear 

toc