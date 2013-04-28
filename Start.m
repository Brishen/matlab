% eps=8.433E-20 % J
% sig=2.321 %Angstroms = 2.321E-10 m
% alpha=2.867 %Angstroms = 2.866E-10 m


% column 1: initial x position
% column 2: initial y position
% column 3: initial x velocity
% column 4: initial y velocity
% column 5: mass

% PARTICLE BOX SIZE 11x5
% L0=10sig

close all
clear all
clc

%Absolute Units
boltzmann=1.38065E-23;          %J/K
epsilon=8.433E-20;              %J
sigma=2.866E-10;                %M
strainrate=10^9;                %sec^-1
M=9.273E-23;                    %g, atom mass

%Reduced Units, divide by factors to get dimensionless value
TempR=6801;                     %Factor for Temperature
Tempr=0.0491;                   
V1r=0.313369;                   %Average speed of all particles, random direction.
V2r=0.003005453;                %Velocity of particles moving to the right, constant in x direction
VR=953.6;                       %Factor for velocity
kbR=1;                          %Benefit of Reduced Units!

%% Because of Reduced Units, we can simplify the LJ potential.
% the force acting on a particle i is now described by
% f(i,j)=48*(r(i,j)^(-14)-0.5*r(i,j)^(-8))rcoeff(i,j)
% r(i,j)=sqrt( (x(i)-x(j))^2+(y(i)-y(j))^2)



%Assigns Initial Positions on a grid with spacing Sep
Nx=11;      %So L0=10
Ny=5;
Sep=1;      %Reduced, equal to 1 sigma

n=1:Sep:(Nx)*Sep;
m=1:Sep:(Ny)*Sep;
l_0=(Nx-1)*Sep;

for i=1:Nx
    for j=1+(i-1)*Ny:i*Ny
        x(j,1)=n(i);
        y(j,1)=m(j-((i-1)*Ny));
    end
end

%To set initial velocitys, we will have each atom moving in opposite
%directions so that it will sum to 0. the left most atoms have v=0, but we
%have to account for that in the temperature calculations. Vy=0 in all
%cases

v(1:1:Ny,1)=0;             %Script switches u and v because it is stupid
v(6:1:10,1)=-0.003005453;   %This is different to account for the right end
v(11:1:15,1)=0.313369;
v(16:1:20,1)=-0.313369;
v(21:1:25,1)=0.313369;
v(26:1:30,1)=-0.313369;
v(31:1:35,1)=0.313369;
v(36:1:40,1)=-0.313369;
v(41:1:45,1)=0.313369;
v(46:1:50,1)=-0.313369;
v(51:1:55,1)=0.003005453;

u(1:1:Nx*Ny,1)=0;

m=0;
m(1:1:Nx*Ny,1)=1;

N=55;
variables=5;

clearvars -except x y v u m variables N             %Just to be careful, remove all but needed
lj
