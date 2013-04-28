%This script plots a bunch of particles and their behavior
%The script can read both 1D and 2D files

%%% DATA TO CHANGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%filename = 'input2d.txt'; %name of the data file which should be loaded
sigma = 1;              %particle diameter
epsilon = 1;            %material parameter
dt = .005;        %time step in seconds
time = 1;              %total time in seconds
method = 'Verlet';      %integration method
borders = [0 30 -5 10];           %set [x0 x1 y0 y1] to turn borders on (particles bounce), 0 is off
output = '';            %name of the output file. Set '' for no output
movie = 'output';       %movie output name. Set '' for no movie output
fps = 20;               %FPS for movie
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b=zeros(size(borders));      %border comparator
timesteps = round(time/dt+1); %number of timesteps to be calculated
%[x, y, v, u, m, variables, N] = loadfile( filename, timesteps ); %load the data

V = zeros(N,timesteps); %allocate potential energy array
T = zeros(N,timesteps); %allocate kinetic energy array
f_x = zeros(1,N);       %allocate force matrix x-direction
f_y = zeros(1,N);       %allocate force matrix y-direction
V_j = zeros(1,N);       %allocate LennardJones potential matrix

%%% Calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for t=1:timesteps
    for i=6:N
        for j=6:N %set force/energy interactions for each particle
            if i~=j %don't take self inducting terms into account
           
                %compute the distance
                rij = sqrt( (x(i,t)-x(j,t))^2  +  (y(i,t)-y(j,t))^2  );
                
                %compute unit normal vector, n(1)=n_x, n(2)=n_y
                n = [ x(i,t)-x(j,t) y(i,t)-y(j,t) ]/norm([ x(i,t)-x(j,t) y(i,t)-y(j,t) ]);
                
                %for x direction
                %f_x(j) = n(1)*-4*epsilon*((6*sigma^6)/rij^7 - (12*sigma^12)/rij^13);
                f_x(j) = n(1)*      -48*      (0.5*rij^(-6)    -      rij^(-12));
                %for y direction
                %f_y(j) = n(2)*-4*epsilon*((6*sigma^6)/rij^7 - (12*sigma^12)/rij^13);
                f_y(j) = n(2)*      -48*      (0.5*rij^(-6)    -      rij^(-12));

                %potential energy
                %V_j(j) = 4*epsilon*((sigma/rij)^12-(sigma/rij)^6);
                V_j(j) = 4*(   rij^(-12)   -   rij^(-6)   );
                
                
                
            else %if i=j: force=0, potential=0
                f_x(j) = 0;
                f_y(j) = 0;
                V_j(j) = 0;
            end
        end
        
        %Compute the position and velocity using Euler or Verlet
        if (strcmp('Euler',method) || t==1) %first timestap is always Euler   
            % Forward Euler x
            v(51:1:55,t)=0.003005453;
            x(i,t+1) = x(i,t) + v(i,t)*dt;
            v(i,t+1) = v(i,t) + sum(f_x)/m(i)*dt;
            
            % Forward Euler y
            y(i,t+1) = y(i,t) + u(i,t)*dt;
            u(i,t+1) = u(i,t) + sum(f_y)/m(i)*dt;
            
            %BC
            v(51:1:55,t)=0.003005453;
            %v(1:1:5,t)=0;
            %u(1:1:5,t)=0; 
            %u(1:1:55,t)=0;
            x(1:1:5,t)=1;
            for xyz=1:1:5
                y(xyz,t)=xyz;
                y(xyz+50,t)=xyz;
            end
            
            
        elseif(strcmp('Verlet',method))
            v(51:1:55,t)=0.003005453;
            %Verlet algorithm x
            x(i,t+1) = -x(i,t-1) + 2*x(i,t) + sum(f_x)/m(i)*dt^2;
            v(i,t) = (x(i,t-1) - x(i,t+1))/(2*dt);
            %Verlet algorithm y
            y(i,t+1) = -y(i,t-1) + 2*y(i,t) + sum(f_y)/m(i)*dt^2;
            u(i,t) = (y(i,t-1) - y(i,t+1))/(2*dt);
            
            %BC
            v(51:1:55,t)=0.003005453;
            %v(1:1:5,t)=0;
            %u(1:1:5,t)=0; 
            %u(1:1:55,t)=0;
            x(1:1:5,t)=1;
            for xyz=1:1:5
                y(xyz,t)=xyz;
                y(xyz+50,t)=xyz;
            end
            
        else
            disp('No valid integration method given!');
        end
        
        %%% If borders are on %%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        if isequal(borders,b) == 0
            if (x(i,t+1) < borders(1)+sigma/2) || (x(i,t+1) > borders(2)-sigma/2) %x
                x(i,t+1) = 2*x(i,t) - x(i,t+1);
                if (strcmp('Euler',method))
                    v(i,t+1) = -v(i,t+1);
                elseif(strcmp('Verlet',method))
                    v(i,t) = -v(i,t);
                end     
            end
            if (y(i,t+1) < borders(3)+sigma/2) || (y(i,t+1) > borders(4)-sigma/2) %y
                y(i,t+1) = 2*y(i,t) - y(i,t+1);
                if (strcmp('Euler',method))
                    u(i,t+1) = -u(i,t+1);
                elseif(strcmp('Verlet',method))
                    u(i,t) = -u(i,t);
                end   
            end
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %Energy calculations
        V(i,t) = sum(V_j); %Potential energy: sum of influence of other particles
        T(i,t) = 1/2*m(i)*(v(i,t)^2+u(i,t)^2); %Kinetic energy: 1/2mv^2
    end
end

%%% Calculate and plot energy %%%%%%%%%%%%%%%%%%%%%%%%
V_sum = sum(V)/2; %sum the potential energy of all particles times a half since matrix is mirrored
T_sum = sum(T) ;%1/2.*m'.*(sum(u.^2)+sum(v.^2)); %sum from all particles: 0.5*m*v^2

fig=figure(1);
hold on
plot(0:dt:time,T_sum+V_sum)
plot(0:dt:time,T_sum,'r')
plot(0:dt:time,V_sum,'k')
title('The total energy of the system','FontSize',14);
xlabel('Time [s]','FontSize',12);
ylabel('Energy [J]','FontSize',12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Create file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp('', output) %check if output file is provided
    createfile(x, y, v, u, output);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Set the axis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isequal(borders,b) == 0 %if borders are given set the axis equal to the borders
    ax_x = borders(1:2);
    ax_y = borders(3:4);
else %if not set the axis automatically to the highest traveled distance of the particles
    ax_x = [min(min(x))-0.5*sigma max(max(x))+0.5*sigma]; %x axis
    if variables == 3
        ax_y = [-1 1]; %y axis for 1D
    else
        ax_y = [min(min(y))-0.5*sigma max(max(y))+0.5*sigma]; %y axis for 2D
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%% Draw the circles and make the movie%%%%%%%%%%%%%%%
if ~(strcmp('',movie))
    fig=figure(2);
    mov = VideoWriter(movie, 'MPEG-4');
    mov.FrameRate = fps;
    open(mov);
    %for k=1:round(1/dt/fps):timesteps %pick the points for the correct number of fps
     for k=1:round(1/dt/fps):timesteps %pick the points for the correct number of fps
        circle = linspace(0,2*pi,100); %create a circle
        hold on;
        for n=1:N %plot N circles...
            xx = sigma/2*cos(circle)+x(n,round(k)); %x-location of particle                  
            yy = sigma/2*sin(circle)+y(n,round(k)); %y-location of particle

            if variables == 3
                color = [abs(v(n,round(k)))/max(max(abs(v))) 0 0];
            else
                color = [abs(v(n,round(k)))/max(max(abs(v))) 0 abs(u(n,round(k)))/max(max(abs(u)))];
            end
            %color is based on the max velocity in v and u direction
            %red = high v, blue = high u, or combo

            fill(xx,yy,color,'erasemode','none');
        end        
        
        %add the text with the energy
        text(ax_x(1),ax_y(2),['V = ', num2str(V_sum(k)), '      T = ' ,num2str(T_sum(k))]);
        text((ax_x(1)+ax_x(2))/2,ax_y(2),['     Total = ',num2str(T_sum(k)+V_sum(k))]);

        set(gca,'DataAspectRatio',[1 1 1]); %set aspect ratio x:y to 1:1
        axis([ax_x ax_y]); %set the axis

        %make the movie
        F=getframe(fig);
        writeVideo(mov,F);

        %clear the figure for removing traces
        clf(fig);
    end

    %close the movie when done
    close(mov);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%