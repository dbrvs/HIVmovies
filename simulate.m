%% make movies of SIV trajectories from Perelson 1996 model
% DBR May 2016
% inputs:
%   length of simulation in days
%   movieflag for 2D=2 or 3D=3
% outputs a .gif movie to the current directory

function HIVmovies(tf,movieflag)

tic;
%model parameters
beta = 7.7e-3;   %infectivity [uL/virion/day]
dS   = 0.045;    %susceptible death rate [1/day] Luo big range 0.045->0.45
aS   = 10*dS;    %constant growth rate of T cells 35-760 [cells/uL/day], I converted using body 6L to get [cells/day]
dI   = 1;
g    = 3;        %clearance rate Luo 18.8 [1/day]
p    = g*1e3;    %[viral expansion rate Luo 2.4-9.8 [virions/cell]
ep   = 0.9;      %ART efficacy

bt = beta*(1-ep); %calculate infectivity on ART

eqSIV = [dI*g/bt/p; aS/dI-dS*g/bt/p; aS*p/dI/g-dS/bt]; %calculate viral set point equilibria

options = odeset('RelTol', 1e-10); %tolerances for ODE solver

SIV0=[aS/dS 0 1e-3]; %initial conditions for S, I and V

tPts=1e4; %number of timepoints for simulation
[td, pop]=ode23s(@SIVode,linspace(0,tf,tPts),SIV0,options,[bt aS dI dS g p]);
  
%make movies dependent on what kind we want 2D or 3D
  if movieflag==2
    X=pop(:,1); Y=pop(:,2); %define the populations to plot
    fn='hiv2D.gif'; %filename
    nframes=200; %number of frames, must be divisor of tPts
    SIVmovie2D(td,X,Y,nframes,fn)

  elseif movieflag==3
    X=pop(:,1); Y=pop(:,2); Z=pop(:,3); %define the populations to plot
    fn='hiv3D.gif'; %filename
    SIVmovie3D(td,X,Y,Z,eqSIV,fn)   
  end

end

%% make 2D movies
%inputs:
%  vectors of t,X(t),Y(t),Z(t)
%  filename fn (example.gif)
function SIVmovie2D(t,X,Y,nframes,fn)

close all
fig2D = figure('Position', [100, 100, 500, 500]); %make figure 5x5
  hold on
  xlabel('S (cells per \muL)'); ylabel('I (cells per \muL)') %axes labels
  xlim([0,max(X)]); ylim([0,max(Y)]); %keep limits constant
  %loop over time points making nframes plots
  for i=1:length(t)/nframes:length(t);         
    scatter(X(i),Y(i),'k','filled'); %scatter S,I data
    legend(['day ' num2str(t(i),2)])
    %pause(.01); 

    %make the gif using matlabs built in imwrite function
    F=getframe(fig2D);
      im = frame2im(F);
      [imind,cm] = rgb2ind(im,256);
      if i == 1;
          imwrite(imind,cm,fn,'gif','Loopcount',inf);
      else
          imwrite(imind,cm,fn,'gif','WriteMode','append');
      end
  end %loop over time

  hold off
end %function SIVmovie2D

%% make 3D movie of a trajectory
%inputs:
%  vectors of t,X(t),Y(t),Z(t)
%  equilibrium solution vector eqSIV=[S*,I*,V*]
%  filename fn (example.gif)

function SIVmovie3D(t,X,Y,Z,eqSIV,fn)

close all
fig3D = figure('Position', [100, 100, 500, 500]); %make figure 5x5
hold on
  xlim([0,max(X)]); ylim([0,max(Y)]); zlim([50,1e8]) %keep limits constant
  set(gca,'zscale','log')
  view(-15,30) %change angle of movie
  grid on; box on; %change style of movie
  %loop over time points making 300 plots
  for i=1:length(t)/250:length(t);         
  
    %put in loop for several trajectories at once!!
    
    scatter3(X(i),Y(i),Z(i)*1e3,'k','filled'); 
    title(['day ' num2str(t(i),2)])
    scatter3(eqSIV(1),eqSIV(2),eqSIV(3)*1e3,100,'r'); %put in equilibrium solution

    F=getframe(fig3D);
      im = frame2im(F);
      [imind,cm] = rgb2ind(im,256);
      if i == 1;
          imwrite(imind,cm,fn,'gif', 'Loopcount',inf);
      else
          imwrite(imind,cm,fn,'gif','WriteMode','append');
      end
  end %loop over time
  
  hold off

end

% function that calculates the differential equation using the rates above
% note nothing is time dependent
function dY=SIVode(~,X,pm)
  %pm stands for for parameter in list
  bt=pm(1); aS=pm(2); dI=pm(3); dS=pm(4); g=pm(5); p=pm(6);
  S=X(1); I=X(2); V=X(3); %rename for readability
  dY=zeros(3,1); %initialize vector to make columnar
  
  %ode solutions for SIV model
  dY(1)= aS - dS*S - bt*S*V;
  dY(2)= bt*S*V - dI*I;
  dY(3)= p*I - g*V;
end
