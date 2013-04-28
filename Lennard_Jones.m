function [fij] = Lennard_Jones(x1, x2, y1, y2)
% Computes the force on particle 1 due to LJ potentials.
% Uses reduced units, i.e. sigma=1, epsilon=1.
% Note that the forces calculated can be used for both atoms.

rij=sqrt((x2-x1)^2+(y2-y1)^2);                  % Computes the distance. Order doesn't actually matter since this is always a positive scalar
fij = 48*(rij^-14 -0.5*rij^-8)*[x1-x2 y1-y2];   % Computes the force (not potential) on atoms. This outputs a vector, positive if attracting, negative if repulsing.

end

