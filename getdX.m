function dX = getdX(X, U)
th = X(2,1);
dz = X(3,1);
dth = X(4,1);

global M m L J Dc Dth g

Mmat = [M+m     m*L*cos(th); m*L*cos(th)    J+m*L^2];
F = [ m*L*dth^2*sin(th)-Dc*dz+U; m*g*L*sin(th)-Dth*dth];
ddX = Mmat\F;

dX = [ dz; dth; ddX ];



