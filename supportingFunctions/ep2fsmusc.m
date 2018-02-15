function fs = ep2fsmusc(ep,n)
xdot=ep(2); ydot=ep(4);
theta=ep(5);    thetadot=ep(6);
v=sqrt(xdot^2+ydot^2);                                      % To find v
% Perfectly correct; Dont be doubtful again.
if (ydot>0)
    delta=((-1)^n)*(atan(xdot/ydot)+theta);                 % To find delta
elseif (ydot<0)&&(xdot>=0)
    delta=((-1)^n)*(pi-atan(abs(xdot)/abs(ydot))+theta);    % To find delta
elseif (ydot<0)&&(xdot<0)
    delta=((-1)^n)*(-(pi-atan(abs(xdot)/abs(ydot)))+theta); % To find delta
end
fs(1)=v; fs(2)=delta;
fs(3)=theta; fs(4)=thetadot;