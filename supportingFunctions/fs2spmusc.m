function sp = fs2spmusc(fs,n)
global comXimp comYimp;
v=fs(1); delta=fs(2);
theta=fs(3); thetadot=fs(4);

xdot=((-1)^(n+1))*v*sin(delta+((-1)^n)*theta);
ydot=v*cos(delta+((-1)^n)*theta);

sp(1)=comXimp;
sp(2)=xdot;
sp(3)=comYimp;
sp(4)=ydot;
sp(5)=theta;
sp(6)=thetadot;