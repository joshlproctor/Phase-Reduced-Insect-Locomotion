function [frontxf frontyf midxf midyf rearxf rearyf] = footpP()
global gam1footfront gam1footmid gam1footrear;
global gam2footfront gam2footmid gam2footrear Vdes;
global Vchangefront Vchangemid Vchangerear;
datamusc;

frontgam1f = gam1footfront - 2*(Vdes-0.25)*(Vdes>0.25)*(Vdes<Vchangefront) - ...
    2*(Vchangefront-0.25)*(Vdes>=Vchangefront); % 1.75
midgam1f = gam1footmid - 1.5*(Vdes-0.25)*(Vdes>0.25)*(Vdes<Vchangemid) - ...
    1.5*(Vchangemid-0.25)*(Vdes>=Vchangemid); % 2.25
reargam1f = gam1footrear - 1.28*(Vdes-0.25)*(Vdes>0.25)*(Vdes<Vchangerear) - ...
    1.28*(Vchangerear-0.25)*(Vdes>=Vchangerear); % 2.5

frontgam2f = gam2footfront + 1.75*(Vdes-Vchangefront)*(Vdes>=Vchangefront); % 1.25, 1.75
midgam2f = gam2footmid + 2.5*(Vdes-Vchangemid)*(Vdes>=Vchangemid); % 2.5, 0.75
reargam2f = gam2footrear + 3*(Vdes-Vchangerear)*(Vdes>=Vchangerear); % 3, 0.64


frontqmag = sqrt( l1front^2 + l2front^2 - 2*l1front*l2front*cos(frontgam2f) );
midqmag = sqrt( l1mid^2 + l2mid^2 - 2*l1mid*l2mid*cos(midgam2f) );
rearqmag = sqrt( l1rear^2 + l2rear^2 - 2*l1rear*l2rear*cos(reargam2f) );

frontphi2 = acos( ( l1front^2 + frontqmag^2 - l2front^2 )/( 2*l1front*frontqmag ) );
midphi2 = acos( ( l1mid^2 + midqmag^2 - l2mid^2 )/( 2*l1mid*midqmag ) );
rearphi2 = acos( ( l1rear^2 + rearqmag^2 - l2rear^2 )/( 2*l1rear*rearqmag ) );

frontlenxf = -frontqmag*sin(frontgam1f+frontphi2);
frontlenyf = frontqmag*cos(frontgam1f+frontphi2);
midlenxf = midqmag*sin(midgam1f+midphi2);
midlenyf = midqmag*cos(midgam1f+midphi2);
rearlenxf = -rearqmag*sin(reargam1f+rearphi2);
rearlenyf = rearqmag*cos(reargam1f+rearphi2);

frontxf = -d1front + frontlenxf;
frontyf = d2front + frontlenyf;
midxf = d1mid + midlenxf;
midyf = d2mid + midlenyf;
rearxf = -d1rear + rearlenxf;
rearyf = d2rear + rearlenyf;

end