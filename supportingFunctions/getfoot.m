function getfoot(sixset,x,y,ind)
global xffrontleft yffrontleft xfmidright yfmidright xfrearleft yfrearleft;
global xffrontright yffrontright xfmidleft yfmidleft xfrearright yfrearright;

theta=sixset(5);
Rtheta=[cos(theta) -sin(theta); sin(theta) cos(theta)];

[frontxf frontyf midxf midyf rearxf rearyf]=footpP;
switch ind
    case 1
        n=0;
        front=Rtheta*[((-1)^n)*frontxf; frontyf];
        xffrontleft=x+front(1,1); yffrontleft=y+front(2,1);    
    case 2
        n=0;
        mid=Rtheta*[((-1)^n)*midxf; midyf];
        xfmidright=x+mid(1,1); yfmidright=y+mid(2,1);
    case 3
        n=0;
        rear=Rtheta*[((-1)^n)*rearxf; rearyf];
        xfrearleft=x+rear(1,1); yfrearleft=y+rear(2,1);
    case 4
        n=1;
        front=Rtheta*[((-1)^n)*frontxf; frontyf];
        xffrontright=x+front(1,1); yffrontright=y+front(2,1);
    case 5
        n=1;
        mid=Rtheta*[((-1)^n)*midxf; midyf];
        xfmidleft=x+mid(1,1); yfmidleft=y+mid(2,1);
    case 6
        n=1;
        rear=Rtheta*[((-1)^n)*rearxf; rearyf];
        xfrearright=x+rear(1,1); yfrearright=y+rear(2,1);
end

end