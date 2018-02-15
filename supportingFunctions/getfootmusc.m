function [xfootfr,yfootfr,xfootm,yfootm,xfoothind,yfoothind] = getfootmusc(fourset,x,y,n)
% global gaussfoot turn;


%% Final decision
theta=fourset(3);
Rtheta=[cos(theta) -sin(theta); sin(theta) cos(theta)];

% if gaussfoot==0
    [frontxf frontyf midxf midyf rearxf rearyf]=footpP;
    


% elseif gaussfoot==1
%     [frontxf frontyf midxf midyf rearxf rearyf]=footpP;
%     
%     frontxf=dev_frontx +frontxf;    frontyf=dev_fronty + frontyf;
%     midxf=dev_midx + midxf;         midyf=dev_midy + midyf;
%     rearxf=dev_rearx + rearxf;      rearyf=dev_reary + rearyf;
% end

%% Turning ((n==6)|(n==8)|(n==10))
% if ((turn==1)&&(mod(n,2)==0)&&(n>5)&&(n<15))
%     frontxf=frontxf+0.003;
%     frontyf=frontyf+0.003;
%     rearxf=rearxf-0.002;
%     rearyf=rearyf+0.002;
%     n
% end



%% Final assignment
front=Rtheta*[((-1)^n)*frontxf; frontyf];
    xfootfr=x+front(1,1); yfootfr=y+front(2,1);
mid=Rtheta*[((-1)^n)*midxf; midyf];
    xfootm=x+mid(1,1); yfootm=y+mid(2,1);
rear=Rtheta*[((-1)^n)*rearxf; rearyf];
    xfoothind=x+rear(1,1); yfoothind=y+rear(2,1);
    