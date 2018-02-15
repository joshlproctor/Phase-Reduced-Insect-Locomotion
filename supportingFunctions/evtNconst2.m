function [value,isterminal,dircn] = evtNconst2(t,x,n,tendlast,TDindex)

%global Ffrontleft Fmidright Frearleft Ffrontright Fmidleft Frearright
global Raghu_MNPhases


Phase1=mod(x(1,1),1);
Phase2=mod(x(26,1),1);
Phase_MN=mod(x(2:25,1),1);

value(1)=x(1,1)-.55;   isterminal(1)=1;   dircn(1)=+1;  %%% Indicates Left Touchdown
value(2)=x(26,1)-.05;   isterminal(2)=1;   dircn(2)=+1;  %%% Indicates Left Touchdown
value(3)=x(1,1)-1; isterminal(3)=1; dircn(3)=+1;  %%% Indicates Left Touchdown


for j=1:24

value(j+3)=x(j+1,1)-1; isterminal(j+3)=1; dircn(j+3)=1;  %%% Indicates Left Touchdown

end

value(28)=x(26)-1; isterminal(28)=1; dircn(28)=1;




end

