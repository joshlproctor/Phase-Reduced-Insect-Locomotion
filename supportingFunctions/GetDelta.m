function [Delta_FBi,Delta_FBe]= GetDelta(kneeTorque,hipTorque,TDindex)

Delta_FBi=zeros(24,1);
Delta_FBe=zeros(24,1);

feedbackparam

for j=1:6
   if TDindex(j)==1
      
       HipFBi=ChangeSlope(j)*SlopeHipi(j)*abs(hipTorque(j)-Hip(j))+yintHipi(j);
       if hipTorque(j)>Hip(j)
           Delta_FBi(1+4*(j-1))=(HipFBi>1)*1+(HipFBi>0)*(HipFBi<1)*HipFBi;
           
       else 
           Delta_FBi(2+4*(j-1))=(HipFBi>1)*1+(HipFBi>0)*(HipFBi<1)*HipFBi;
           
       end
       
       KneeFBi=ChangeSlope(j)*SlopeKneei(j)*abs(kneeTorque(j)-Knee(j))+yintKneei(j);
       if kneeTorque(j)<Knee(j)
           Delta_FBi(3+4*(j-1))=(KneeFBi>1)*1+(KneeFBi>0)*(KneeFBi<1)*KneeFBi;
       else 
           Delta_FBi(4+4*(j-1))=(KneeFBi>1)*1+(KneeFBi>0)*(KneeFBi<1)*KneeFBi;
       end
       
       HipFBe=ChangeSlope(j)*SlopeHipe(j)*abs(hipTorque(j)-Hip(j))+yintHipe(j);
       if hipTorque(j)<Hip(j)
           Delta_FBe(1+4*(j-1))=(HipFBe>1)*1+(HipFBe>0)*(HipFBe<1)*HipFBe;
       else 
           Delta_FBe(2+4*(j-1))=(HipFBe>1)*1+(HipFBe>0)*(HipFBe<1)*HipFBe;
       end
       
       KneeFBe=ChangeSlope(j)*SlopeKneee(j)*abs(kneeTorque(j)-Knee(j))+yintKneee(j);
       if kneeTorque(j)>Knee(j)
           Delta_FBe(3+4*(j-1))=(KneeFBe>1)*1+(KneeFBe>0)*(KneeFBe<1)*KneeFBe;
       else 
           Delta_FBe(4+4*(j-1))=(KneeFBe>1)*1+(KneeFBe>0)*(KneeFBe<1)*KneeFBe;
       end
          
   end
    
end

end

