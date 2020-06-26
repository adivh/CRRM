function [bew] = bewer(Opt,r,u,d,S,n,strike,Ausz)


Wert=zeros(1,n+1);
bew=zeros(1,n);
Wert2=zeros(1,n);
for ud = 1:n+1
        
        Wert(ud)=S*u^(ud-1)*d^(n+1-ud);
        
end
for ud=1:n
    Wert2(ud)=S*u^(ud-1)*d^(n-ud);
end

if(isequal(Opt,'Call'))
    
    for i = 1:n
        delta= (Ausz(i+1)-Ausz(i))/(Wert(i+1)-Wert(i));
        bew(i)=-(Wert(i+1)*delta-Ausz(i+1))/(1+r)+delta*Wert2(i);
        %Bei amerikanischem Call
%         if(Wert2(i)-strike>bew(i))
%             bew(i)= strike-Wert2(i);
%         end
    end


elseif(isequal(Opt,'Put'))
     for i = 1:n
        delta= (Ausz(i+1)-Ausz(i))/(Wert(i+1)-Wert(i));
        bew(i)=-(Wert(i+1)*delta-Ausz(i+1))/(1+r)+delta*Wert2(i);
        %Bei amerikanischem Put
%         if(strike-Wert2(i)>bew(i))
%             bew(i)=Wert2(i)-strike;
%         end
    end
end