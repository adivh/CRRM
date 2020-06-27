% %function [aus]= bewertung%(art,vol,strike,n)
% n=10;   %Zeitschritte
% art='Put';
% S = 11911.35; %Aktueller Stand des Underlyings
% T=4/12; %Laufzeit von 3 Monaten
% vol = 0.3086;
% r = -0.04;                 %Mit festem u,d
% Opt=art;
% strike = 12000;
% 
% u = exp(vol*sqrt(T/n));
% d = 1/u;
% 
% mu=0;
% p=0.5+0.5*mu/vol*sqrt(T/n);
% 
% 
% 
% if(isequal(Opt,'Put'))
%    Ausz=zeros(1,n+1);
%    for ud = 1:n+1
%        Ausz(ud)=max(-S*u^(ud-1)*d^(n+1-ud)+strike,0);
%    end
%    for i=1:n
%        [Ausz]=bewer('Put',r,u,d,S,n-i+1,strike,Ausz)
%    end
% elseif(isequal(Opt,'Call'))
%     Ausz=zeros(1,n+1);
%     for ud = 1:n+1
%         Ausz(ud)=max(S*u^(ud-1)*d^(n+1-ud)-strike,0);
%     end
%     for i=1:n
%         [Ausz]=bewer('Call',r,u,d,S,n-i+1,strike,Ausz)
%     end
% else
%     "Bitte geben Sie als Optionsart 'Call' oder 'Put' ein"
% end
% 



%Black-Scholes
n=100;   %Zeitschritte
S = 11911.35; %Aktueller Stand des Underlyings
T=4/12; %Laufzeit von 3 Monaten
vol = 0.3086;
r = -0.04;                 %Mit festem u,d
Opt=art;
strike = 12000;

d1=(log(S/strike)+(r+vol^2/2)^T)/(vol*sqrt(T));
preis = S*normcdf(d1)-strike*(1+r)^(-T)*normcdf(d1-vol*sqrt(T));
fprintf("Der Preis nach dem Black-Scholes Modell ist: %d",preis);

