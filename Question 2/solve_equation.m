%Solve [P]'(t)=k3[ES]
%[ES]'(t)=k1[E][S]-k2[ES]-k3[ES]
%[E]'(t)=k2[ES]-k1[E][S]+K3[ES]
%[S]'(t)=k2[ES]-k1[E][S]
P0=0;
ES0=0;
E0=1; %1uM
S0=10; %10uM
h=0.000005; %min, time step
t=0:h:0.3; %t goes from 0 to 0.05 min
k1=100;
k2=600;
k3=150;
%initialize the concentration of P,E,ES,S at each time point.
Pstar=zeros(size(t));
ESstar=zeros(size(t));
Estar=zeros(size(t));
Sstar=zeros(size(t));
Pstar(1)=P0;
ESstar(1)=ES0;
Estar(1)=E0;
Sstar(1)=S0;
for i=1:(length(t)-1)
    Sp1=k3*ESstar(i);  %S means the slope of the line
    Ses1=k1*Estar(i)*Sstar(i)-k2*ESstar(i)-k3*ESstar(i);
    ES1=ESstar(i)+Ses1*h/2;
    Se1=k2*ESstar(i)-k1*Estar(i)*Sstar(i)+k3*ESstar(i);
    E1=Estar(i)+Se1*h/2;
    Ss1=k2*ESstar(i)-k1*Estar(i)*Sstar(i);
    S1=Sstar(i)+Ss1*h/2;
    
    
    Sp2=k3*ES1;
    Ses2=k1*E1*S1-k2*ES1-k3*ES1;
    ES2=ESstar(i)+Ses2*h/2;
    Se2=k2*ES1-k1*E1*S1+k3*ES1;
    E2=Estar(i)+Se2*h/2;
    Ss2=k2*ES1-k1*E1*S1;
    S2=Sstar(i)+Ss2*h/2;
    
    Sp3=k3*ES2;
    Ses3=k1*E2*S2-k2*ES2-k3*ES2;
    ES3=ESstar(i)+Ses3*h;
    Se3=k2*ES2-k1*E2*S2+k3*ES2;
    E3=Estar(i)+Se3*h;
    Ss3=k2*ES2-k1*E2*S2;
    S3=Sstar(i)+Ss3*h;
    
    Sp4=k3*ES3;
    Ses4=k1*E3*S3-k2*ES3-k3*ES3;
    Se4=k2*ES3-k1*E3*S3+k3*ES3;
    Ss4=k2*ES3-k1*E3*S3;
    
    
    Pstar(i+1)=Pstar(i)+(Sp1+2*Sp2+2*Sp3+Sp4)*h/6;
    ESstar(i+1)=ESstar(i)+(Ses1+2*Ses2+2*Ses3+Ses4)*h/6;
    Estar(i+1)=Estar(i)+(Se1+2*Se2+2*Se3+Se4)*h/6;
    Sstar(i+1)=Sstar(i)+(Ss1+2*Ss2+2*Ss3+Ss4)*h/6;
end
figure(1)
h1=plot(t,Pstar);
xlabel('t');
ylabel('[P]');
hold on;
figure(2)
h2=plot(t,ESstar);
xlabel('t');
ylabel('[ES]');
hold on;
figure(3)
h3=plot(t,Estar);
xlabel('t');
ylabel('[E]')
hold on;
figure(4)
h4=plot(t,Sstar);
xlabel('t');
ylabel('[S]');
hold on;
figure(5)
h5=plot(Sstar,k3*ESstar);
xlabel('[s]')
ylabel('v')
hold on;
[maxr,index]=max(k3*ESstar);
i=Sstar(index);
text(i,maxr,'o','color','r')
text(i,maxr,['(',num2str(i),',',num2str(maxr),')'],'color','k')
