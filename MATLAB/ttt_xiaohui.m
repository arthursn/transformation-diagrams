clc;
clear;
close all;
% C=0.09;Mn=1.84;Si=0.36;Al=0.05;Mo=0.01;Cr=0.02;Cu=0.03;S=0.005;P=0.01;Ni=0;Co=0;
% Composition of SA508 Hamelin:
% C=0.2;Si=0.25;Mn=1.4;Ni=0.8;Cr=0;Mo=0.006;V=0.006;Al=0.027;Cu=0.031;Co=0;Ti=0.002;
% 
% % composition from Krajewski
C=0.093;Si=0.209;Mn=1.87;Ni=0.001;Cr=0;Mo=0.006;V=0.006;Al=0.027;Cu=0.031;Co=0;Ti=0.002;
% composition from CCT figure in Akbari
% C=0.2;Si=0.35;Mn=0.55;P=0.014;S=0.001;Cu=0;Cr=1.4;Ni=3.3;Al=0;
% Mo=0.2;V=0;Co=0;Ti=0;


% Calculation of FC, PC and BC. In Hamelin's paper, BC=FC, which could be
% an error. Here the calculation of BC is from the original work of Li. 
FC=exp((1.0+6.31*C+1.78*Mn+0.31*Si+1.12*Ni+2.7*Cr+4.06*Mo));
PC=exp(-4.25+4.12*C+4.36*Mn+0.44*Si+1.71*Ni+3.33*Cr+5.19*sqrt(Mo));
BC=exp(-10.23+10.18*C+0.85*Mn+0.55*Ni+0.9*Cr+0.36*Mo);

% Grange's calculation for Ae1 and Ae3. Ae3=1464, which is too high. 
Ae1_Grange=1333-25*Mn+40*Si-26*Ni-42*Cr;
Ae3_Grange=1570-323*C-25*Mn+80*Si-32*Ni-3*Cr;
% Andrews'es calculation for Ae1 and Ae3
Ae1_Andrews=712-17.8*Mn+20.1*Si-19.1*Ni+11.9*Cr+9.8*Mo;
Ae3_Andrews=871-254.4*sqrt(C)+51.7*Si-14.2*Ni;
Ae1 = 590;
% Ae1 = 715;
Ae3 = 970;
% Ae1 = 854;
%Ae1 and Ae3 based on CCT diagrams
% Ae1=715;Ae3=750;
% Bs calculation from Li's work
Bs=637-58*C-35*Mn-15*Ni-34*Cr-41*Mo;
% Bs based on CCT diagram
% Bs=600;
Ms=539-423*C-30.4*Mn-17.7*Ni-12.1*Cr-7.5*Mo+10*Co- 7.5*Si;
G=6;R=8.314459;
para=struct('Ae1',Ae1,'Ae3',Ae3,'Bs',Bs,'Ms',Ms,'G',G,'FC',FC,'BC',BC,'PC',PC);
%% plot TTT diagram
% plot ferrite part 
T1=[Bs:3:Ae3];
for i=1:1:length(T1)
    T=T1(i);
    tf1(i)=FC/(2^(0.41*G)*(Ae3-T)^3*exp(-27500*4.184/(R*(T+273.15))))*0.1043;
end 
semilogx(tf1,T1,'LineWidth',1.5);
hold on;

% plot pearlite part 
tf1 =0;
T1=[Bs:3:Ae3];
for i=1:1:length(T1)
    T=T1(i);
    tf1(i)=PC/(2^(0.32*G)*(Ae1-T)^3*exp(-27500*4.184/(R*(T+273.15))))*0.1043;
end 
semilogx(tf1,T1,'LineWidth',1.5);

% plot bainite part 
tf1 =0;
T1=[Ms:3:Bs];
for i=1:1:length(T1)
    T=T1(i);
    tf1(i)=BC/(2^(0.29*G)*(Bs-T)^2*exp(-27500*4.184/(R*(T+273.15))))*0.1043;
end 
semilogx(tf1,T1,'LineWidth',1.5);


xlimit=[0.01 1e8];
semilogx(xlimit,[Ms Ms],'LineWidth',1.5);
xlim(xlimit);ylim([300,1000]);
xlabel('Time (s)');ylabel('Temperature ^\circC');

%% plot CCT figure
T1=1000;
for j=1:1:Ae3-1-Ms
    T2=Ae3-1-j;
    dtime=12000;inc=1e-2;
    for i=1:1:10
        error(i)=clcCCT(T1,T2,dtime,Ae3,Ms,0.32,G,3,FC);
        if (error(i) < 0.001)
            break
        end 
        error1=clcCCT(T1,T2,dtime+inc,Ae3,Ms,0.32,G,3,FC);
        slope=(error1-error(i))/inc;
        dtime=dtime-error(i)/slope;
    end 
    Tf(j)=T2;
    stimef(j)=dtime;
end 

for j=1:1:Ae1-1-Ms
    T2=Ae1-1-j;
    dtime=12000;inc=1e-2;
    for i=1:1:50
        error(i)=clcCCT(T1,T2,dtime,Ae1,Ms,0.32,G,3,PC);
        if (error(i) < 0.001)
            break
        end 
        error1=clcCCT(T1,T2,dtime+inc,Ae1,Ms,0.32,G,3,PC);
        slope=(error1-error(i))/inc;
        dtime=dtime-error(i)/slope;
    end 
    Tp(j)=T2;
    stimep(j)=dtime;
end 

for j=1:1:Bs-1-Ms
    T2=Bs-1-j;
    dtime=12000;inc=1e-2;
    for i=1:1:10
        error(i)=clcCCT(T1,T2,dtime,Bs,Ms,0.29,G,2,BC);
        if (error(i) < 0.001)
            break
        end 
        error1=clcCCT(T1,T2,dtime+inc,Bs,Ms,0.29,G,2,BC);
        slope=(error1-error(i))/inc;
        dtime=dtime-error(i)/slope;
    end 
    Tb(j)=T2;
    stimeb(j)=dtime;
end 

for index1=length(stimef):-1:1
    timeb=interp1(Tb,stimeb,Tf(index1),'linear','extrap');
    if stimef(index1)<timeb
%         index1 = index1+1;
        break
    end 
end 
for index2=length(stimep):-1:1
    timeb=interp1(Tb,stimeb,Tp(index2),'linear','extrap');
    if stimep(index2)<timeb
%         index2 = index2+1;
        break
    end 
end 
figure();
semilogx(stimef(1:index1),Tf(1:index1),'LineWidth',1.5);hold on;
semilogx(stimep(1:index2),Tp(1:index2),'LineWidth',1.5);
semilogx(stimeb,Tb,'LineWidth',1.5);
xlim([0.01 1e8]);
semilogx([0.01 10000],[Ms Ms],'LineWidth',1.5);
ylim([300,1000]);
xlabel('Time (s)');ylabel('Temperature ^\circC');

%% importdata and plot temperature line in CCT figure
data=importdata('abaqus.rpt',' ',4);
data1=data.data;
[maxT,maxTi]=max(data1(:,2));
data=data1(maxTi:end,:);
index=data(:,2)<=T1;
t0=interp1(data(:,2),data(:,1),T1);
data=[t0,T1;data(data(:,2)<=T1,:)];
data(:,1)=data(:,1)-t0;data(1,1)=data(1,1)+1e-8;
semilogx(data(:,1),data(:,2),'LineWidth',1.5);

%% calcualte the phase evolution curves according to the temperature evolution curve

data2(:,2)=[T1:-1:200];
data2(:,1)=interp1(data(:,2),data(:,1),data2(:,2));

austenite0=struct('f',1,'incub',0,'v',0,'t0',0,'T',0,'t',0);
pearlite0=struct('f',0,'incub',0,'v',0,'t0',0,'T',0,'t',0);
ferrite0=struct('f',0,'incub',0,'v',0,'t0',0,'T',0,'t',0);
bainite0=struct('f',0,'incub',0,'v',0,'t0',0,'T',0,'t',0);
martensite0=struct('f',0,'incub',0,'v',0,'t0',0,'T',0,'t',0);
material0=[austenite0,ferrite0,pearlite0,bainite0,martensite0];
% material=material0;

for i=2:1:length(data2(:,1))
    dt=data2(i,1)-data2(i-1,1);
    T1=data2(i-1,2);
    T2=data2(i,2);
    T=T1+(T2-T1)/2;
%     dT=1;n=abs(T2-T1)/dT;dT=(T2-T1)/n;
    % Ae3=para.Ae3;Ae1=para.Ae1;Bs=para.Bs;Ms=para.Ms;
    % FC=para.FC;PC=para.PC;BC=para.BC;G=para.G;
    if i==2
        ferrite = material0(2);
        pearlite = material0(3);
        bainite = material0(4);
        austenite = material0(1);
    end 
    
    if T<Ae3 && T>=Ms && bainite.incub<1 
        ferrite=clcphase(ferrite,austenite,T,Ae3,dt,0.41,G,3,FC);
       
    end 

    if T<Ae1 && T>=Ms && bainite.incub<1     
        pearlite=clcphase(pearlite,austenite,T,Ae1,dt,0.32,G,3,PC);
    end 

    if T<Bs && T>=Ms
        bainite=clcphase(bainite,austenite,T,Bs,dt,0.29,G,2,BC);
    end 
    austenite.f=1-ferrite.f-bainite.f;
    allf(i,:)=[austenite.f,ferrite.f,pearlite.f,bainite.f];
    allincub(i,:)=[austenite.incub,ferrite.incub,pearlite.incub,bainite.incub];
end 
allf(1,:)=allf(2,:);
allincub(1,:)=allincub(2,:);
figure();
for i=1:1:4
    hold on;
    plot(data2(:,1),allf(:,i),'LineWidth',1.5);
end 
xlim([0 20]);

% function to update phase fraction, which is appplicable for
% ferrite,bainite and pearlite
function phase2=clcphase(phase,austenite,T,Ts,dt,n2,G,nT,Fe)
R=8.314459;
t0=Fe/(2^(n2*G)*(Ts-T)^nT*exp(-27500*4.184/(R*(T+273.15))))*0.1043;
fac=phase.f+austenite.f;
if phase.incub<=1
    t0=Fe/(2^(n2*G)*(Ts-T)^nT*exp(-27500*4.184/(R*(T+273.15))))*0.1043;
    phase.incub = phase.incub + dt/t0;
    phase.t0=t0;
    phase.T=T;
else 
    phase.f=max(phase.f,0.01*fac);
    x=phase.f/fac;
    phase.v=(2^(n2*G)*(Ts-T)^nT*exp(-27500*4.184/(R*(T+273.15))))...
        *x^(0.4*(1-x))*(1-x)^(0.4*x)/Fe*fac;
    phase.f=phase.f+min(phase.v*dt,austenite.f);
    phase.T=T;phase.t0=t0;
end 
phase2=phase;
end 

%% function to calcualte the integrated incubation value
function error=clcCCT(T1,T2,dtime,Ts,Te,n2,G,nT,Fe)
    n=abs(T2-T1)/1;
    dT=(T2-T1)/n;
    incub=0;dt=dtime/n;
    for i=1:1:n
        T=T1+dT*(i-1)+dT/2;
        if T<Ts && T>=Te
            t0=Fe/(2^(n2*G)*(Ts-T)^nT*exp(-27500*4.184/(8.314*(T+273.15))))*0.1043;
            incub = incub + dt/t0;
        end 
    end
    error=incub-1;
end 