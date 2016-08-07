
%%%  TLM 2D serie method MATLab implementation for a waveguide
%%%  Author:  Vinícius Pletsch

%%

clear all
clc 
close all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Características do meio
co=3e+8; 
Eps0=8.85*1E-12;
Mu0=4*pi*1E-7;  
velocidade=1/sqrt(Eps0*Mu0);
Z0=120*pi;               % Impedancia do espaço livre
Ztl=Z0/sqrt(2);

%% Dimensões física do guia de onga e modos de propagação
a=19.05E-3;
b=9.525E-3; 

c=52;                  % Número de pontos no eixo x  (+ 2 colunas)
l=27;                  % Número de pontos no eixo y  (+ 2 linhas) 

%%%% Dimensões maiores para não zerar indíces de i-1=0


deltal = 0.381e-3;                % Distancia entre os nós)
deltat=0.898e-12;

%% Modo de propagação 
m=1;             
n=0;

fc=(co/(2*pi))*sqrt((((m*pi)/a)^2)+(((n*pi)/b)^2));   %Frequência de corte
lambdac=co/fc;
cj=0.0+j*1.0;
T=(1/fc);        %período de oscilação de cada modo

Tk=round((T/2)/deltat);         %T/2 correspondente ao número de interações



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Inicialização dos vetores e matrizes 


ninter=3000;           % Número de interações
ttotal=ninter*deltat;  % Tempo total de simulação
pg(ninter)=0;          % Vetor Pulso Gaussiano/Interações 

mvi(l,c,4,ninter)=0;   % Matriz de Vis todos pontos repetida 4x (4portas) para cada interação
mvr(l,c,4,ninter)=0;   % Matriz de Vrs todos pontos repetida 4x (4portas) para cada interação
Ex(l,c,ninter)=0;      % Matriz de Ex para todos os pontos para cada interação
Ey(l,c,ninter)=0;      % Matriz de Ey para todos os pontos para cada interação
Hz(l,c,ninter)=0;      % Matriz de Hz para todos os pontos para cada interação
Hzs(ninter)=0;         % Vetor com valores de Hz (ponto específico) para cada interação
Exs(ninter)=0;         % Vetor com valores de Ex (ponto específico) para cada interação
Eys(ninter)=0;         % Vetor com valores de Ey (ponto específico) para cada interação
fEy(l,c,ninter)=0;
fHz(l,c,ninter)=0;
fftHzs(ninter)=0;
fftEys(ninter)=0;

%%  Formatação do Pulso Gausiano

Eo=1;             %Amplitude do pulso 
D=10*deltat;      % Duração do pulso 
L=10*deltat;      % Largura do pulso 


% pg=Eo*(exp(-18*((t-D)/L)^2));

%% Ínico da interação no tempo (loop principal)

for k=1:ninter
    t=deltat*k;
    pg(k)=((Eo*(exp(-18*(((t-D)/L)^2))))*((0.5*deltal))); % Valores de V2=V4


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Excitação forçada (Sobrescreve os valores de tensão pulso enquanto ele existir)   
    if k<=19
   
    mvi(2:26,26,2,k)= - pg(k);
    mvi(2:26,26,4,k)= - pg(k);
    mvi(2:26,26,1,k)= 0;
    mvi(2:26,26,3,k)= 0;
   
    
    end
    
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% Início de inteações de linhas e colunas 
  for i=2:l-1
     for j=2:c-1
 
 
%%%% Cálculo dos campos
 Ex(i,j,k)=-((mvi(i,j,1,k)+mvi(i,j,3,k))/deltal);
 Ey(i,j,k)=-((mvi(i,j,2,k)+mvi(i,j,4,k))/deltal);
 Hz(i,j,k)=((mvi(i,j,1,k)+mvi(i,j,4,k)-mvi(i,j,3,k)-mvi(i,j,2,k))/(2*Ztl))/deltal; 
 
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% DTF para todos os pontos da malha para um frequência de corte 
 if k>1
    fEy(i,j,k)= fEy(i,j,k-1) + Ey(i,j,k).*exp(-2*cj*pi*fc*t);
    fHz(i,j,k)= fHz(i,j,k-1) + Hz(i,j,k).*exp(-2*cj*pi*fc*t);
 else
    fEy(i,j,k)=  Ey(i,j,k).*exp(-2*cj*pi*fc*t);
    fHz(i,j,k)=  Hz(i,j,k).*exp(-2*cj*pi*fc*t);
 
 end
 

%% Matriz Vs refletidas no mesmo K 
 
 mvr(i,j,1,k)= 0.5*(mvi(i,j,1,k)+ mvi(i,j,2,k)+ mvi(i,j,3,k)- mvi(i,j,4,k));
 mvr(i,j,2,k)= 0.5*(mvi(i,j,1,k)+ mvi(i,j,2,k)- mvi(i,j,3,k)+ mvi(i,j,4,k));
 mvr(i,j,3,k)= 0.5*(mvi(i,j,1,k)- mvi(i,j,2,k)+ mvi(i,j,3,k)+ mvi(i,j,4,k));
 mvr(i,j,4,k)= 0.5*(-mvi(i,j,1,k)+ mvi(i,j,2,k)+ mvi(i,j,3,k)+ mvi(i,j,4,k));

 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %% Conexão com nó seguinte 
        
             mvi(i+1,j,3,k+1)=mvr(i,j,1,k);
             mvi(i-1,j,1,k+1)=mvr(i,j,3,k);
             mvi(i,j-1,4,k+1)=mvr(i,j,2,k);    
             mvi(i,j+1,2,k+1)=mvr(i,j,4,k);

         

 
     end
  end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Condição de contorno 
             mvi(26,:,1,k+1)=-mvr(26,:,1,k);
             mvi(2,:,3,k+1)=-mvr(2,:,3,k);
             mvi(:,2,2,k+1)=-mvr(:,2,2,k);   
             mvi(:,51,4,k+1)=-mvr(:,51,4,k); 
  
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 % Saídas para um ponto específico da malha 
 Hzs(k)=Hz(13,29,k);
 Exs(k)=Ex(13,29,k);
 Eys(k)=Ey(13,29,k);   
 

end






vetorf=linspace(0,30,3000);



vetort=(0:1:ninter-1)*deltat;
% vetorf=[0:S3:30E9-S3];

% figure(1)
% plot(pg)

figure(2)
plot(vetort,Hzs);
grid
legend('Hz')

figure(3)
plot(vetort,Eys);
grid
legend('Ey')

% figure(4)
% plot(vetorf,abs(fft(Hzs,3000)));
% grid
% 
% figure(5)
% plot(vetorf,abs(fft(Eys,3000)));
% grid


figure(6)
surf(abs(fEy(2:26,2:51,1500)));
legend('fEy')


figure(7)
surf(abs(fHz(2:26,2:51,1500)));
legend('fHz');

figure(8)
quiver(fHz(2:26,2:51,71),fEy(2:26,2:51,71))
axis([0 51 0 26])
xlabel('Nx')
ylabel('Ny')




