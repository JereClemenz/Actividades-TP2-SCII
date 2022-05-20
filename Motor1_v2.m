clear all;clc;close all;

% Parametros del motor

Laa=366e-6;J=5e-9;Ra=55.6;Bm=0;Ki=6.49e-3;Km=6.53e-3;

% [A,B,C,D] 

Mat_A=[-Ra/Laa -Km/Laa 0;Ki/J -Bm/J 0;0 1 0];
Mat_B=[1/Laa;0;0];
Mat_C=[0 0 1];

%Matriz Controlabilidad----------------------------------------------------
Mat_M=[Mat_B Mat_A*Mat_B Mat_A^2*Mat_B];%
%el calculo va hasta 2 pan n-1=2 y A tiene dimensiones de nxn
%Comprobacion si es controlable
rank(Mat_M);% el rango es 3 asi que si puede controlar

%Cálculo del LQR-----------------------------------------------------------
Q=diag([1/10 1/10000 40 ]);R=0.00001;
%Construcción del Hamiltoniano para el cálculo del controlador
Ha=[Mat_A -Mat_B*inv(R)*Mat_B'; -Q -Mat_A'];
[n,va]=size(Ha);%tamaño de Ha
[V,D]=eig(Ha);%returns diagonal matrix D of eigenvalues and matrix V whose columns are the corresponding right eigenvectors, so that A*V = V*D.
MX1X2=[];%extraigo autovectores negativos
for ii=1:n
    if real(D(ii,ii))<0
        MX1X2=[MX1X2 V(:,ii)];
    end
end
%divido la matriz
MX1=MX1X2(1:n/2,:);%M 
MX2=MX1X2(n/2+1:end,:);%PM 
% P=abs(MX2*inv(MX1));
P=real(MX2*inv(MX1));%PM*inv(M)=P
K=inv(R)*Mat_B'*P;
%Fin cálculo del controlador-----------------------------------------------
disp('Controlador ampliado en ')
eig(Mat_A-Mat_B*K)

%Cálculo del Observador---------------------------------------------------
% A_d=A';
% B_d=C';
% C_d=B';
% Qo=1e5*diag([1 1]);Ro=1e1;
% %Construcción del Hamiltoniano para el cálculo del observador
% Ho=[A_d -B_d*inv(Ro)*B_d'; -Qo -A_d'];
% [no,va]=size(Ho);
% [V,D]=eig(Ho);MX1X2=[];
% for ii=1:no
%  if real(D(ii,ii))<0
%  MX1X2=[MX1X2 V(:,ii)];
%  end
% end
% MX1=MX1X2(1:no/2,:); MX2=MX1X2(no/2+1:end,:);
% Po=real(MX2*inv(MX1));
% Ko=(inv(Ro)*B_d'*Po)';
% disp('Observador en ')
% eig(A-Ko*C)
%FIN cálculo del Observador------------------------------------------

%Calculo de ganancia de prealimentación de la referencia.
Gj=-inv(Mat_C*inv(Mat_A-Mat_B*K)*Mat_B);

%Verifico que todos los polos estén en el semiplano izquierdo
eig(Mat_A-Mat_B*K)


%Definciones para simular
h=1e-6;%preliminar
tiempo=round(2/h);%5 seg de simulacion
t=0:h:(tiempo*h-h);

%entrada de referencia que varia entre pi/2 a -pi/2 cada 0.3 seg
%array de entrada
ref=zeros(1,200000);
i=1;%empiieza a valer pi/2
ver=false;%sirve para asignarle 2 valores
contador = 300000;
while(i<(tiempo+1))
    
    if contador==0
       contador=300000;%0.3/h  
       ver=~ver;
    end
    
    if ver==true
        ref(1,i)=-pi/2;
        contador=contador-1;
    end
    
    if ver==false
        ref(1,i)=pi/2;
        contador=contador-1;
    end
    i=i+1;
end
% figure (1)
% plot(t,ref)


%entrada de torque se aplica a pi/2, mismo tiempos q la ref
TL_in=zeros(1,2000000);
i=1;
ver=false;%sirve para asignarle 2 valores
contador = 300000;
while(i<(tiempo+1))
    
    if contador==0
       contador=300000;%0.3/h  
       ver=~ver;
    end
    
    if ver==true
        TL_in(1,i)=0;
        contador=contador-1;
    end
    
    if ver==false
        TL_in(1,i)=1.15*10^-3;
        contador=contador-1;
    end
    i=i+1;
end
% figure
% plot(t,TL_in)

% %Condiciones iniciales
x1(1)=0; x2(1)=0; x3(1)=0; u(1)=0;
 
i=2;
while(i<(tiempo+1))
    
 estado=[x1(i-1); x2(i-1); x3(i-1)];
 u(i)=-K*estado+ref(i-1)*Gj;
 
 x1_p=-Ra*x1(i-1)/Laa-Km*x2(i-1)/Laa+u(i)/Laa;
 x2_p=Ki*x1(i-1)/J-Bm*x2(i-1)/J-TL_in(i-1)/J;
 x3_p=x2(i-1);
 

 x1(i)=x1(i-1)+h*x1_p;
 x2(i)=x2(i-1)+h*x2_p;
 x3(i)=x3(i-1)+h*x3_p;
 y_sal(i)=Mat_C*estado;
 i=i+1;
end
% t=0:h:tiempo*h-h;
% figure(1);hold on;
% plot(t,x3);grid on; title('ángulo');


figure (1)
hold on;grid on;
plot(t,ref);title('Referencia');
xlabel('Tiempo');ylabel('Radianes');
hold off;

figure (2)
hold on;grid on;
plot(t,TL_in);title('Perturbacion de torque');
xlabel('Tiempo');ylabel('Torque');
hold off;

figure (3)
hold on;grid on;
plot(t,ref);plot(t,x3,'r');title('Salida con referencia');
xlabel('Tiempo');ylabel('Torque');
hold off;

figure (4)
hold on;grid on;
plot(t,x1,'r');title('Corriente');
xlabel('Tiempo');ylabel('Amper');
hold off;

