 clear all, close all,clc;

 % las matrices son
 % A=[-a a 0 0;0 0 1 0;w^2 -w^2 0 0;c 0 0 0]
 %B=[0 0 b*w^2 0]'
 %C=[1 0 0 0;0 0 0 1];D=0
 
% Ítem [3] Para el caso del avión, emplear un tiempo de integración por Euler adecuado y un tiempo 
% de simulación de 70seg. Los parámetros son a=0.05; =3; b=5; c=100, hallar un controlador para 
% que los polos de lazo cerrado se ubican en i=-1515j; -0.50.5j, para referencias de 100 y -100 
% metros en altura, ambas con alturas iniciales de -500 y 500.

w=3;a=0.05;b=5;c=100;

Dt=10^-3;
% ts=5;timpo de simulacion

tiempo=(70/Dt);%
t=0:Dt:(tiempo*Dt-Dt);%definimos el paso y el valor maximo de t
u=linspace(0,0,tiempo+1);

%reeferencias y Cond inciales
hInc=-500;hRef=100;

%Versión linealizada con el avion recto;
A=[-a, a, 0, 0;0, 0, 1, 0;w^2, -w^2, 0, 0;c, 0, 0, 0];
B=[0 0 b*w^2 0]';
C=[0 0 0 1];D=0;

%------Controlador--------------
Mat_M = [B A*B A^2*B A^3*B];       %matriz de controlabilidad 
rank(Mat_M)%chequea el rango de la matris M, tiene que ser 4

%---------ubicacion de los polos a lazo cerrado-------------
auto_val=eig(A);
c_ai=conv(conv(conv([1 -auto_val(1)],[1 -auto_val(2)]),[1 -auto_val(3)]),[1 -auto_val(4)]);
Mat_W=[;c_ai(4) c_ai(3) c_ai(2) 1;c_ai(3) c_ai(2) 1 0;c_ai(2) 1 0 0;1 0 0 0];
Mat_T=Mat_M*Mat_W;
A_controlable=inv(Mat_T)*A*Mat_T;%Verificación de que T esté bien
%Ubicación de los polos de lazo cerrado en mui:
mui(1)=-15+15*i;mui(2)=conj(mui(1));mui(3)=-0.5+0.5*i;mui(4)=conj(mui(3));
alfa_ia=conv(conv(conv([1 -mui(3)],[1 -mui(4)]),[1 -mui(2)]),[1 -mui(1)]);
K=(alfa_ia(2:5)-c_ai(2:5))*inv(Mat_T);
Gj=-inv(C*inv(A-B*K)*B);
eig(A-B*K) %Verifico que todos los polos estén en el semiplano izquierdo

%------condiciones iniciales-------
alfa(1)=0;
fi(1)=0;
fi_p(1)=0;
h(1)=hInc;
u(1)=0;

Xop=[0 0 0 0]';
x=[alfa(1) fi(1) fi_p(1) h(1)]';

ref(1)=hRef;

i=1;
while(i<(tiempo+1))
    ref(i)=hRef;%referencia fija
    estado=[alfa(i);fi(i);fi_p(i);h(i)];
    u(i)=-K*estado+Gj*ref(i);
    
    %--------sistema lineal--------
    xp=A*(x-Xop)+B*u(i);
    x=x+xp*Dt;
    alfa(i+1)=x(1);
    fi(i+1)  =x(2);
    fi_p(i+1)=x(3);
    h(i+1)   =x(4);

    i=i+1;
end

t=0:Dt:(tiempo*Dt);

figure(1)
subplot(3,1,1);
plot(t,alfa);hold on;title('\alpha_t');grid on;
subplot(3,1,2);
plot(t,fi);hold on;title('\phi_t');grid on;
subplot(3,1,3);
plot(t,h);hold on;grid on;title('altura (h)');
xlabel('Tiempo.[Seg]');


figure (2)
plot(t,u);hold on;title('u [Acción de control]');
grid on;xlabel('Tiempo.[Seg]');
