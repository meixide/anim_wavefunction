% Representamos a evolucion temporal de psi en 0 < x < L
clear;
s=2; % sigma
x0=20; % posicion inicial
m=1; % masa
hb=1; % cte de planck con barra
L = 100; % lonxitude do intervalo
N = 400; % numero de puntos
x = linspace(0,L,N)'; % coordenadas
	dx = x(2) - x(1); % paso
	k0 = 2; 
	dk = 2*pi/L; % paso do espazo de momentos
	km=N*dk; % borde do espazo de momentos
	k=linspace(0,km,N)'; % espazo de momentos

NT = 200; % numero de instantes
TF = 29; T = linspace(0,TF,NT); % linha do tempo
%inicializamos os vectores que conteran a evolucion temporal dos parametros
media=T;
dispersion2=T; 
mediap=T;
dispersion2p=T;

for t = 1:NT % para cada instante
% integrando no instante T
% salvo eikx
phi=exp(1i*(-(k-k0)*x0-((hb*(k.^2)*T(t))/(2*m)))-((s^2)*(k-k0).^2)/2);
psi = ifft(phi); % multiplicamos por eikx e integramos en todo o espazo de momentos 
%(Inverse Fast Fourier Transform)

psi = psi/sqrt(psi'*psi*dx); % normalize the psi(x,0)
	rho = conj(psi).*psi; % amplitude de probabilidade (modulo ao cadrado)
	plot(x,rho,'k'); 
	axis([0 L 0 0.3]); 
	xlabel('x [m]', 'FontSize', 10);
	ylabel('densidade de prob. [1/m]','FontSize', 10);
	pause(0.05);
	media(t)=trapz(x,x.*rho); %aproveitamos para calcular a media
	dispersion2(t)=trapz(x,(x.^2).*rho)-media(t)^2; %e o momento de orde 2-media^2 = disp^2
	%en todo instante
	
	phi = phi/sqrt(phi'*phi*dk); %normalizamos a funcion de onda en representacion p
rhop = conj(phi).*phi; %atopamos a densidade de probabilidade
mediap(t)=trapz(k,k.*rhop); %aproveitamos para calcular a media
dispersion2p(t)=trapz(k,(k.^2).*rhop)-mediap(t)^2;
%e o momento de orde 2-media^2 = disp^2

end

plot(T,media,'k'); 
xlabel('t [s]', 'FontSize', 20);
ylabel('<x>','FontSize', 20);
saveas(gcf,'mediax.png')

plot(T,dispersion2,'k'); 
xlabel('t [s]', 'FontSize', 20);
ylabel('(deltax)^2','FontSize', 20);
saveas(gcf,'dispx2.png')

plot(T,mediap,'k'); 
xlabel('t [s]', 'FontSize', 20);
ylabel('<k>','FontSize', 20);
saveas(gcf,'mediap.png')
plot(T,dispersion2p,'k'); 
xlabel('t [s]', 'FontSize', 20);
ylabel('(deltak)^2','FontSize', 20);
saveas(gcf,'dispk2.png')

