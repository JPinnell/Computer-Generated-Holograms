function [Holo,alpha,beta,gamma] = GenLGHolo3(P,L,weights,w0,Nx,Ny,shiftx,shifty,Arrizon,H,V)
%Generates the hologram for modulation of a flat beam to
%an LG mode; 
% P -> radial number, L -> topological charge, w0 -> beam waist, 
% weights -> mode weights (weights = 0 -> Phase only hologram but L must be a scalar)
% shiftx,y -> shift hologram, Nx,Ny are number of grating lines,  
% Arrizon refers to type of amplitude encoding in Arrizon (if 1 -> Type 1 Sinc, 0 Type 3 Bessel) 
% H,V are the column (x) and row (y) screen resolutions
% SNR best for Arrizon type 3 if w0 = 35*8e-6; 

%Not the most efficient implementation but this function highlights all the
%different modifications

% For SLM, H = 1920 (or 960 if want 2 holograms on 1 screen) and V = 1080
Ux = H*8e-6;
Uy = V*8e-6;
x=linspace(-Ux/2,Ux/2,H);
y=linspace(-Uy/2,Uy/2,V);
[X,Y]=meshgrid(x,y); X = X + shiftx; Y = Y + shifty;
% [theta,rho] = cart2pol(X,Y); %if needed

if weights == 0 %doing a phase-only hologram
    Phase = L.*atan2(Y,X);
    F = 1;
    gamma = sum(sum(abs(LG2(X,Y,P,L,1,w0).^2)));
    alpha = 1; beta = 1;
else %Doing 1 of 3 possible complex amplitude modulations
    mode = LG2(X,Y,P,L,weights,w0);
    Phase=angle(mode);
    A=abs(mode); alpha = max(max(A));
    A = A/alpha;
    beta = sum(sum(conj(A.*exp(1i*Phase)).*(A.*exp(1i*Phase))));
    gamma = sum(sum(abs(mode.^2)));
    
    load(strcat('fx',num2str(Arrizon),'.mat')); %load Arrizon
    aux = round(A.*(length(fx)-1)+1);
    for mh = 1:V
        for nh = 1:H
            temp = aux(mh,nh);
            F(mh,nh) = fx(temp);
        end
    end
end
Gx = Nx/(H*8e-6); Gy = Ny/(V*8e-6);

%Make hologram
if Arrizon == 1
    Holo = F.*mod(Phase+2*pi.*(Gx.*X+Gy.*Y),2*pi);
elseif Arrizon == 2
    Holo = Phase + F.*sin(Phase+2*pi.*(Gx.*X+Gy.*Y));
elseif Arrizon == 3
    Holo = F.*sin(Phase+2*pi.*(Gx.*X+Gy.*Y));
end  
Holo = Holo - min(Holo(:));
Holo = Holo/max(Holo(:))*255;

end