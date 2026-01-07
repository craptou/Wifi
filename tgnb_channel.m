function [h, H, Heff, delais_ech] = tgnb_channel(P, seed)
%TGNB_CHANNEL  Génère une réalisation stationnaire du canal TGn-B.
%
% Sorties :
%   h         : réponse impulsionnelle discrète du canal
%   H         : réponse fréquentielle H[k] sur Nfft points (FFT de h)
%   Heff      : H[k] "sécurisé" (évite division par 0 pour ZF)
%   delais_ech: retards en nombre d'échantillons

if nargin < 2 || isempty(seed)
    seed = P.seed_canal;
end
rng(seed);

% Conversion retards (ns) -> retards en échantillons
delais_ech = round((P.tau_ns * 1e-9) / P.Te);
Lh = max(delais_ech) + 1;

% Réponse impulsionnelle : somme de trajets sur les mêmes indices de retard
h = zeros(Lh,1);
for k = 1:length(P.tau_ns)
    % coefficient complexe gaussien : CN(0, p_lin(k))
    beta = sqrt(P.p_lin(k)/2) * (randn + 1j*randn);
    h(delais_ech(k)+1) = h(delais_ech(k)+1) + beta;
end

% Réponse fréquentielle sur Nfft
H = fft(h, P.Nfft);

% Protection contre les valeurs très petites (division en ZF)
epsH = 1e-12;
Heff = H;
Heff(abs(Heff) < epsH) = epsH;

end
