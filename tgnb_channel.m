function [h, H, Heff] = tgnb_channel(P, seed)
%TGNB_CHANNEL  Génère une réalisation stationnaire du canal TGn-B.
% Renvoie h, H=fft(h,N) et Heff (avec epsilon pour éviter /0).

if nargin < 2 || isempty(seed)
    seed = P.seedChannel;
end

rng(seed);

delay_samp = round((P.tau_ns*1e-9)/P.Ts);
Lh = max(delay_samp) + 1;

h = zeros(Lh,1);
for kk = 1:length(P.tau_ns)
    beta = sqrt(P.p_lin(kk)/2) * (randn + 1j*randn);
    h(delay_samp(kk)+1) = h(delay_samp(kk)+1) + beta;
end

H = fft(h, P.N);

epsH = 1e-12;
Heff = H;
Heff(abs(Heff) < epsH) = epsH;

end
