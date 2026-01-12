function [h, H, Heff, delais_ech] = tgnb_channel(P, seed)
% Génère une réalisation du canal TGn-B

rng(seed);

delais_ech = round((P.tau_ns * 1e-9) / P.Te);
Lh = max(delais_ech) + 1;

h = zeros(Lh, 1);
for k = 1:length(P.tau_ns)
    beta = sqrt(P.p_lin(k)/2) * (randn + 1j*randn);
    h(delais_ech(k)+1) = h(delais_ech(k)+1) + beta;
end

H = fft(h, P.Nfft);
Heff = H;

end
