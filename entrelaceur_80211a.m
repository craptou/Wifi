function bits_entrelaces = entrelaceur_80211a(bits_entree, Ncbps, Nbpsc)
%ENTRELACEUR_80211A  Entrelacement double permutation IEEE 802.11a (section 2.3)
%
% Entrées :
%   bits_entree : vecteur colonne de bits codés (longueur = n * Ncbps)
%   Ncbps       : nombre de bits codés par symbole OFDM (Nsd * Nbpsc)
%   Nbpsc       : nombre de bits par sous-porteuse (log2(M))
%
% Sortie :
%   bits_entrelaces : vecteur colonne de bits entrelacés (même longueur)

% Vérification que la longueur est un multiple de Ncbps
if mod(length(bits_entree), Ncbps) ~= 0
    error('La longueur de bits_entree doit être un multiple de Ncbps=%d', Ncbps);
end

% Nombre de symboles OFDM à entrelacer
nb_symboles = length(bits_entree) / Ncbps;

% Calcul du paramètre s selon la norme IEEE 802.11a
s = max(Nbpsc/2, 1);

% Pré-calcul des indices de permutation (indices MATLAB : base 1)
% k varie de 0 à Ncbps-1 (indices norme IEEE 802.11a)
k = (0:Ncbps-1)';

% Première permutation : k -> i (espacement fréquentiel)
% Objectif : séparer les bits adjacents sur différentes sous-porteuses
i = (Ncbps/16) * mod(k, 16) + floor(k/16);

% Deuxième permutation : i -> j (diversité dans le symbole)
% Objectif : éviter que 2 bits consécutifs soient sur les mêmes positions MSB/LSB
j = s * floor(i/s) + mod(i + Ncbps - floor(16*i/Ncbps), s);

% Table de permutation : entrée[k+1] -> sortie[j+1]
idx_sortie = j + 1;  % indices MATLAB base 1

% Application de l'entrelacement symbole par symbole OFDM
bits_entrelaces = zeros(size(bits_entree));

for n = 1:nb_symboles
    % Extraction du bloc de Ncbps bits du n-ième symbole OFDM
    debut = (n-1)*Ncbps + 1;
    fin = n*Ncbps;
    bloc_entree = bits_entree(debut:fin);
    
    % Permutation selon les indices calculés
    % Le bit à la position k (entrée) va à la position j (sortie)
    bloc_sortie = zeros(Ncbps, 1);
    for kk = 1:Ncbps
        bloc_sortie(idx_sortie(kk)) = bloc_entree(kk);
    end
    
    % Stockage du bloc entrelacé
    bits_entrelaces(debut:fin) = bloc_sortie;
end

end
