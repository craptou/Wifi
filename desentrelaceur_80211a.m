function bits_desentrelaces = desentrelaceur_80211a(bits_entree, Ncbps, Nbpsc)
%DESENTRELACEUR_80211A  Désentrelacement inverse IEEE 802.11a (section 2.3)
%
% Entrées :
%   bits_entree : vecteur colonne de bits entrelacés (longueur = n * Ncbps)
%   Ncbps       : nombre de bits codés par symbole OFDM (Nsd * Nbpsc)
%   Nbpsc       : nombre de bits par sous-porteuse (log2(M))
%
% Sortie :
%   bits_desentrelaces : vecteur colonne de bits dans l'ordre original

% Vérification que la longueur est un multiple de Ncbps
if mod(length(bits_entree), Ncbps) ~= 0
    error('La longueur de bits_entree doit être un multiple de Ncbps=%d', Ncbps);
end

% Nombre de symboles OFDM à désentrelacer
nb_symboles = length(bits_entree) / Ncbps;

% Calcul du paramètre s selon la norme IEEE 802.11a
s = max(Nbpsc/2, 1);

% Pré-calcul des indices de DÉSENTRELACEMENT (inverse exact de l'entrelacement)
% On reconstruit la table k -> j de l'entrelaceur, puis on l'inverse

% Recalcul des permutations directes (identique à l'entrelaceur)
k = (0:Ncbps-1)';

% Première permutation : k -> i
i = (Ncbps/16) * mod(k, 16) + floor(k/16);

% Deuxième permutation : i -> j
j = s * floor(i/s) + mod(i + Ncbps - floor(16*i/Ncbps), s);

% Construction de la table inverse : pour chaque position j en entrée du désentrelaceur,
% on doit retrouver la position k d'origine
idx_inverse = zeros(Ncbps, 1);
for kk = 1:Ncbps
    jj = j(kk) + 1;  % position de sortie de l'entrelaceur (base 1)
    idx_inverse(jj) = kk;  % cette position jj vient de la position kk
end

% Application du désentrelacement symbole par symbole OFDM
bits_desentrelaces = zeros(size(bits_entree));

for n = 1:nb_symboles
    % Extraction du bloc de Ncbps bits du n-ième symbole OFDM
    debut = (n-1)*Ncbps + 1;
    fin = n*Ncbps;
    bloc_entree = bits_entree(debut:fin);
    
    % Permutation inverse : le bit à la position j (entrée) revient à la position k (sortie)
    bloc_sortie = zeros(Ncbps, 1);
    for jj = 1:Ncbps
        kk = idx_inverse(jj);
        bloc_sortie(kk) = bloc_entree(jj);
    end
    
    % Stockage du bloc désentrelacé
    bits_desentrelaces(debut:fin) = bloc_sortie;
end

end
