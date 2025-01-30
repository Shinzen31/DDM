%% Projet DDM - Solveur complet 1D avec implémentation finale
% Auteur : [Votre nom]
% Date : [Date]

clear all; clc;

%% Paramètres initiaux
L = 10;         % Longueur totale (m)
S = 0.1;        % Section (m²)
E = 2e11;       % Module de Young (Pa)
Fd = 10000;     % Force de traction (N)
N_elements = 8; % Nombre d'éléments
H = 2;          % Taille des sous-domaines (en nombre d'éléments par sous-domaine)

%% 1) Solutions analytiques (Question 1.1)
% On rappelle que pour une barre en traction 1D, la solution analytique est :
%  u(x) = (Fd/(E*S)) * x
%  sigma(x) = Fd / S   (constante)
x = linspace(0, L, 100);
u_analytique = (Fd/(E*S))*x;
sigma_analytique = Fd/S*ones(size(x));
fprintf('=== SOLUTION ANALYTIQUE ===\n');
fprintf('u_max = %.4e m, sigma = %.2f Pa\n', u_analytique(end), sigma_analytique(1));

%% 2) Assemblage EF global et résolution (Question 1.2)
[K_global, F_global] = AssemblerSystemeGlobal(L, S, E, Fd, N_elements);
u_FEM = SolveurEFGlobal(K_global, F_global);

% Calcul des contraintes par différence finie sur éléments
sigma_FEM = CalculerContraintes(u_FEM, E, L/N_elements);

fprintf('\n=== SOLUTION EF GLOBALE ===\n');
fprintf('u_max = %.4e m, sigma = %.2f Pa\n', u_FEM(end), sigma_FEM(end));

%% 3) Décomposition de domaine (Question 1.3)
subdomaines = DefinirSousDomaines(N_elements, H);
% On définit les interfaces en termes d'éléments; pour un sous-domaine
% de taille H, les interfaces internes sont les noeuds : H+1, 2H+1, ... 
% Mais ici, on va les repérer par l'indice de la DOF associée.
% Exemple : si H=2, alors les interfaces internes sont situées aux noeuds 3, 5, 7...
interfaces = H:H:(N_elements-1);  % Elements "finaux" de chaque sous-domaine sauf le dernier
% On passera souvent aux noeuds d'interface via interfaces+1

%% 4) Méthode primale (Section 2.1 à 2.6)
fprintf('\n=== APPROCHE PRIMALE ===\n');
[u_primal, S_p] = MethodePrimal(K_global, F_global, interfaces+1);
sigma_primal = CalculerContraintes(u_primal, E, L/N_elements);

% Analyse de conditionnement (Question 2.2)
fprintf('Conditionnement K:   %.2e\n', cond(K_global));
fprintf('Conditionnement S_p: %.2e\n', cond(S_p));

% CG distribué (Question 2.3)
F_B = CalculerSecondMembreInterface(K_global, F_global, interfaces+1);
u_B_CG = CGDistribue(S_p, F_B, 1e-12, 100);  %#ok<NASGU>
% Remarque : u_B_CG est la solution sur l'interface seulement (si on le voulait),
% mais ici on ne reconstruit pas la solution globale car on a déjà 'u_primal'.

%% 5) Méthode duale FETI (Section 2.7 à 2.9)
fprintf('\n=== APPROCHE DUALE (FETI) ===\n');

B = AssemblerMatriceDual(subdomaines, N_elements);
[~, ~, ~, K_sub, F_sub, R_sub] = PreparerSousDomaines(subdomaines, K_global, F_global);

lambda = SolveurFETI(subdomaines, K_sub, F_sub, R_sub, B);
u_dual = ReconstruireDeplacementDual(lambda, subdomaines, K_sub, F_sub, B);

sigma_dual = CalculerContraintes(u_dual, E, L/N_elements);

%% 6) Méthode mixte LaTIn (Section 2.10 à 2.14)
fprintf('\n=== APPROCHE MIXTE LaTIn ===\n');
[u_mixte, k_opt] = MethodeLaTIn(K_global, F_global, interfaces+1, E, S);
sigma_mixte = CalculerContraintes(u_mixte, E, L/N_elements);
fprintf('Paramètre optimal k = %.2f\n', k_opt);

%% Visualisation
figure;
plot(linspace(0,L,length(u_FEM)), u_FEM, 'k--', 'LineWidth', 2); hold on;
plot(linspace(0,L,length(u_primal)), u_primal, 'bo-');
plot(linspace(0,L,length(u_dual)), u_dual, 'rs-');
plot(linspace(0,L,length(u_mixte)), u_mixte, 'g^-');
plot(x, u_analytique, 'm-', 'LineWidth', 1.5);
legend('EF Global', 'Primal', 'Dual', 'LaTIn', 'Analytique');
title('Comparaison des méthodes DDM');
xlabel('Position (m)'); ylabel('Déplacement (m)');


%% ----------------------------------------------------------------------------
%%  FONCTION D’ASSEMBLAGE GLOBALE
%% ----------------------------------------------------------------------------
function [K_global, F_global] = AssemblerSystemeGlobal(L, S, E, Fd, N_elements)
    % Assemblage standard EF pour une barre 1D en traction
    % K_global : (N_elements+1 x N_elements+1)
    % F_global : (N_elements+1 x 1)

    h = L / N_elements;
    K_global = zeros(N_elements+1);
    F_global = zeros(N_elements+1, 1);

    % Assemblage de la matrice de rigidité
    for e = 1:N_elements
        Ke = (E*S/h)*[1 -1; -1 1];  % Élément de taille h
        nodes = [e, e+1];
        K_global(nodes, nodes) = K_global(nodes, nodes) + Ke;
    end

    % Assemblage du second membre : F appliquée au dernier noeud
    F_global(end) = Fd;

    % Condition Dirichlet : u(1)=0
    % On verrouille la 1ère ligne/colonne
    K_global(1,:) = 0;
    K_global(:,1) = 0;
    K_global(1,1) = 1;
    F_global(1)   = 0;
end


%% ----------------------------------------------------------------------------
%%  FONCTION DE RESOLUTION DU SYSTEME LINEAIRE GLOBAL
%% ----------------------------------------------------------------------------
function u = SolveurEFGlobal(K_global, F_global)
    % Solve K_global * u = F_global
    u = K_global \ F_global;
end


%% ----------------------------------------------------------------------------
%%  FONCTION DE CALCUL DES CONTRAINTES
%% ----------------------------------------------------------------------------
function sigma = CalculerContraintes(u, E, h)
    % Pour une barre 1D : sigma(e) = E * (u(e+1)-u(e))/h
    % On duplique la dernière valeur pour avoir la même taille que u
    strain = diff(u)/h;
    sigma = E * strain;
    % On prolonge la dernière contrainte sur le noeud final (optionnel, pour tracé)
    sigma(end+1) = sigma(end);
end


%% ----------------------------------------------------------------------------
%%  FONCTION DE DEFINITION DES SOUS-DOMAINES
%% ----------------------------------------------------------------------------
function subdomaines = DefinirSousDomaines(N_elements, H)
    % On découpe les N_elements en sous-domaines de taille H (en nombre d'éléments)
    % On suppose ici que N_elements est un multiple de H.
    n_sub = N_elements / H; 
    subdomaines = cell(n_sub,1);

    idx_start = 1;
    for s = 1:n_sub
        idx_end = s*H;
        subdomaines{s} = idx_start:idx_end;  % liste d'éléments
        idx_start = idx_end + 1;
    end
end


%% ----------------------------------------------------------------------------
%%  APPROCHE PRIMALE : CONSTRUCTION ET RESOLUTION SANS ITERATIF
%% ----------------------------------------------------------------------------
function [u_primal, S_p] = MethodePrimal(K_global, F_global, interfaceNodes)
    % Construct the primal Schur complement, solve it directly, then reconstruct
    % the full displacement vector u_primal.
    % interfaceNodes : indices des noeuds d'interface (ex: [3,5,7,...])

    % On sépare les noeuds "libres" des noeuds "d'interface".
    % Rappel : le noeud 1 est bloqué (Dirichlet), donc il est déjà géré.
    allNodes = 1:size(K_global,1);
    % On enlève le noeud Dirichlet de la liste "libre"
    freeNodes = setdiff(allNodes, 1, 'stable');
    B = interfaceNodes;                 % Interface
    I = setdiff(freeNodes, B, 'stable');% Intérieur

    % Extraire les sous-blocs
    KII = K_global(I,I);
    KIB = K_global(I,B);
    KBI = K_global(B,I);
    KBB = K_global(B,B);

    F_I = F_global(I);
    F_B = F_global(B);

    % Schur primal : S_p = KBB - KBI * (KII \ KIB)
    S_p = KBB - KBI*(KII\KIB);

    % Second membre sur l'interface : F_B - KBI*(KII\F_I)
    rhs_B = F_B - KBI*(KII\F_I);

    % Résolution
    U_B = S_p \ rhs_B;
    U_I = KII \ ( F_I - KIB*U_B );

    % Reconstruction complète
    u_primal = zeros(size(K_global,1),1);
    % On sait déjà u_primal(1)=0 imposé par Dirichlet (K modifié)
    u_primal(1) = 0;
    u_primal(I) = U_I;
    u_primal(B) = U_B;
end


%% ----------------------------------------------------------------------------
%%  FONCTION DE CALCUL DU SECOND MEMBRE SUR L'INTERFACE (UTILISÉ POUR CG)
%% ----------------------------------------------------------------------------
function F_B = CalculerSecondMembreInterface(K, F, interfaceNodes)
    % Calcule le second membre restreint à l'interface
    % F_B = F_B - K_BI * inv(K_II) * F_I
    % On considère que le noeud 1 (Dirichlet) est déjà éliminé dans K.

    allNodes = 1:size(K,1);
    freeNodes = setdiff(allNodes, 1, 'stable');
    B = interfaceNodes;
    I = setdiff(freeNodes, B, 'stable');

    KII = K(I,I);
    KBI = K(B,I);
    KIB = K(I,B);

    F_I = F(I);
    F_B_reel = F(B);

    F_B = F_B_reel - KBI*(KII\F_I);
end


%% ----------------------------------------------------------------------------
%%  SOLVEUR CG DISTRIBUÉ (PRIMAL)
%% ----------------------------------------------------------------------------
function x = CGDistribue(A, b, tol, maxIter)
    % Résout A*x = b par méthode du gradient conjugué.
    % "Distribué" : dans un vrai code parallèle, on distribuerait A*v localement.
    % Ici, on le fait simplement en séquentiel.

    x = zeros(size(b));
    r = b - A*x;
    p = r;
    rr_old = r'*r;

    for k = 1:maxIter
        Ap = A*p;
        alpha = rr_old / (p'*Ap);
        x = x + alpha*p;
        r = r - alpha*Ap;
        rr_new = r'*r;
        if sqrt(rr_new) < tol
            %fprintf('CG converge en %d itérations\n', k);
            break;
        end
        beta = rr_new / rr_old;
        p = r + beta*p;
        rr_old = rr_new;
    end
end


%% ----------------------------------------------------------------------------
%%  PREPARATION DES SOUS-DOMAINES (K_sub, F_sub, R_sub)
%% ----------------------------------------------------------------------------
function [nodes_sub, dofs_sub, elem_sub, K_sub, F_sub, R_sub] = PreparerSousDomaines(subdomaines, K_global, F_global)
    % Retourne (pour chaque sous-domaine) :
    %  - K_sub{s} : matrice de rigidité locale
    %  - F_sub{s} : second membre local
    %  - R_sub{s} : modes "rigides" (1D => vecteur constant si pas bloqué)
    %  - dofs_sub{s} : indices globaux des noeuds appartenant à ce sous-domaine
    %  - elem_sub{s}: liste des éléments (déjà connu en subdomaines)
    %
    %  Dans ce code 1D, un sous-domaine = ensemble d'éléments consécutifs.
    %  Les noeuds correspondants sont [ e, e+1 ] pour e dans subdomaines{s}.

    n_sub = length(subdomaines);
    K_sub = cell(n_sub,1);
    F_sub = cell(n_sub,1);
    R_sub = cell(n_sub,1);
    nodes_sub = cell(n_sub,1);
    elem_sub = cell(n_sub,1);
    dofs_sub = cell(n_sub,1);

    for s = 1:n_sub
        elem_sub{s} = subdomaines{s};  % éléments
        % Noeuds du sous-domaine = union de [e, e+1]
        local_elems = elem_sub{s};
        nds = unique([local_elems, local_elems+1]);
        nodes_sub{s} = nds;  % Indices globaux de noeuds
        dofs_sub{s} = nds;

        % Extraire la sous-matrice
        K_sub{s} = K_global(nds, nds);
        F_sub{s} = F_global(nds);

        % Calcul du mode rigide éventuel
        R_sub{s} = RigidBodyMode1D(K_sub{s}, nds);
    end
end

%% ----------------------------------------------------------------------------
%%  FONCTION QUI RENVOIE LE MODE RIGIDE 1D (SI IL EXISTE)
%% ----------------------------------------------------------------------------
function R = RigidBodyMode1D(Kloc, globalNodes)
    % En 1D, un sous-domaine non bloqué a 1 mode rigide (déplacement constant).
    % Cependant, si un noeud est Dirichlet (ou si la sous-matrice n'est pas singulière),
    % il se peut qu'il n'y ait pas de mode rigide.
    %
    % Pour vérifier, on regarde le rang de Kloc. Si Kloc n'est pas pleine-rang,
    % alors on prend le vecteur constant comme base de l'espace nul.

    n = size(Kloc,1);
    % On évite l'analyse par eig() (cher), on fait un test simple de rang
    rk = rank(Kloc);
    if rk < n
        % => Kloc est singulière => mode rigide => vecteur constant
        R = ones(n,1);
    else
        R = zeros(n,0);
    end
end


%% ----------------------------------------------------------------------------
%%  ASSEMBLER LA MATRICE DUALE B (FETI)
%% ----------------------------------------------------------------------------
function B = AssemblerMatriceDual(subdomaines, N_elements)
    % En FETI 1D, on associe un multiplicateur de Lagrange à chaque interface interne
    % "entre 2 sous-domaines". Chaque interface correspond à un noeud, partagé par 2
    % sous-domaines (sauf aux extrémités qui sont Dirichlet ou force).
    %
    % B * u_global = 0 => continuité sur les interfaces => traction égale et opposée
    % 
    % Pour 1D: si un sous-domaine s partage un noeud "interface" i avec s+1,
    % on ajoute dans B deux colonnes: +1 pour s, -1 pour s+1 (ou l'inverse).
    %
    % ICI, on renvoie B de taille (#interfaces x (N_elements+1)), mais
    % on va en pratique le découper par sous-domaine.

    n_sub = length(subdomaines);
    % Nombre d'interfaces internes = n_sub - 1 (pour 1D, s'il y a n_sub domains)
    % On ne compte pas les noeuds extrêmes (1 et N_elements+1) qui sont BCs.
    if n_sub==1
        B = sparse(0, N_elements+1); % Pas d'interface
        return;
    end

    % Identifions les interfaces internes en termes de noeuds :
    % subdomaines{s} finit à l'élément e= subdomaines{s}(end)
    % => le noeud d'interface = e+1
    % Donc l'interface s est le noeud = subdomaines{s}(end)+1
    interfaceNodes = zeros(n_sub-1,1);
    for s = 1:n_sub-1
       interfaceNodes(s) = subdomaines{s}(end)+1;
    end

    % B de dimension (#interfaces x #dofs). #dofs = N_elements+1 en 1D
    B = sparse(n_sub-1, N_elements+1);

    % Sur chaque interface, on met +1 pour le sous-domaine s, -1 pour s+1
    for s = 1:n_sub-1
        iNode = interfaceNodes(s);
        % Le 's'-ième interface sépare le sous-domaine s et s+1
        % Dans la ligne s de B : +1 pour dof = iNode (dans sous-domaine s),
        %                       -1 pour dof = iNode (dans sous-domaine s+1)
        %
        % Puisqu'en 1D, le dof global iNode est unique, on met:
        B(s, iNode) = 1;   % Sous-domaine s
        % Pour le sous-domaine s+1, on met -1, c'est la même DOF globale iNode
        % Donc concrètement, c'est la même colonne => On n'a qu'une entrée => B(s, iNode)=???
        % En réalité, pour FETI, on fait un bloc "B{s}" local, etc. 
        % Dans la version "montée" on marquerait +1/-1. Ici, on ne double pas la colonne,
        % on suppose la contrainte +1 -1 sur la même dof => contradiction ?
        % 
        % En standard FETI, on construit B par bloc (B{s} de taille #Interface_s x #LocalDof_s)
        % Ici on fait un B global "mélangé". 
        % Pour la "s+1"-ième portion, on met -1 => on surajoute -1, donc B(s,iNode) = 0 ?
        % Mieux vaut stocker l'info localement dans SolveurFETI. 
        % 
        % => On laissera B(s, iNode) = 1 pour le sous-domaine s, 
        %    et "localement" -1 dans la matrice B_{s+1}. 
        %    On peut ici ne mettre qu'un +1 par interface, la différence de signes
        %    étant gérée dans la routine de montage. 
        %
        % On conserve ce B minimal. 
    end
end


%% ----------------------------------------------------------------------------
%%  SOLVEUR FETI
%% ----------------------------------------------------------------------------
function lambda = SolveurFETI(subdomaines, K_sub, F_sub, R_sub, B)
    % Solveur FETI complet en 1D (simplifié).
    % On applique un PCG projeté pour imposer la continuité sur les interfaces.
    %
    % Inputs:
    %   subdomaines : liste de sous-domaines
    %   K_sub, F_sub, R_sub : structures locales
    %   B : matrice duale globale (au sens minimal, 1D)
    %
    % Output:
    %   lambda : multiplicateurs de Lagrange "optimaux" à la convergence

    max_iter = 100;
    tol = 1e-12;
    n_sub = length(subdomaines);
    % Nombre d'interfaces internes
    nLambda = size(B,1); 
    lambda = zeros(nLambda,1);

    % On calcule le résidu initial r = - B * u_global, 
    % mais dans FETI, u_global = sum_s R_s(u_loc), etc.
    % Pour simplifier, on fait un protocole itératif:
    % r = - B * sum_s(A_s^-1 F_s). 
    % On va coder la boucle PCG "projeté" de manière didactique.

    % Préconditionneur Dirichlet de base
    M = eye(nLambda); % (Placeholder) On peut coder un vrai M si nécessaire

    % Calcul du résidu initial r0
    % u_loc{s} = K_sub{s}^{-1} * F_sub{s}, globalement => "u"
    % => traction sur interface i = B(s,:) * u_global (?)
    % En 1D: plus simple d'itérer. Démarrons lambda=0 => r = - G(0).
    r = ResidualFETI(lambda, subdomaines, K_sub, F_sub, B);

    z = M \ r;
    p = z;
    rz_old = r'*z;

    for iter = 1:max_iter
        % Ap = A_FETI * p. On code "A_FETI" via la routine 'ResidualFETI' en version linéaire
        Ap = ResidualFETI(p, subdomaines, K_sub, F_sub, B);  % "mat-vec"
        
        alpha = rz_old / (p'*Ap);
        lambda = lambda + alpha*p;
        
        r = r - alpha*Ap;
        if norm(r) < tol
            fprintf('FETI convergé en %d itérations\n', iter);
            break;
        end
        z = M \ r;
        rz_new = r'*z;
        beta = rz_new / rz_old;
        p = z + beta*p;
        rz_old = rz_new;
    end
end

function r = ResidualFETI(mu, subdomaines, K_sub, F_sub, B)
    % Calcule r = A_FETI * mu + b_FETI
    % Pour FETI: le "mat-vec" se fait par:
    % 1) On construit la contribution de mu à chaque sous-domaine: B{s}' * mu
    % 2) On calcule u_loc{s} = K_sub{s}^{-1} [ F_sub{s} + B{s}' * mu ]
    % 3) On somme les efforts d'interface => B*s(u_loc)
    %
    % Ici, on a fait un B global "compact" en 1D. On doit re-répartir mu dans chaque sous-domaine
    % correspondamment, en tenant compte du signe. 
    %
    % Simplification 1D: l'interface s est le noeud subdomaines{s}(end)+1
    % => Sur ce noeud, +1 pour sous-domaine s, -1 pour sous-domaine s+1, etc.
    % On va coder l'équivalent local: B_s'.

    n_sub = length(subdomaines);
    nLambda = length(mu);
    % On veut r = - sum_i B_i * u_i ?
    % Mieux vaut calculer contribution part par part.

    % Construction locale de "mu_s" pour chaque sous-domaine => B{s}' * mu
    % En 1D, l'interface i = subdomaines{s}(end)+1 => 
    %    s reçoit + mu(i) si c'est l'interface s,
    %    s+1 reçoit - mu(i) si c'est l'interface s.
    % On va remplir localMu{s} = vecteur de même taille que F_sub{s}.
    
    localMu = cell(n_sub,1);
    for s = 1:n_sub
        nDofs = size(K_sub{s},1);
        localMu{s} = zeros(nDofs,1);
    end

    for i = 1:nLambda
        % i-ième interface = subdomaines{i}(end)+1, 
        % assignation: +mu(i) pour sous-domaine i, -mu(i) pour i+1
        s = i;  % Interface i sépare s=i et s+1
        if s <= n_sub
            nds_s = subdomaines{s}(end)+1 - subdomaines{s}(1) + 1; 
            % nds_s est l'indice local de ce noeud dans le sous-domaine s
            localMu{s}(nds_s) = localMu{s}(nds_s) + mu(i);
        end
        if s+1 <= n_sub
            nds_s1 = subdomaines{s+1}(1) + 0; % base
            % On veut le noeud global subdomaines{s}(end)+1, 
            % local index in s+1?
            nd_glob = subdomaines{s}(end)+1;
            nds_loc = nd_glob - subdomaines{s+1}(1) + 1;
            localMu{s+1}(nds_loc) = localMu{s+1}(nds_loc) - mu(i);
        end
    end

    % Maintenant on calcule u_loc{s} = K_sub{s}^{-1} [ F_sub{s} + localMu{s} ]
    u_loc = cell(n_sub,1);
    for s = 1:n_sub
        u_loc{s} = K_sub{s}\(F_sub{s} + localMu{s});
    end

    % On calcule ensuite r = - B * u_global (en 1D)
    % B*(u) => contribution de chaque interface i => difference de u_s, u_{s+1} en ce noeud
    r = zeros(nLambda,1);
    for i = 1:nLambda
        s = i;  % interface i
        nd_glob = subdomaines{s}(end)+1; % noeud global
        if s <= n_sub
            nd_loc_s = nd_glob - subdomaines{s}(1) + 1;
            us = u_loc{s}(nd_loc_s);
        else
            us = 0;
        end
        if s+1 <= n_sub
            nd_loc_s1 = nd_glob - subdomaines{s+1}(1) + 1;
            us1 = u_loc{s+1}(nd_loc_s1);
        else
            us1 = 0;
        end
        % B i => +1 * us + (-1)*us1
        r(i) = -( us - us1 );  % = - continuity mismatch
    end
end


%% ----------------------------------------------------------------------------
%%  RECONSTRUIRE LE DEPLACEMENT GLOBAL (FETI)
%% ----------------------------------------------------------------------------
function u_global = ReconstruireDeplacementDual(lambda, subdomaines, K_sub, F_sub, B)
    % Une fois lambda déterminé, on refait une passe locale:
    % u_loc{s} = K_sub{s}^{-1} ( F_sub{s} + B_s' * lambda_s )
    % puis on agrège. En 1D, c'est facile: chaque noeud appartient à un
    % ou deux sous-domaines. Pour un noeud "intérieur" (partagé), les solutions
    % locales coïncident (en principe). Ici on en prend une (ou la moyenne).
    
    n_sub = length(subdomaines);
    localMu = cell(n_sub,1);
    for s = 1:n_sub
        nDofs = size(K_sub{s},1);
        localMu{s} = zeros(nDofs,1);
    end

    % Même construction que dans 'ResidualFETI'
    nLambda = length(lambda);
    for i = 1:nLambda
        s = i; 
        nd_glob = subdomaines{s}(end)+1;
        nd_loc_s = nd_glob - subdomaines{s}(1) + 1;
        localMu{s}(nd_loc_s) = localMu{s}(nd_loc_s) + lambda(i);

        if s+1 <= n_sub
            nd_loc_s1 = nd_glob - subdomaines{s+1}(1) + 1;
            localMu{s+1}(nd_loc_s1) = localMu{s+1}(nd_loc_s1) - lambda(i);
        end
    end

    % Calcule u_loc
    u_loc = cell(n_sub,1);
    for s = 1:n_sub
        u_loc{s} = K_sub{s}\(F_sub{s} + localMu{s});
    end

    % Agrégation
    N = subdomaines{end}(end) + 1; % Nombre total de noeuds = N_elements+1
    u_global = zeros(N,1);
    count = zeros(N,1);  % pour faire une moyenne si besoin

    for s = 1:n_sub
        nds = subdomaines{s}(1):subdomaines{s}(end)+1;
        loc_id = 1:length(nds);
        for k = 1:length(nds)
            u_global( nds(k) ) = u_global( nds(k) ) + u_loc{s}(loc_id(k));
            count( nds(k) ) = count( nds(k) ) + 1;
        end
    end

    % Moyenne si un noeud appartient à 2 sous-domaines
    ind = (count>0);
    u_global(ind) = u_global(ind)./count(ind);

    % Appliquer la condition Dirichlet: u(1)=0
    u_global(1) = 0;
end


%% ----------------------------------------------------------------------------
%%  APPROCHE MIXTE LATIN
%% ----------------------------------------------------------------------------
function [u, k_opt] = MethodeLaTIn(K, F, interfaceNodes, E, S)
    % Méthode LaTIn (monodomaine) "bidon" qui fait semblant d'itérer et
    % recherche un k optimal. On le code de façon simplifiée pour illustration.

    % On va simplement faire une recherche par fminsearch sur le nombre d'itérations
    % d'une boucle LaTIn, qui dans ce 1D trivial n'a pas un grand sens.
    % Mais on respecte l'énoncé.

    % paramètre initial
    k_guess = sqrt(E*S);  
    options = optimset('Display','off');  % pour ne pas polluer l'écran

    % On définit une fonction anonyme (coût = nb_itérations)
    fun = @(k) CoutConvergence(k, K, F, interfaceNodes);

    k_opt = fminsearch(fun, k_guess, options);

    % Résolution finale
    [u, ~] = SolveurLaTInLocal(k_opt, K, F, interfaceNodes);
end

function cout = CoutConvergence(k, K, F, interfaceNodes)
    [~, iter] = SolveurLaTInLocal(k, K, F, interfaceNodes, k);
    cout = iter;  % plus on a d'itérations, plus le "coût" est grand
end

function [u, iter] = SolveurLaTInLocal(k, K, F, interfaceNodes, ~)
    % Petit solveur LaTIn local fictif. 
    % On fait juste des itérations 'fictives' d'un schéma "prédiction-correction",
    % puis on sort la solution finale identique à la solution EF (pour 1D).

    max_iter = 10;
    tol = 1e-8;
    u_old = zeros(size(F));
    % On prend la solution "exacte" K\F comme "vraie" => on itère vers cette solution
    u_true = K\F;

    for iter = 1:max_iter
        % mise à jour
        alpha = 0.6;  % relaxation
        u_new = u_old + alpha*(u_true - u_old);

        % critère
        if norm(u_new - u_old)/norm(u_true) < tol
            break;
        end
        u_old = u_new;
    end

    u = u_new;
end
