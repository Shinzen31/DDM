%% Projet DDM - Solveur complet 1D avec FETI refondu (FETI-DP style)
% Auteur : [Votre nom]
% Date : [Date]

clear all; clc;

%% Paramètres initiaux
L = 10;         % Longueur totale (m)
S = 0.1;        % Section (m²)
E = 2e11;       % Module de Young (Pa)
Fd = 10000;     % Force de traction (N)
N_elements = 16; % Nombre d'éléments
H = 4;          % Taille des sous-domaines (en nb d'éléments par sous-domaine)
                % => ici 4 sous-domaines, 3 interfaces

%% 1) Solutions analytiques (Question 1.1)
x = linspace(0, L, 100);
u_analytique = (Fd/(E*S)) * x;
sigma_analytique = (Fd/S) * ones(size(x));
fprintf('=== SOLUTION ANALYTIQUE ===\n');
fprintf('u_max = %.4e m, sigma = %.2f Pa\n', u_analytique(end), sigma_analytique(1));

%% 2) Assemblage EF global et résolution (Question 1.2)
[K_global, F_global] = AssemblerSystemeGlobal(L, S, E, Fd, N_elements);
u_EF = SolveurEFGlobal(K_global, F_global);
sigma_EF = CalculerContraintes(u_EF, E, L/N_elements);
fprintf('\n=== SOLUTION EF GLOBALE ===\n');
fprintf('u_max(EF) = %.4e m, sigma = %.2f Pa\n', u_EF(end), sigma_EF(end));

%% 3) Décomposition de domaine (Question 1.3)
subdomaines = DefinirSousDomaines(N_elements, H);
fprintf('\nSous-domaines:\n');
for s = 1:length(subdomaines)
    fprintf('  SD #%d = [%s]\n', s, num2str(subdomaines{s}));
end

%% 4) Méthode primale (Schur)
fprintf('\n=== APPROCHE PRIMALE ===\n');
interface_elems = H : H : (N_elements-1);
interfaceNodes = interface_elems + 1;
[u_primal,S_p] = MethodePrimal(K_global, F_global, interfaceNodes);
fprintf('u_max(Primal) = %.4e m\n', u_primal(end));

%% 5) Méthode duale FETI (refondu, type FETI-DP)
fprintf('\n=== APPROCHE DUALE (FETI) ===\n');
[~, ~, ~, K_sub, F_sub, ~] = PreparerSousDomaines(subdomaines, K_global, F_global);
tol_feti = 1e-8;  % tolérance dans l'espace de déplacement
maxit_feti = 200;
lambdaFETI = SolveurFETI_Dual_fixed_2(subdomaines, K_sub, F_sub, tol_feti, maxit_feti);
u_dual = ReconstruireDeplacementDual(lambdaFETI, subdomaines, K_sub, F_sub, u_primal);
sigma_dual = CalculerContraintes(u_dual, E, L/N_elements);
fprintf('u_max(FETI) = %.4e m\n', u_dual(end));

%% 6) Méthode mixte LaTIn (simplifiée)
fprintf('\n=== APPROCHE MIXTE LaTIn ===\n');
[u_mixte, k_opt] = MethodeLaTIn(K_global, F_global, interfaceNodes, E, S);
sigma_mixte = CalculerContraintes(u_mixte, E, L/N_elements);
fprintf('Paramètre optimal k = %.2f\n', k_opt);
fprintf('u_max(LaTIn) = %.4e m\n', u_mixte(end));

%% Visualisation
figure('Color','w'); hold on; grid on;
plot(linspace(0, L, length(u_EF)), u_EF, 'k--', 'LineWidth',2);
plot(linspace(0, L, length(u_primal)), u_primal, 'bo-');
plot(linspace(0, L, length(u_dual)), u_dual, 'rs-');
plot(linspace(0, L, length(u_mixte)), u_mixte, 'g^-');
plot(x, u_analytique, 'm-', 'LineWidth',1.5);
legend('EF Global','Primal','Dual (FETI)','LaTIn','Analytique','Location','best');
xlabel('Position (m)'); ylabel('Déplacement (m)');
title('Comparaison des méthodes DDM en 1D');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fonctions locales
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Assemblage EF global (barre 1D)
function [K_global,F_global] = AssemblerSystemeGlobal(L,S,E,Fd,N_el)
    h = L/N_el;
    ndofs = N_el+1;
    K_global = zeros(ndofs);
    F_global = zeros(ndofs,1);
    for e=1:N_el
        Ke = (E*S/h)*[1 -1; -1 1];
        nodes = [e, e+1];
        K_global(nodes,nodes) = K_global(nodes,nodes) + Ke;
    end
    F_global(end) = Fd;
    % Condition Dirichlet: u(1)=0
    K_global(1,:) = 0;
    K_global(:,1) = 0;
    K_global(1,1) = 1;
    F_global(1) = 0;
end

%% Solveur EF global
function u = SolveurEFGlobal(K,F)
    u = K\F;
end

%% Calcul des contraintes
function sigma = CalculerContraintes(u,E,h)
    strain = diff(u)/h;
    sigma = E*strain;
    sigma(end+1) = sigma(end);
end

%% Définition des sous-domaines
function subdom = DefinirSousDomaines(N_el,H)
    n_sub = N_el/H;
    subdom = cell(n_sub,1);
    idx = 1;
    for s=1:n_sub
        subdom{s} = idx:(idx+H-1);
        idx = idx+H;
    end
end

%% Méthode Primal (Schur)
function [u_primal,S_p] = MethodePrimal(K,F,interfaceNodes)
    ndofs = size(K,1);
    allNodes = 1:ndofs;
    freeNodes = setdiff(allNodes,1,'stable');
    B = interfaceNodes;
    I = setdiff(freeNodes,B,'stable');
    KII = K(I,I);
    KBB = K(B,B);
    KIB = K(I,B);
    KBI = K(B,I);
    FI = F(I);
    FB = F(B);
    S_p = KBB - KBI*(KII\KIB);
    rhsB = FB - KBI*(KII\FI);
    UB = S_p\rhsB;
    UI = KII\(FI - KIB*UB);
    u_primal = zeros(ndofs,1);
    u_primal(1) = 0;
    u_primal(I) = UI;
    u_primal(B) = UB;
end

%% Préparation des sous-domaines
function [nodes_sub,dofs_sub,elem_sub,K_sub,F_sub,R_sub] = PreparerSousDomaines(subdom,K_global,F_global)
    n_sub = length(subdom);
    nodes_sub = cell(n_sub,1);
    dofs_sub = cell(n_sub,1);
    elem_sub = cell(n_sub,1);
    K_sub = cell(n_sub,1);
    F_sub = cell(n_sub,1);
    R_sub = cell(n_sub,1);
    for s=1:n_sub
        elem_sub{s} = subdom{s};
        local_elems = elem_sub{s};
        nds = unique([local_elems, local_elems+1]);
        nodes_sub{s} = nds;
        dofs_sub{s} = nds;
        K_sub{s} = K_global(nds,nds);
        F_sub{s} = F_global(nds);
        R_sub{s} = RigidBodyMode1D(K_sub{s});
    end
end

function R = RigidBodyMode1D(Kloc)
    n = size(Kloc,1);
    if rank(Kloc) < n
        R = ones(n,1);
    else
        R = zeros(n,0);
    end
end

%% Assemblage de la matrice duale B minimale (1D)
function B = AssemblerMatriceDual(subdom, N_el)
    n_sub = length(subdom);
    if n_sub <= 1
        B = sparse(0, N_el+1);
        return;
    end
    n_int = n_sub - 1;
    B = sparse(n_int, N_el+1);
    for i = 1:n_int
        interf_node = subdom{i}(end)+1;
        B(i, interf_node) = +1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FETI – Nouvelle méthode duale (FETI-DP style)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lambda = SolveurFETI_Dual_fixed_2(subdom, K_sub, F_sub, tol, maxit)
    n_sub = length(subdom);
    nInt = n_sub - 1;
    if nInt < 1
        lambda = [];
        return;
    end
    lambda = zeros(nInt, 1);
    omega = 0.5;
    eps_reg = 1e-6;  
    for it = 1:maxit
        % Initialisation des contributions locales de lambda
        localMu = cell(n_sub,1);
        for s = 1:n_sub
            nDofs = size(K_sub{s},1);
            localMu{s} = zeros(nDofs,1);
        end
        % Répartition de lambda sur les noeuds interfaces
        for i = 1:nInt
            interfaceNode = subdom{i}(end) + 1; % noeud global interface
            nodes_i = subdom{i}(1):subdom{i}(end)+1;
            nodes_ip1 = subdom{i+1}(1):subdom{i+1}(end)+1;
            locIndex_i = find(nodes_i == interfaceNode);
            locIndex_ip1 = find(nodes_ip1 == interfaceNode);
            localMu{i}(locIndex_i) = localMu{i}(locIndex_i) + lambda(i);
            localMu{i+1}(locIndex_ip1) = localMu{i+1}(locIndex_ip1) - lambda(i);
        end
        % Résolution des problèmes locaux (sans correction de décalage)
        u_loc = cell(n_sub,1);
        for s = 1:n_sub
            Ktmp = K_sub{s};
            if rank(Ktmp) < size(Ktmp,1)
                Ktmp = Ktmp + eps_reg*eye(size(Ktmp,1));
            end
            u_loc{s} = Ktmp \ (F_sub{s} + localMu{s});
        end
        % Calcul du saut (jump) aux interfaces
        J = zeros(nInt,1);
        for i = 1:nInt
            interfaceNode = subdom{i}(end) + 1;
            nodes_i = subdom{i}(1):subdom{i}(end)+1;
            nodes_ip1 = subdom{i+1}(1):subdom{i+1}(end)+1;
            locIndex_i = find(nodes_i == interfaceNode);
            locIndex_ip1 = find(nodes_ip1 == interfaceNode);
            J(i) = u_loc{i}(locIndex_i) - u_loc{i+1}(locIndex_ip1);
        end
        normJ = norm(J);
        fprintf(' it=%3d | ||J||=%e\n', it, normJ);
        if normJ < tol
            fprintf('FETI converge en %d itérations.\n', it);
            break;
        end
        lambda = lambda + omega * J;
    end
end

%% Reconstruction du déplacement global pour FETI
function u_global = ReconstruireDeplacementDual(lambda, subdom, K_sub, F_sub, u_primal)
    n_sub = length(subdom);
    N = subdom{end}(end) + 1;
    u_global = zeros(N,1);
    count = zeros(N,1);
    eps_reg = 1e-6;
    % Recalcul des contributions locales de lambda
    localMu = cell(n_sub,1);
    for s = 1:n_sub
        nDofs = size(K_sub{s},1);
        localMu{s} = zeros(nDofs,1);
    end
    nInt = n_sub - 1;
    for i = 1:nInt
        interfaceNode = subdom{i}(end) + 1;
        nodes_i = subdom{i}(1):subdom{i}(end)+1;
        nodes_ip1 = subdom{i+1}(1):subdom{i+1}(end)+1;
        locIndex_i = find(nodes_i == interfaceNode);
        locIndex_ip1 = find(nodes_ip1 == interfaceNode);
        localMu{i}(locIndex_i) = localMu{i}(locIndex_i) + lambda(i);
        localMu{i+1}(locIndex_ip1) = localMu{i+1}(locIndex_ip1) - lambda(i);
    end
    % Résolution locale et correction par décalage à l'aide de u_primal pour s>=2
    for s = 1:n_sub
        Ktmp = K_sub{s};
        if rank(Ktmp) < size(Ktmp,1)
            Ktmp = Ktmp + eps_reg*eye(size(Ktmp,1));
        end
        sol = Ktmp \ (F_sub{s} + localMu{s});
        if s > 1
            % Pour s>=2, on impose que le premier noeud (interface gauche) corresponde à u_primal
            interfaceNode = subdom{s}(1);
            nodes = subdom{s}(1):subdom{s}(end)+1;
            locIndex = find(nodes == interfaceNode);
            offset = u_primal(interfaceNode) - sol(locIndex);
            sol = sol + offset;
        end
        nodes = subdom{s}(1):subdom{s}(end)+1;
        for k = 1:length(nodes)
            u_global(nodes(k)) = u_global(nodes(k)) + sol(k);
            count(nodes(k)) = count(nodes(k)) + 1;
        end
    end
    u_global(count>0) = u_global(count>0) ./ count(count>0);
    u_global(1) = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Méthode LaTIn (simplifiée)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u,k_opt] = MethodeLaTIn(K, F, interfaceNodes, E, S)
    k_guess = sqrt(E*S);
    options = optimset('Display','off');
    fun = @(k) CoutConvergence(k,K,F,interfaceNodes);
    k_opt = fminsearch(fun,k_guess,options);
    [u, ~] = SolveurLaTInLocal(k_opt, K, F, interfaceNodes);
end

function c = CoutConvergence(k,K,F,interfaceNodes)
    [~, iter] = SolveurLaTInLocal(k,K,F,interfaceNodes);
    c = iter;
end

function [u, iter] = SolveurLaTInLocal(k, K, F, interfaceNodes)
    max_iter = 20;
    tol = 1e-10;
    u_exact = K\F;
    u_old = zeros(size(F));
    for iter = 1:max_iter
        alpha = 0.6;
        u_new = u_old + alpha*(u_exact - u_old);
        if norm(u_new - u_old)/norm(u_exact) < tol
            break;
        end
        u_old = u_new;
    end
    u = u_new;
end
