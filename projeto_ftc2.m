function projeto_ftc_gauss()
    % =========================================================================
    % SIMULAÇÃO TÉRMICA 2D - CÉLULA LFP (ADI) COM CURRENT CROWDING (GAUSSIANA)
    % =========================================================================
    clc; clear; close all;

    % --- Parâmetros Termofísicos e Geométricos ---
    Lx = 0.10; Ly = 0.15; rho = 2240; cp = 800; k = 3.0;
    alpha = k / (rho * cp); T_inf = 25.0; h = 64.0;

    % --- Discretização da Malha ---
    Nx = 40; Ny = 60;
    dx = Lx / (Nx - 1); dy = Ly / (Ny - 1);
    t_final = 1800; dt = 1.0; N_steps = t_final / dt;

    Fx = (alpha * dt) / (2 * dx^2); Fy = (alpha * dt) / (2 * dy^2);
    Bix = (h * dx) / k; Biy = (h * dy) / k;

    % --- MODELAGEM DA GERAÇÃO DE CALOR NÃO-HOMOGÊNEA (GAUSSIANA) ---
    Q_gen = zeros(Nx, Ny);
    S = zeros(Nx, Ny);
    q_fundo = 30000;
    Q_tab = 150000;
    sigma = 0.015;

    % Coordenadas das abas (no topo da célula, y = Ly)
    x_tab1 = 0.02;
    x_tab2 = 0.08;
    y_tab = Ly;

    x_coord = linspace(0, Lx, Nx);
    y_coord = linspace(0, Ly, Ny);

    % Construção da Matriz de Geração de Calor
    for i = 1:Nx
        for j = 1:Ny
            dist_tab1 = (x_coord(i) - x_tab1)^2 + (y_coord(j) - y_tab)^2;
            dist_tab2 = (x_coord(i) - x_tab2)^2 + (y_coord(j) - y_tab)^2;

            gauss1 = exp(-dist_tab1 / (2 * sigma^2));
            gauss2 = exp(-dist_tab2 / (2 * sigma^2));

            Q_gen(i,j) = q_fundo + Q_tab * gauss1 + Q_tab * gauss2;
            S(i,j) = (Q_gen(i,j) * dt) / (2 * rho * cp);
        end
    end

    % --- Inicialização das Matrizes ---
    T = T_inf * ones(Nx, Ny);
    T_half = zeros(Nx, Ny);
    T_max_hist = zeros(N_steps, 1);
    time_vec = linspace(dt, t_final, N_steps);

    fprintf('Iniciando simulação ADI com Fonte Gaussiana (%d passos)...\n', N_steps);

    for n = 1:N_steps
        if mod(n, 300) == 0
            fprintf('Progresso: %d%% concluído...\n', round((n/N_steps)*100));
        end

        % ETAPA 1: Implícito em X
        for j = 1:Ny
            a = zeros(Nx,1); b = zeros(Nx,1); c = zeros(Nx,1); d = zeros(Nx,1);
            for i = 1:Nx
                if j == 1
                    RHS = Fy*T(i,j+1) + (1-2*Fy)*T(i,j) + Fy*(T(i,j+1) + 2*Biy*(T_inf - T(i,j))) + S(i,j);
                elseif j == Ny
                    RHS = Fy*(T(i,j-1) + 2*Biy*(T_inf - T(i,j))) + (1-2*Fy)*T(i,j) + Fy*T(i,j-1) + S(i,j);
                else
                    RHS = Fy*T(i,j-1) + (1-2*Fy)*T(i,j) + Fy*T(i,j+1) + S(i,j);
                end

                if i == 1
                    b(i) = 1 + 2*Fx + 2*Fx*Bix; c(i) = -2*Fx; d(i) = RHS + 2*Fx*Bix*T_inf;
                elseif i == Nx
                    a(i) = -2*Fx; b(i) = 1 + 2*Fx + 2*Fx*Bix; d(i) = RHS + 2*Fx*Bix*T_inf;
                else
                    a(i) = -Fx; b(i) = 1 + 2*Fx; c(i) = -Fx; d(i) = RHS;
                end
            end
            T_half(:, j) = thomas_algorithm(a, b, c, d, Nx);
        end

        % ETAPA 2: Implícito em Y
        for i = 1:Nx
            a = zeros(Ny,1); b = zeros(Ny,1); c = zeros(Ny,1); d = zeros(Ny,1);
            for j = 1:Ny
                if i == 1
                    RHS = Fx*T_half(i+1,j) + (1-2*Fx)*T_half(i,j) + Fx*(T_half(i+1,j) + 2*Bix*(T_inf - T_half(i,j))) + S(i,j);
                elseif i == Nx
                    RHS = Fx*(T_half(i-1,j) + 2*Bix*(T_inf - T_half(i,j))) + (1-2*Fx)*T_half(i,j) + Fx*T_half(i-1,j) + S(i,j);
                else
                    RHS = Fx*T_half(i-1,j) + (1-2*Fx)*T_half(i,j) + Fx*T_half(i+1,j) + S(i,j);
                end

                if j == 1
                    b(j) = 1 + 2*Fy + 2*Fy*Biy; c(j) = -2*Fy; d(j) = RHS + 2*Fy*Biy*T_inf;
                elseif j == Ny
                    a(j) = -2*Fy; b(j) = 1 + 2*Fy + 2*Fy*Biy; d(j) = RHS + 2*Fy*Biy*T_inf;
                else
                    a(j) = -Fy; b(j) = 1 + 2*Fy; c(j) = -Fy; d(j) = RHS;
                end
            end
            T(i, :) = thomas_algorithm(a, b, c, d, Ny)';
        end
        T_max_hist(n) = max(max(T));
    end

    fprintf('Simulação Finalizada!\n');

    % --- Gráficos ---
    [X, Y] = meshgrid(x_coord, y_coord);

    % Figura 1: Mapa de Calor da Temperatura
    figure(1);
    contourf(X, Y, T', 30, 'LineColor', 'none'); colormap('jet'); colorbar;
    title(['Distribuição de Temperatura (t = ', num2str(t_final), 's)']);
    xlabel('Lx (m)'); ylabel('Ly (m)');

    % Figura 2: Gráfico de Temperatura vs Tempo
    figure(2);
    plot(time_vec, T_max_hist, 'r', 'LineWidth', 2);
    title('Temperatura Máxima vs Tempo');
    xlabel('Tempo (s)'); ylabel('T_{max} (C)'); grid on;

    % Figura 3: Mapa da Geração de Calor (Para provar o modelo)
    figure(3);
    contourf(X, Y, Q_gen', 30, 'LineColor', 'none'); colormap('hot'); colorbar;
    title('Distribuição Espacial da Geração de Calor (W/m^3)');
    xlabel('Lx (m)'); ylabel('Ly (m)');
end

% --- Função Auxiliar: Algoritmo de Thomas (Aninhada Internamente) ---
function x = thomas_algorithm(a, b, c, d, n)
    c_star = zeros(n, 1); d_star = zeros(n, 1); x = zeros(n, 1);
    c_star(1) = c(1) / b(1);
    d_star(1) = d(1) / b(1);
    for i = 2:n
        m = 1.0 / (b(i) - a(i) * c_star(i-1));
        if i < n, c_star(i) = c(i) * m; end
        d_star(i) = (d(i) - a(i) * d_star(i-1)) * m;
    end
    x(n) = d_star(n);
    for i = n-1:-1:1, x(i) = d_star(i) - c_star(i) * x(i+1); end
end
