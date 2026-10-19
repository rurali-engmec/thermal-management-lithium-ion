function projeto_ftc()
    % =========================================================================
    % SIMULAÇÃO TÉRMICA 2D - CÉLULA LFP (ADI) COM INDICADOR DE PROGRESSO
    % =========================================================================
    clc; close all;

    % --- Parâmetros ---
    Lx = 0.10; Ly = 0.15; rho = 2240; cp = 800; k = 3.0;
    alpha = k / (rho * cp); q_gen = 65000; T_inf = 25.0; h = 166.0;

    Nx = 10; Ny = 20;
    dx = Lx / (Nx - 1); dy = Ly / (Ny - 1);
    t_final = 1800; dt = 1.0; N_steps = t_final / dt;

    Fx = (alpha * dt) / (2 * dx^2); Fy = (alpha * dt) / (2 * dy^2);
    S  = (q_gen * dt) / (2 * rho * cp);
    Bix = (h * dx) / k; Biy = (h * dy) / k;

    T = T_inf * ones(Nx, Ny);
    T_half = zeros(Nx, Ny);
    T_max_hist = zeros(N_steps, 1);
    time_vec = linspace(dt, t_final, N_steps);

    fprintf('Iniciando simulação ADI (%d passos)...\n', N_steps);

    for n = 1:N_steps
        % --- INDICADOR DE PROGRESSO NO TERMINAL ---
        if mod(n, 300) == 0
            fprintf('Progresso: %d%% concluído...\n', round((n/N_steps)*100));
        end

        % ETAPA 1: Implícito em X
        for j = 1:Ny
            a = zeros(Nx,1); b = zeros(Nx,1); c = zeros(Nx,1); d = zeros(Nx,1);
            for i = 1:Nx
                if j == 1
                    RHS = Fy*T(i,j+1) + (1-2*Fy)*T(i,j) + Fy*(T(i,j+1) + 2*Biy*(T_inf - T(i,j))) + S;
                elseif j == Ny
                    RHS = Fy*(T(i,j-1) + 2*Biy*(T_inf - T(i,j))) + (1-2*Fy)*T(i,j) + Fy*T(i,j-1) + S;
                else
                    RHS = Fy*T(i,j-1) + (1-2*Fy)*T(i,j) + Fy*T(i,j+1) + S;
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
                    RHS = Fx*T_half(i+1,j) + (1-2*Fx)*T_half(i,j) + Fx*(T_half(i+1,j) + 2*Bix*(T_inf - T_half(i,j))) + S;
                elseif i == Nx
                    RHS = Fx*(T_half(i-1,j) + 2*Bix*(T_inf - T_half(i,j))) + (1-2*Fx)*T_half(i,j) + Fx*T_half(i-1,j) + S;
                else
                    RHS = Fx*T_half(i-1,j) + (1-2*Fx)*T_half(i,j) + Fx*T_half(i+1,j) + S;
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

    % Gráficos
    [X, Y] = meshgrid(linspace(0, Lx, Nx), linspace(0, Ly, Ny));
    figure(1);
    contourf(X, Y, T', 20, 'LineColor', 'none'); colormap('jet'); colorbar;
    title(['Mapa de Calor (t = ', num2str(t_final), 's)']);
    xlabel('Lx (m)'); ylabel('Ly (m)');

    figure(2);
    plot(time_vec, T_max_hist, 'r', 'LineWidth', 2);
    title('Temperatura Máxima vs Tempo');
    xlabel('Tempo (s)'); ylabel('T_m_a_x (C)'); grid on;
end

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
