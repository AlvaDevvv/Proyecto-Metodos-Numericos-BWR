% =========================================================================
% Archivo: Escenario3_Falla.m
% Descripción: Falla de refrigeración (LOCA). Cae drásticamente el 
% coeficiente convectivo (h_inf) provocando un sobrecalentamiento rápido.
% =========================================================================
clear; clc; close all;

% 1. Cargar la configuración base (Plantilla)
Parametros_y_Propiedades;

% 2. Sobrescritura Local de Parámetros (Escenario 3)
q_gen = 3e8;    % Sigue generando calor a tope
h_inf = 1000;   % Cae drásticamente, el agua ya no roba calor eficientemente
T_inf = 560;    
T_ini = 600;    

% 3. Parámetros de Discretización
Nr = 50; 
r = linspace(0, R3, Nr)';
dr = r(2) - r(1);

dt = 0.5;        
t_final = 200;   % Se simula más tiempo (200s) para ver el desastre térmico 
t_vec = 0:dt:t_final;
Nt = length(t_vec);

idx_c = r <= R1;                       
idx_h = (r > R1) & (r <= R2);          
idx_v = r > R2;                        

% 4. Inicialización
T_ant = T_ini * ones(Nr, 1); 
T_hist = zeros(Nr, Nt);      
T_hist(:,1) = T_ant;
tol_picard = 1e-4; 

% 5. Bucle Temporal
for n = 2:Nt
    T_iter = T_ant; 
    error_picard = 1;
    iteraciones = 0;
    
    while error_picard > tol_picard && iteraciones < 50
        iteraciones = iteraciones + 1;
        
        a_diag = zeros(Nr, 1); b_diag = zeros(Nr, 1);
        c_diag = zeros(Nr, 1); d_vec  = zeros(Nr, 1);
        k_eval = zeros(Nr, 1); rhoC_eval = zeros(Nr, 1); q_eval = zeros(Nr, 1);
        
        for i = 1:Nr
            if idx_c(i)
                k_eval(i) = k_c(T_iter(i)); rhoC_eval(i) = rhoC_c(T_iter(i)); q_eval(i) = q_gen; 
            elseif idx_h(i)
                k_eval(i) = k_h(T_iter(i)); rhoC_eval(i) = rhoC_h(T_iter(i)); q_eval(i) = 0;
            else
                k_eval(i) = k_v(T_iter(i)); rhoC_eval(i) = rhoC_v(T_iter(i)); q_eval(i) = 0;
            end
        end
        
        for i = 1:Nr
            if i == 1 
                k_int = (2 * k_eval(1) * k_eval(2)) / (k_eval(1) + k_eval(2));
                b_diag(1) = rhoC_eval(1)/dt + 4*k_int/(dr^2);
                c_diag(1) = 4*k_int/(dr^2);
                d_vec(1)  = rhoC_eval(1)/dt * T_ant(1) + q_eval(1);
            elseif i == Nr 
                k_int = (2 * k_eval(Nr-1) * k_eval(Nr)) / (k_eval(Nr-1) + k_eval(Nr));
                A_cond = k_int * (r(Nr) - dr/2) / (r(Nr) * dr);
                A_conv = h_inf / dr;
                a_diag(Nr) = A_cond;
                b_diag(Nr) = rhoC_eval(Nr)/(2*dt) + A_cond + A_conv;
                d_vec(Nr)  = rhoC_eval(Nr)/(2*dt) * T_ant(Nr) + A_conv * T_inf;
            else 
                r_mas  = r(i) + dr/2; r_menos = r(i) - dr/2;
                k_mas   = (2 * k_eval(i) * k_eval(i+1)) / (k_eval(i) + k_eval(i+1));
                k_menos = (2 * k_eval(i-1) * k_eval(i)) / (k_eval(i-1) + k_eval(i));
                coef_E = (k_mas * r_mas) / (r(i) * dr^2); coef_W = (k_menos * r_menos) / (r(i) * dr^2);
                a_diag(i) = coef_W; c_diag(i) = coef_E;
                b_diag(i) = rhoC_eval(i)/dt + coef_W + coef_E;
                d_vec(i)  = rhoC_eval(i)/dt * T_ant(i) + q_eval(i);
            end
        end
        T_new = Solver_TDMA(-a_diag, b_diag, -c_diag, d_vec);
        error_picard = norm(T_new - T_iter, inf);
        T_iter = T_new; 
    end
    T_ant = T_new; T_hist(:,n) = T_new;
end

% 6. Gráficos
figure;
plot(r*1000, T_hist(:,1), '--k', 'DisplayName', 't = 0 s'); hold on;
plot(r*1000, T_hist(:, floor(Nt/4)), 'b', 'DisplayName', ['t = ', num2str(t_vec(floor(Nt/4))), ' s']);
plot(r*1000, T_hist(:, floor(Nt/2)), 'g', 'DisplayName', ['t = ', num2str(t_vec(floor(Nt/2))), ' s']);
plot(r*1000, T_hist(:, end), 'r', 'LineWidth', 1.5, 'DisplayName', ['t = ', num2str(t_final), ' s']);

xline(R1*1000, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');
xline(R2*1000, '--', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');
max_T = max(T_hist(:));
text(R1*1000/2, max_T*0.75, 'Pastilla UO_2', 'HorizontalAlignment','center');
text((R1+R2)*1000/2, max_T*0.85, 'Huelgo', 'HorizontalAlignment','center');
text((R2+R3)*1000/2, max_T*0.75, 'Vaina', 'HorizontalAlignment','center');

xlabel('Radio (mm)'); ylabel('Temperatura (K)');
title('Escenario 3: Falla de Refrigerante (Sobrecalentamiento)');
legend('Location', 'best'); grid on;