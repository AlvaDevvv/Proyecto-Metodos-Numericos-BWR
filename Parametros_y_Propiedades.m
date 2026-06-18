% =========================================================================
% Archivo: Parametros_y_Propiedades.m
% Descripción: Definición de la geometría, condiciones iniciales y 
% funciones de propiedades termofísicas (Tabla 2).
% =========================================================================

% 1. Parámetros Geométricos (en metros)
R1 = 4.3815e-3; % Radio exterior de la pastilla
R2 = 4.4705e-3; % Radio interior de la vaina
R3 = 5.131e-3;  % Radio exterior de la vaina

% 2. Condiciones de Operación Nominales
q_gen = 3e8;    % Generación de calor en la pastilla [W/m^3] 
h_inf = 50000;  % Coeficiente de convección [W/m^2K] 
T_inf = 560;    % Temperatura del refrigerante [K] 
T_ini = 600;    % Temperatura inicial del sistema [K]

% 3. Funciones de Propiedades (Dependientes de la Temperatura T en Kelvin)

% --- PASTILLA (UO2) ---
% Usamos real() para ignorar cualquier error de redondeo complejo
k_c = @(T) real(3.82502e3 ./ (T + 129.411) + 6.08011e-11 .* T.^3);

rhoC_c = @(T) real(8.510322e11 .* exp(535.285./T) ./ (T .* (exp(535.285./T) - 1).^2) + ...
               2.434842e2 .* T + 1.660985e16 .* exp(-18970.61./T) ./ T.^2);

% --- HUELGO (Helio) ---
k_h = @(T) real(2.517e-3 .* T.^0.72);
rhoC_h = @(T) 1e-6; 

% --- VAINA (Zircaloy-2) ---
k_v = @(T) real(7.51 + 2.09e-2.*T - 1.45e-5.*T.^2 + 7.67e-9.*T.^3);

rhoC_v = @(T) real(1.820453e6 + 3.038627e5.*((T - 300)/200) - ...
              1.063741e5.*((T - 300)/200).^2 + ...
              2.810287e4.*((T - 300)/200).^3 - ...
              2.723618e3.*((T - 300)/200).^4);