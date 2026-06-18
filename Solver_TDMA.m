function T_new = Solver_TDMA(a, b, c, d)
% =========================================================================
% Archivo: Solver_TDMA.m
% Descripción: Algoritmo de Thomas para resolver sistemas tridiagonales.
% a: diagonal inferior (subdiagonal)
% b: diagonal principal
% c: diagonal superior (superdiagonal)
% d: vector de términos independientes
% =========================================================================
    n = length(d);
    c_star = zeros(n, 1);
    d_star = zeros(n, 1);
    T_new = zeros(n, 1);

    % Forward sweep
    c_star(1) = c(1) / b(1);
    d_star(1) = d(1) / b(1);
    
    for i = 2:n
        denom = b(i) - a(i) * c_star(i-1);
        if i < n
            c_star(i) = c(i) / denom;
        end
        d_star(i) = (d(i) - a(i) * d_star(i-1)) / denom;
    end

    % Back substitution
    T_new(n) = d_star(n);
    for i = n-1:-1:1
        T_new(i) = d_star(i) - c_star(i) * T_new(i+1);
    end
end