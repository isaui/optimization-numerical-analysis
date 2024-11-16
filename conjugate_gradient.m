function conjugate_gradient()
    % Parameter
    x0 = 2;  % Tebakan awal
    max_iter = 2;  % Jumlah iterasi
    h = 1e-4;  % Step size untuk central difference
    
    fprintf('Tebakan awal:\n');
    fprintf('x_0 = %.10f\n', x0);
    fprintf('f(x_0) = %.10f\n', objective(x0));
    fprintf('f''(x_0) = %.10f\n', central_diff_first(x0, h));
    fprintf('|f''(x_0)| = %.10f\n\n', abs(central_diff_first(x0, h)));

    x = x0;
    g = central_diff_first(x, h);  % g_k
    p = -g;  % Arah pencarian awal sama dengan steepest descent
    
    for i = 1:max_iter
        % Line search untuk mencari α
        alpha = line_search_binary(x, p, h);
        
        % Update x
        x_new = x + alpha*p;
        
        % Hitung gradient baru (g_k+1)
        g_new = central_diff_first(x_new, h);
        
        % Hitung y_k = g_k+1 - g_k
        y_k = g_new - g;
        
        % Hitung beta (menggunakan formula dari slide: β_k = g_k+1^T * y_k / p_k^T * y_k)
        beta = (g_new' * y_k) / (p' * y_k);
        
        % Update arah pencarian untuk iterasi berikutnya
        p_new = -g_new + beta * p;
        
        % Hitung error optimal dari x yang baru
        error = abs(central_diff_first(x_new, h));
        
        fprintf('Iterasi CG %d:\n', i);
        fprintf('g_k = %.10f\n', g);
        fprintf('g_k+1 = %.10f\n', g_new);
        fprintf('y_k = %.10f\n', y_k);
        fprintf('β_k = %.10f\n', beta);
        fprintf('p_k = %.10f\n', p);
        fprintf('p_k+1 = %.10f\n', p_new);
        fprintf('α = %.10f\n', alpha);
        fprintf('x_%d = %.10f\n', i, x_new);
        fprintf('f(x_%d) = %.10f\n', i, objective(x_new));
        fprintf('|f''(x_%d)| = %.10f\n\n', i, error);
        
        % Update untuk iterasi berikutnya
        x = x_new;
        g = g_new;  % Update g_k menjadi g_k+1
        p = p_new;  % Update p_k menjadi p_k+1
    end
end

function alpha = line_search_binary(x0, p, h)
    % Binary search untuk line search
    a1 = 0;
    a2 = 1;
    tol = 1e-3;  % Tolerance untuk binary search
    
    fprintf('\nLine Search Iterasi Manual:\n');
    % Tampilkan 2 iterasi pertama manual
    for i = 1:2
        m = (a1 + a2)/2;
        
        % Hitung x di berbagai titik untuk central difference
        x1p = x0 + (a1 + h)*p;
        x1m = x0 + (a1 - h)*p;
        xmp = x0 + (m + h)*p;
        xmm = x0 + (m - h)*p;
        
        % Nilai x di titik a1, m, dan a2
        x_a1 = x0 + a1*p;
        x_m = x0 + m*p;
        x_a2 = x0 + a2*p;
        
        df1 = (objective(x1p) - objective(x1m))/(2*h);
        dfm = (objective(xmp) - objective(xmm))/(2*h);
        
        fprintf('Iterasi Line Search %d:\n', i);
        fprintf('a1 = %.10f, x(a1) = %.10f\n', a1, x_a1);
        fprintf('f(x(a1)) = %.10f\n', objective(x_a1));
        fprintf('f''(x(a1)) = %.10f\n', central_diff_first(x_a1, h));
        
        fprintf('x1p = %.10f\n', x1p);
        fprintf('x1m = %.10f\n', x1m);
        fprintf('f(x1p) = %.10f\n', objective(x1p));
        fprintf('f(x1m) = %.10f\n', objective(x1m));
        
        fprintf('m = %.10f, x(m) = %.10f\n', m, x_m);
        fprintf('f(x(m)) = %.10f\n', objective(x_m));
        fprintf('f''(x(m)) = %.10f\n', central_diff_first(x_m, h));
        
        fprintf('xmp = %.10f\n', xmp);
        fprintf('xmm = %.10f\n', xmm);
        fprintf('f(xmp) = %.10f\n', objective(xmp));
        fprintf('f(xmm) = %.10f\n', objective(xmm));
        
        fprintf('a2 = %.10f, x(a2) = %.10f\n', a2, x_a2);
        fprintf('f(x(a2)) = %.10f\n', objective(x_a2));
        fprintf('f''(x(a2)) = %.10f\n', central_diff_first(x_a2, h));
        
        fprintf('df1 = %.10f\n', df1);
        fprintf('dfm = %.10f\n', dfm);
        fprintf('df1*dfm = %.10f\n\n', df1*dfm);
        
        if df1*dfm < 0
            a2 = m;
        else
            a1 = m;
        end
    end
    
    % Lanjutkan binary search sampai mencapai toleransi
    while (a2 - a1) > tol
        m = (a1 + a2)/2;
        
        x1p = x0 + (a1 + h)*p;
        x1m = x0 + (a1 - h)*p;
        df1 = (objective(x1p) - objective(x1m))/(2*h);
        
        xmp = x0 + (m + h)*p;
        xmm = x0 + (m - h)*p;
        dfm = (objective(xmp) - objective(xmm))/(2*h);
        
        if df1*dfm < 0
            a2 = m;
        else
            a1 = m;
        end
    end
    
    alpha = (a1 + a2)/2;
    x_final = x0 + alpha*p;
    fprintf('Line search selesai dengan:\n');
    fprintf('α = %.10f\n', alpha);
    fprintf('x = %.10f\n', x_final);
    fprintf('f(x) = %.10f\n', objective(x_final));
    fprintf('f''(x) = %.10f\n\n', central_diff_first(x_final, h));
end

function obj = objective(x)
    % A(1,3) dan B(5,-1)
    a = 1; b = 3;
    c = 5; d = -1;
    % v_udara : v_kaca = 3 : 2
    obj = (1/3)*sqrt((x-a)^2 + b^2) + (1/2)*sqrt((x-c)^2 + d^2);
end

function df = central_diff_first(x, h)
    df = (objective(x + h) - objective(x - h))/(2*h);
end