function newton()
    % Main script untuk metode Newton
    x0 = 2;  % Tebakan awal
    max_iter = 2;  % Jumlah iterasi yang diminta
    h = 1e-6; % Step size untuk central difference

    % Debug print untuk tebakan awal
    fprintf('Tebakan awal:\n');
    fprintf('x_0 = %.10f\n', x0);
    fprintf('f(x_0) = %.10f\n', objective(x0));
    fprintf('f''(x_0) = %.10f\n', derivative(x0));
    fprintf('|f''(x_0)| = %.10f\n\n', abs(derivative(x0)));

    x = x0;  % x sekarang
    for i = 1:max_iter
        % Hitung turunan pertama (eksak) dan kedua (aproksimasi dari turunan pertama)
        df = derivative(x);
        d2f = (derivative(x + h) - derivative(x - h))/(2*h);
        
        % Update x menggunakan metode Newton
        x_new = x - df/d2f;
        
        % Hitung error optimal dari x yang baru menggunakan turunan eksak
        error = abs(derivative(x_new));
        
        % Debug print
        fprintf('Iterasi %d:\n', i);
        fprintf('f''(x) = %.10f\n', df);
        fprintf('f''''(x) = %.10f\n', d2f);
        fprintf('x_%d = %.10f\n', i, x_new);
        fprintf('f(x_%d) = %.10f\n', i, objective(x_new));
        fprintf('f''(x_%d) = %.10f\n', i, derivative(x_new));
        fprintf('|f''(x_%d)| = %.10f\n\n', i, error);
        
        x = x_new;
    end
end

function obj = objective(x)
    % A(1,3) dan B(5,-1)
    a = 1; b = 3;
    c = 5; d = -1;
    % v_udara : v_kaca = 3 : 2
    obj = (1/3)*sqrt((x-a)^2 + b^2) + (1/2)*sqrt((x-c)^2 + d^2);
end

function df = derivative(x)
    % Turunan pertama yang sudah dihitung secara analitik
    a = 1; b = 3;
    c = 5; d = -1;
    df = (x-a)/(3*sqrt((x-a)^2 + b^2)) + (x-c)/(2*sqrt((x-c)^2 + d^2));
end