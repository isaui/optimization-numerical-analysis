function LagrangeNewton
    % Inisialisasi
    clear; clc;
    
    x0 = 2*sqrt(2);
    y0 = 0;
    lambda0 = 1;
    z0 = [x0; y0; lambda0];
    
    max_iter = 12;
    errors = zeros(max_iter, 1);
    
    z = z0;
    fprintf('Iterasi 0: x = %.4f, y = %.4f, lambda = %.4f\n', z(1), z(2), z(3));
    
    grad_init = gradientL(z(1), z(2), z(3));
    error_init = norm(grad_init);
    
    % Debug error awal per komponen
    fprintf('Error awal total = %.6f\n', error_init);
    fprintf('Error awal per arah:\n');
    fprintf('  dL/dx = %.6f\n', grad_init(1));
    fprintf('  dL/dy = %.6f\n', grad_init(2));
    fprintf('  dL/dλ = %.6f\n', grad_init(3));
    fprintf('Kendala = %.6f\n\n', constraintValue(z(1), z(2)));
    
    % Debug Hessian awal
    H_init = hessianL(z(1), z(2), z(3));
    fprintf('Hessian awal:\n');
    fprintf('[%.6f %.6f %.6f]\n', H_init(1,1), H_init(1,2), H_init(1,3));
    fprintf('[%.6f %.6f %.6f]\n', H_init(2,1), H_init(2,2), H_init(2,3));
    fprintf('[%.6f %.6f %.6f]\n\n', H_init(3,1), H_init(3,2), H_init(3,3));
    
    for k = 1:max_iter
        grad = gradientL(z(1), z(2), z(3));
        Hk = hessianL(z(1), z(2), z(3));
        
        fprintf('Iterasi %d:\n', k);
        
        % Debug Hessian SEBELUM update z
        fprintf('Hessian:\n');
        fprintf('[%.6f %.6f %.6f]\n', Hk(1,1), Hk(1,2), Hk(1,3));
        fprintf('[%.6f %.6f %.6f]\n', Hk(2,1), Hk(2,2), Hk(2,3));
        fprintf('[%.6f %.6f %.6f]\n\n', Hk(3,1), Hk(3,2), Hk(3,3));
        
        z_new = z - Hk\grad;
        
        grad_new = gradientL(z_new(1), z_new(2), z_new(3));
        error = norm(grad_new);
        errors(k) = error;
        
        fprintf('x = %.6f, y = %.6f, lambda = %.6f\n', z_new(1), z_new(2), z_new(3));
        fprintf('Error total = %.6f\n', error);
        fprintf('Error per arah:\n');
        fprintf('  dL/dx = %.6f\n', grad_new(1));
        fprintf('  dL/dy = %.6f\n', grad_new(2));
        fprintf('  dL/dλ = %.6f\n', grad_new(3));
        fprintf('Kendala = %.6f\n\n', constraintValue(z_new(1), z_new(2)));
        
        z = z_new;
    end
    
    % Plot konvergensi
    figure;
    semilogy(0:max_iter, [error_init; errors], 'b-o');
    grid on;
    xlabel('Iterasi');
    ylabel('||∇L|| (log scale)');
    title('Konvergensi Metode Newton - Norm Gradien');
end

% Fungsi-fungsi pendukung tetap sama
function grad = gradientL(x, y, lambda)
    dL_dx = -4*(-2*x + y + 20)/5 + lambda*x/4;
    dL_dy = 2*(-2*x + y + 20)/5 + lambda*y/16;
    dL_dlambda = x^2/8 + y^2/32 - 1;
    
    grad = [dL_dx; dL_dy; dL_dlambda];
end

function H = hessianL(x, y, lambda)
    % Hitung turunan kedua untuk blok kiri atas 2x2
    d2f_dx2 = 8/5;
    d2f_dxdy = -4/5;
    d2f_dy2 = 2/5;
    
    d2g_dx2 = 1/4;
    d2g_dy2 = 1/16;
    d2g_dxdy = 0;
    
    % Turunan kedua terhadap lambda dan x,y
    d2L_dldx = x/4;    % ∂/∂λ(∂L/∂x) = ∂/∂x(∂L/∂λ) = x/4
    d2L_dldy = y/16;   % ∂/∂λ(∂L/∂y) = ∂/∂y(∂L/∂λ) = y/16
    
    H = [d2f_dx2 + lambda*d2g_dx2,  d2f_dxdy + lambda*d2g_dxdy,  d2L_dldx;
         d2f_dxdy + lambda*d2g_dxdy, d2f_dy2 + lambda*d2g_dy2,   d2L_dldy;
         d2L_dldx,                   d2L_dldy,                    0      ];
end
function c = constraintValue(x, y)
    c = x^2/8 + y^2/32 - 1;
end