function second_derivative_analysis()
    % Titik-titik yang diprediksi oleh berbagai metode
    x_newton = 4.367354;
    x_steepest = 2.697613;
    x_conjugate = 2.368752;
    x_quasi = 4.421619;
    
    % Step size untuk aproksimasi turunan kedua
    h = 1e-6;
    
    % Hitung turunan kedua untuk setiap titik
    f2_newton = calculate_second_derivative(x_newton, h);
    f2_steepest = calculate_second_derivative(x_steepest, h);
    f2_conjugate = calculate_second_derivative(x_conjugate, h);
    f2_quasi = calculate_second_derivative(x_quasi, h);
    
    % Tampilkan hasil
    fprintf('Analisis Turunan Kedua:\n');
    fprintf('1. Metode Newton (x = %.6f):\n', x_newton);
    fprintf('   f\"(x) = %.6f\n\n', f2_newton);
    
    fprintf('2. Steepest Descent (x = %.6f):\n', x_steepest);
    fprintf('   f\"(x) = %.6f\n\n', f2_steepest);
    
    fprintf('3. Conjugate Gradient (x = %.6f):\n', x_conjugate);
    fprintf('   f\"(x) = %.6f\n\n', f2_conjugate);
    
    fprintf('4. Quasi Newton (x = %.6f):\n', x_quasi);
    fprintf('   f\"(x) = %.6f\n\n', f2_quasi);
end

function f2 = calculate_second_derivative(x, h)
    % Menghitung turunan kedua menggunakan aproksimasi
    % f''(x) ≈ (f'(x+h) - f'(x-h))/(2h)
    f2 = (derivative(x + h) - derivative(x - h))/(2*h);
end

function df = derivative(x)
    % Turunan pertama yang sudah dihitung analitik
    a = 1; b = 3; % titik A(1,3)
    c = 5; d = -1; % titik B(5,-1)
    df = (x-a)/(3*sqrt((x-a)^2 + b^2)) + (x-c)/(2*sqrt((x-c)^2 + d^2));
end