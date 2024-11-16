function compare_snell_law()
   % Titik-titik x dari berbagai metode
   x_newton = 4.367354;
   x_steepest = 2.697613;
   x_conjugate = 2.368752;
   x_quasi = 4.421619;
   
   % Koordinat titik-titik tetap
   a = 1; b = 3;    % A(a,b)
   c = 5; d = -1;   % B(c,d)
   
   % Hitung dan tampilkan untuk setiap metode
   fprintf('Analisis Hukum Snell:\n');
   analyze_point('Newton', x_newton, a, b, c, d);
   analyze_point('Steepest Descent', x_steepest, a, b, c, d);
   analyze_point('Conjugate Gradient', x_conjugate, a, b, c, d);
   analyze_point('Quasi-Newton', x_quasi, a, b, c, d);
end

function analyze_point(method_name, x, a, b, c, d)
   % Untuk θ₁:
   % sin(θ₁) = opposite/hypotenuse = (x-a)/√((x-a)² + b²)
   sin_theta1 = abs(x-a)/sqrt((x-a)^2 + b^2);
   
   % Untuk θ₂:
   % sin(θ₂) = opposite/hypotenuse = (c-x)/√((c-x)² + d²)
   sin_theta2 = abs(c-x)/sqrt((c-x)^2 + d^2);
   
   % Hitung rasio
   ratio = sin_theta1/sin_theta2;
   
   % Tampilkan hasil
   fprintf('\n%s Method:\n', method_name);
   fprintf('sin(θ₁) = %.6f\n', sin_theta1);
   fprintf('sin(θ₂) = %.6f\n', sin_theta2);
   fprintf('sin(θ₁)/sin(θ₂) = %.6f\n', ratio);
   fprintf('Error dari 3/2 = %.6f%%\n', abs(ratio - 1.5)/1.5 * 100);
end