function ValidateMinimum
    % Blok kiri atas 2x2 dari Hessian terakhir
    H_reduced = [2.680386 -0.800000;
                -0.800000 0.670097];
    
    fprintf('Blok kiri atas 2x2 Hessian:\n');
    fprintf('[%.6f %.6f]\n', H_reduced(1,1), H_reduced(1,2));
    fprintf('[%.6f %.6f]\n\n', H_reduced(2,1), H_reduced(2,2));
    
    % Cek eigenvalue
    eig_vals = eig(H_reduced);
    fprintf('Eigenvalues:\n');
    fprintf('λ₁ = %.6f\n', eig_vals(1));
    fprintf('λ₂ = %.6f\n\n', eig_vals(2));
    
    if all(eig_vals > 0)
        fprintf('Semua eigenvalue positif → Hessian definit positif → Minimum lokal\n');
    else
        fprintf('Ada eigenvalue ≤ 0 → Hessian tidak definit positif → Bukan minimum lokal\n');
    end
end