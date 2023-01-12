function [f_hat_mat, u_hat_mat] = fPCAManifold(X, mass, stiff, loglambdaseq, n_iter, indices_NAN)
    
  X_clean = X;
  X_clean(:,indices_NAN) = [];
  
  np = size(X_clean,2);
  %X X_clean
  nbasis = size(X,2);
  
  nlambda = length(loglambdaseq);
  f_hat_mat = zeros(nbasis, nlambda);
  u_hat_mat = zeros(size(X,1), nlambda);
  
  for ilambda = 1:nlambda
      NW = spdiags(ones(size(mass,2),1), 0, size(mass,2), size(mass,2));
      NW(indices_NAN,indices_NAN) = 0;
      SW = (10^(loglambdaseq(ilambda)))*stiff;
      NE = (10^(loglambdaseq(ilambda)))*stiff;
      SE = -(10^(loglambdaseq(ilambda)))*mass;
      A = [NW NE ; SW SE];
      [U,S,V] = svd(X_clean, 'econ');
      u_hat = U(:,1);
      for iter = 1:n_iter  
        b = [X' * u_hat; zeros(nbasis,1)];
        b(indices_NAN) = 0;
        sol = (A \ b);
        f = sol(1:nbasis);
        g = sol((1+nbasis):2*nbasis);
        %disp(b(1:3)-f(1:3))
        fp = f;
        fp(indices_NAN) = [];
        u_hat = X_clean*fp;
        u_hat = u_hat/sqrt(sum(u_hat.^2));
        %disp(f(1:3))
      end
      u_hat_mat(:,ilambda) = u_hat;
      f_hat_mat(:,ilambda) = f;
  end
end
