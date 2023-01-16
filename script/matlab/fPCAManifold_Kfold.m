50function [CVseq] = fPCAManifold_Kfold(X, mass, stiff, loglambdaseq, K, n_iter, indices_NAN)
    
  X_clean = X;
  X_clean(:,indices_NAN) = [];
  
  %np = size(X_clean,2);
  %X X_clean
  nbasis = size(X,2);
  
  nlambda = length(loglambdaseq);
  %f_hat_mat = zeros(nbasis, nlambda);
  CVseq = zeros(1,nlambda);
  
  % Compute solution for each lambda
  parfor ilambda = 1:nlambda
      NW = spdiags(ones(size(mass,2),1), 0, size(mass,2), size(mass,2));
      NW(indices_NAN,indices_NAN) = 0;
      SW = (10^(loglambdaseq(ilambda)))*stiff;
      NE = (10^(loglambdaseq(ilambda)))*stiff;
      SE = -(10^(loglambdaseq(ilambda)))*mass;
      A = [NW NE ; SW SE]; %Compute left hand side linear system paper
      
    for ifold = 1:K
        %Split dataset into training set and test set
        length_chunk = floor(size(X,1)/K);
        indices_valid = (ifold-1)*length_chunk+(1:length_chunk);
        indices_train = setdiff(1:size(X,1), indices_valid);
        X_train = X(indices_train,:); 
        X_train = X_train - ones(size(X_train,1),1)*mean(X_train);
        %X_valid = X(indices_valid,:);
        %X_valid = X_valid - ones(size(X_valid,1),1)*mean(X_valid);
        X_clean_train = X_clean(indices_train,:);
        X_clean_train = X_clean_train - ones(size(X_clean_train,1),1)*mean(X_clean_train);
        X_clean_valid = X_clean(indices_valid,:);
        X_clean_valid = X_clean_valid - ones(size(X_clean_valid,1),1)*mean(X_clean_valid);
        [U,~,~] = svd(X_clean_train, 'econ');
        u_hat = U(:,1);
        fp = 0;
        g = 0;
          %Compute solution fo each training set
          for iter = 1:n_iter  
            b = [X_train' * u_hat; zeros(nbasis,1)];
            b(indices_NAN) = 0;
            %sol = U\(L\b);
            sol = (A \ b);
            f = sol(1:nbasis);
            g = sol((1+nbasis):2*nbasis);
            %disp(b(1:3)-f(1:3))
            fp = f;
            fp(indices_NAN) = []; % questa i loadings output della fPCA
            u_hat = X_clean_train*fp;
            u_hat = u_hat/sqrt(sum(u_hat.^2)); %questa è il risultato della fPCA normale (gli score output)
            %disp(f(1:3))
          end
      
      u_hat_const = sum(fp.^2) + 10^(loglambdaseq(ilambda))* (g' *mass*g);
      u_hat_valid = (X_clean_valid*fp)./u_hat_const;
      CV_local = min([sum(sum((X_clean_valid - u_hat_valid*(fp')).^2))/numel(X_clean_valid) sum(sum((X_clean_valid + u_hat_valid*(fp')).^2))/numel(X_clean_valid)]);
      CVseq(ilambda) = CVseq(ilambda) + CV_local;
    end
    disp(['     CV index computed for lambda st log(lambda)=' num2str(loglambdaseq(ilambda))])
  end
end
