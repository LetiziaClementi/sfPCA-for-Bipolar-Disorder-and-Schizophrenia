function [norm] = getL2Norm(f, R0)
       % f vector of ceofficients
       % R0, R0 matrix of where the solution has been computed
    f = reshape(f,[],1);
	norm = sqrt(f'*R0*f);
end