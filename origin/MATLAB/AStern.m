function ASternf = AStern(f)
% Operator A

ASternf =   radon_adjoint(f, 0:1:179);

end