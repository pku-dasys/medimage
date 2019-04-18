
function fupdate = minimize_in_IMAGE_variable_SD(g,f,v,A,AStar,alpha,epsilon,IM_TOL, IterImage)
% This function minimizes the Ambrosio-Tortorelli functional for a fixed
% v in the f variable.
%
% The functional is approximately minimized by the method of steepest descent
% with a fixed number of iteration steps.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input/Variables
% g                      l x k measurement
% f                      n x m reconstruction
% v                      n x m edge indicator
% alpha                  weight for smoothnes penalty
% A                      linear operator given as a function handle
% AStar                  adjoint operator of A given as a function handle
% IterImage              number of descend steps
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% fupdate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialize values:
global imageHight;
global imageWidth;

% Matrix for zero-padding of image and one-padding for edgeindicator
% at the boundaries
tempZeros = zeros(imageHight+2, imageWidth+2);
padded_v  = ones(imageHight+2, imageWidth+2);
padded_v(2:imageHight+1,2:imageWidth+1)  = v; %one padding
fupdate = f;


for k = 1:IterImage
    
    fold = fupdate;
    
    
    % d = - \nabla_f AT(f,v);
    % d = A*(g - Af) + 2alpha div(v^2 \nabla f) +0.01epsilon^2 LaPlace f:
    tempZeros(2:imageHight+1,2:imageWidth+1) = fold; %zero padding
    d = AStar(g - A(fold));
    for i=2:imageHight+1
        for j=2:imageWidth+1
            d(i-1,j-1) =d(i-1,j-1) + (alpha)*(...
                +padded_v(i,j)^2 *(tempZeros(i+1,j) -tempZeros(i,j))...
                +padded_v(i,j)^2 *(tempZeros(i,j+1) -tempZeros(i,j))...
                -padded_v(i-1,j)^2 *(tempZeros(i,j) -tempZeros(i-1,j))...
                -padded_v(i,j-1)^2 *(tempZeros(i,j) -tempZeros(i,j-1))...
                +0.01*epsilon^2 *...
                (tempZeros(i+1,j)+tempZeros(i,j+1) -4*tempZeros(i,j)+tempZeros(i-1,j)+tempZeros(i,j-1)));
        end
    end
    
    
    % Compute descent stepsize c_1
    % c_1 = <d,d>/(||Ad||^2 + alpha *||(v \nabla d)||^2 + alpha *||(sqrt(k_eps) \nabla d)||^2):
    
    tempZeros(2:imageHight+1,2:imageWidth+1) = d; %zero padding
    nV_KepsNablad       = 0;
    for i=2:imageHight+1
        for j=2:imageWidth+1
            nV_KepsNablad       = nV_KepsNablad +(padded_v(i,j)^2 + 0.01*epsilon^2)*(...
                ((tempZeros(i,j)- tempZeros(i-1,j)))^2 ...
                + ((tempZeros(i,j)- tempZeros(i,j-1)))^2 );
        end
    end
    
    dSkalard    = sum(sum(d.*d));
    normAd      = sum(sum(A(d).^2));
    c_1         = dSkalard/(normAd + alpha*nV_KepsNablad);
    

    
    
    % Descent step
    fupdate           = fold + c_1 * d; 
    
end

end
