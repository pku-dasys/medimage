function [f,v]=MumfordShah_reg(g,alpha,beta,epsilon,initialImage,initialEdge,IterAlt,IterImage,IterEdge,IM_TOL,ED_TOL)

% [f,v]= mumfordshah_reg(g,initialImage,initialEdge,Tag,A,AStar,alpha,beta,epsilon,AlternateIt,IM_TOL,ED_TOL,ItImage,ItEdge)
% computes a local minimizer of the Ambrosio-Tortorelli functional
% AT_epsilon(f,v) = ||Af - g||^2
%      + alpha*int_(\Omega) (v^2+0.01epsilon^2)|\nabla f|^2 dx
%      + beta *int_(\Omega) epsilon|\nabla v|^2 + 1/(4*epsilon)*(1-v)^2 dx
%
% by an alternating minimization in the image variable f and the edge
% indicator v. If epsilon -> 0, then global minimizers of AT_epsilon
% approximate minimizers of the Mumford-Shah functional
% MS(f,K) = ||Af - g||^2
%           + alpha* int_(\Omega/K)|\nabla f|^2 dx
%           + beta*  Length(K)
%
% where f is an image and K the edge set.
%
% The algorithm alternatingly keeps one variable of (f,v) fixed and
% minimizes the functional in the other one. This is done for a fixed number of
% alternating steps.
% For the minimization in the image variable the method of conjugate
% gradients is used with a fixed number of descent steps.
% For the minimization in the edge variable the method of steepest
% descent with a fixed number of descent steps is used. After each descent
% step the current edge is projected to the admissible set of functions
% bound by 0 and 1 from below and above respectively.
% For both methods the stepsize is chosen by an exact line search.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Input
% g                      l x k matrix measurement
% A                      linear operator given as a function handle
% AStar                  adjoint operator of A given as a function handle
% alpha                  weight for smoothnes penalty
% beta                   weight for edge penalty
%
% Optional Input
% initialImage           n x m matrix
% initialEdge            n x m matrix 
%                        else initial image = zeros; initial edge = ones.
% epsilon                gamma convergence parameter:
%                        AT_epsilon Gamma -> MS as espilon ->0
% IterAlt                number of alternating steps
% IM_TOL                 tolerance for each minimization in image variable
% ED_TOL                 tolerance for each minimization in edge variable
% ItImage                max number of cg-descend steps in image variable
% ItEdge                 max number of descend steps in edge variablen
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
% f                     n x m matrix reconstruction
% v                     n x m matrix edge indicator
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example for segmentation (A = I):
%  g =  phantom(256) + 0.05*randn(256);
%  I = @ (x) 1*x;
%  [f,v] = MumfordShah_reg(g,I,I,15,0.000005,0.001,g,[],5,15,15,0.001,0.001);
%  subplot(1,3,1)
%  imagesc(g)
%  subplot(1,3,2)
%  imagesc(f)
%  subplot(1,3,3)
%  imagesc(v)
%  colormap gray
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Initialization:
Temp       = AStern(g);
imageHight = size(Temp,1);
imageWidth = size(Temp,2);
dataHight  = size(g,1);
dataWidth  = size(g,2); 

%Check optional inputs
if ~exist('epsilon','var') || isempty(epsilon)
    epsilon = 0.0001; 
end

if ~exist('initialImage','var')  || isempty(initialImage)
     f   = zeros(imageHight,imageWidth);
else
     f  = initialImage;
end

if ~exist('initialEdge','var')  || isempty(initialEdge)  
    v = ones(imageHight,imageWidth);
else
    v = initialEdge;
end 

if ~exist('IterAlt','var') || isempty(IterAlt)  
    IterAlt = 10; 
end

if ~exist('IterImage','var') || isempty(IterImage)  
    IterImage = 10;  
end

if ~exist('IterEdge','var') || isempty(IterEdge)  
    IterEdge = 10;  
end

if ~exist('IM_TOL','var')   || isempty(IM_TOL)
    IM_TOL = 0.0001;  
end

if ~exist('ED_TOL','var')   || isempty(ED_TOL)
    ED_TOL = 0.0001;  
end


% The Alternat minimization
fprintf('\n Begin alternative minimization...\n')
for i=1:1%IterAlt
    fprintf('Iteration %d...\n',i);
    
fid=fopen('my_g.txt','wt');
for i=1:128
    for j=1:512
        fprintf(fid,'%f ',g(i,j));
    end
    fprintf(fid,'\n');
end
    f   = minimize_in_IMAGE_variable(g,f,v,imageHight,imageWidth,alpha,epsilon,IM_TOL,IterImage);
    
    v   = minimize_in_EDGE_variable(f,v,imageHight,imageWidth,alpha,beta,epsilon,ED_TOL,IterEdge);
    
end
fprintf('\n Leaving MumfordShah_reg...\n')
return;
end

