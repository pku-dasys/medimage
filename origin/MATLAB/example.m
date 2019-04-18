  % Before first run, compile
  mex radon.cpp
  mex radon_adjoint.cpp
  mex minimize_in_EDGE_variable.cpp
  mex minimize_in_IMAGE_variable.cpp

  %Parameters
  IterAlt = 1;
  IterImage = 100;
  IterEdge = 10;
  beta = 0.0002;
  alpha = 1000;
  epsilon = 0.0001;

  % Data
  f = phantom(512);
  g      = Anm(f);
  c      = 0.05;
  n      = randn(size(g));
  gDelta = g + c*max(max(abs(g)))*n;
  
  % Mumford-Shah minimization     
  tic
  [I,E]=MumfordShah_reg(gDelta,alpha,beta,0.00001);
  toc

  subplot(1,3,1)
  imagesc(gDelta)
  subplot(1,3,2)
  imagesc(I)
  subplot(1,3,3)
  imagesc(E)
  colormap gray
