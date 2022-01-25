fid = fopen('ssim.txt','w');
b=load('phantom\std_f.dat');
for i=1:50
  a=load(strcat('exper\latest\f_iter_',num2str(i),'.dat'));
  [mmism map] = ssim_index(a,b);
  fprintf(fid,'%f\n',mmism);
endfor
fclose(fid);
