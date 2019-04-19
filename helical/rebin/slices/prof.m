function main()
  
  recon_folder = './SHEPP_LOGAN/slices_1mm/';
  files = dir(strcat(recon_folder,'shepp_gauss_1*'));
  L = length(files)
  for i = 1:L
    files(i);
    if files(i).isdir==0
      disp(files(i).name);
      if strfind(files(i).name,'SCAN')
        prof_image_test(strcat(recon_folder,files(i).name),874,0);
      elseif strfind(files(i).name,'NLST')
        prof_image_test(strcat(recon_folder,files(i).name),771,0);
      else
        prof_image_test(strcat(recon_folder,files(i).name),512,0);
      end
    end
  end
end

%{
function n = prof_image_test(filename,SIZE,offset)
  recon_images = load (filename);
  xyz = size(recon_images);
  disp(xyz);
  mkdir(strcat(filename, '_png'));
  Z = xyz(1)/SIZE;
  for z = offset+1 : Z-offset
    %disp(z);
    recon_image = recon_images((z-1)*SIZE+1:(z)*SIZE,:);
    gray_name = strcat(filename,'_png/',int2str(z-offset),'.txt');
    dlmwrite(gray_name,recon_image);
    recon_image = recon_image/max(max(recon_image));
    recon_image = recon_image';%'
    recon_image= recon_image(:,SIZE:-1:1);
    png_name = strcat(filename,'_png/',int2str(z-offset),'.png');
    imwrite(mat2gray(recon_image),png_name);
  end
end
%}

function n = prof_image_test(filename,SIZE,offset)
  recon_images = load (filename);
  xyz = size(recon_images);
  disp(xyz);
  mkdir(strcat(filename, '_png'));
  Z = xyz(1)/SIZE;
  for z = offset+1 : Z-offset
    %disp(z);
    recon_image = recon_images((z-1)*SIZE+1:(z)*SIZE,:);
    gray_name = strcat(filename,'_png/',int2str(z-offset),'.txt');
    dlmwrite(gray_name,recon_image);
    recon_image = recon_image/max(max(recon_image));
    recon_image = recon_image';%'
    recon_image = recon_image(:,SIZE:-1:1);
    png_name = strcat(filename,'_png/',int2str(z-offset),'.png');
    imwrite(recon_image,png_name);
    end
end


