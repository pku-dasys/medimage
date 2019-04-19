function main()
  %prof_edge_image('shepp_logan_f_9_alpha_0.001_beta_1','shepp_logan_v_9_alpha_0.001_beta_1');
  %prof_edge_image('shepp_logan_f_9_alpha_0.01_beta_1','shepp_logan_v_9_alpha_0.01_beta_1');
  %prof_edge_image('shepp_logan_f_9_alpha_0.1_beta_1','shepp_logan_v_9_alpha_0.1_beta_1');
  %prof_edge_image('shepp_logan_f_9_alpha_1_beta_1','shepp_logan_v_9_alpha_1_beta_1');
  %prof_edge_image('shepp_logan_f_9_alpha_2_beta_1','shepp_logan_v_9_alpha_2_beta_1');
  %prof_edge_image('shepp_logan_f_9_alpha_4_beta_1','shepp_logan_v_9_alpha_4_beta_1');
  %prof_edge_image('shepp_logan_f_9_alpha_8_beta_1','shepp_logan_v_9_alpha_8_beta_1');
  %prof_edge_image('shepp_logan_f_9_alpha_10_beta_1','shepp_logan_v_9_alpha_10_beta_1');
  %prof_edge_image('SCAN3_f_9_alpha_8_beta_4','SCAN3_v_9_alpha_8_beta_4');
  prof_edge_image('SCAN3_f_0_alpha_1_beta_1_detz_2.05239','SCAN3_f_0_alpha_1_beta_1_detz_2.05239');
  %prof_edge_image('SCAN5_f_9_alpha_8_beta_2','SCAN5_v_9_alpha_8_beta_2');
end

function prof_edge_image(image_filename,edge_filename)
  k = strfind(image_filename,'SCAN');
  if length(k) == 0
    prof_image_human(image_filename);
    prof_image_human(edge_filename);
  else
    prof_image_phantom(edge_filename);
    %prof_data_phantom(edge_filename);
    prof_image_phantom(image_filename);
    %prof_data_phantom(image_filename);
  end
end

function n = prof_image_phantom(filename)
  recon_images = load (filename);
  xyz = size(recon_images);
  disp(xyz);
  mkdir(strcat(filename, '_png'));
  Z = xyz(1)/874;
  for z = 11 : Z-10
    disp(z);
    recon_image = recon_images((z-1)*874+1:(z)*874,:);
    recon_image = recon_image/max(max(recon_image));
    recon_image = recon_image';%'
    recon_image = recon_image(:,874:-1:1);
    recon_image = recon_image(182:693,182:693);
    png_name = strcat(filename,'_png/',int2str(z-10),'.png');
    imwrite(recon_image,png_name);
    end
end

function n = prof_data_phantom(filename)
  recon_images = load (filename);
  xyz = size(recon_images);
  disp(xyz);
  mkdir(strcat(filename, '_data'));
  Z = xyz(1)/874;
  for z = 11 : Z-11
    disp(z);
    recon_image = recon_images((z-1)*874+1:(z)*874,:);
    recon_image = recon_image';%'
    recon_image = recon_image(:,874:-1:1);
    recon_image = recon_image(182:693,182:693);
    data_name = strcat(filename,'_data/',int2str(z-10),'_data');
    %imwrite(recon_image,png_name);
    dlmwrite(data_name,recon_image);
    end
end




function n = prof_image_human(filename)
  recon_images = load (filename);
  xyz = size(recon_images);
  disp(xyz);
  mkdir(strcat(filename, '_png'));
  Z = xyz(1)/512;
  for z = 1 : Z-1
    %disp(z);
    recon_image = recon_images((z-1)*512+1:(z)*512,:);
    recon_image = recon_image/max(max(recon_image));
    recon_image = recon_image';%'
    recon_image = recon_image(:,512:-1:1);
    png_name = strcat(filename,'_png/',int2str(z),'.png');
    imwrite(recon_image,png_name);
    end
end
