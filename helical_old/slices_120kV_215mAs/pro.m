for k = 1:50
    slice_name = strcat('f_',int2str(k));
    f=load(slice_name);
    xyz = size(f);
    Z = xyz(1)/512;
    mkdir(strcat('f_',int2str(k),'_png'));
    for z = 1:Z
        slice = f((z-1)*512+1:z*512,:);
        slice=slice/max(max(slice));
        png_name = strcat('f_',int2str(k),'_png\',int2str(z),'.png');
        imwrite(slice,png_name);
    end
end
