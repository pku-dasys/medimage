filename='SCAN1_f_1.5mm_';
for k = 0:2
    slice_name = strcat(filename,int2str(k));
    f=load(slice_name);
    xyz = size(f);
    Z = xyz(1)/512;
    mkdir(strcat(filename,int2str(k),'_png'));
    for z = 1:Z
        slice = f((z-1)*512+1:z*512,:);
        slice=slice/max(max(slice));
        png_name = strcat(filename,int2str(k),'_png\',filename,int2str(k),'_',int2str(z),'.png');
        imwrite(slice,png_name);
    end
end
