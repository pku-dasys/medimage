for i=0:1:1023
    img_data=['IMG_',num2str(i)]
    img_data=img_data/5
    jpg_file=['Fig_',num2str(i),'.jpg']
    imwrite(img_data, jpg_file, 'jpg')
end