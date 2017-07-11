function write_file(A,B)
    
    fid = fopen(B,'w');
    
    for i=1:512
        for j=1:512
            fprintf(fid,'%f ',A(i,j));
        end
        fprintf(fid,'\n');
    end
    
    fclose(fid);

end