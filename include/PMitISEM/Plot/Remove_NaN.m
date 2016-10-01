function Remove_NaN(fname)
    FID = fopen(fname,'r+');
    FID2 = fopen('temp.txt','w+');
    while ~feof(FID)
        s = fgets(FID);
        s = strrep(s, 'NaN s', '--');
        s = strrep(s, 'NaN', '--');
        fprintf(FID2,'%s',s);
    end
    
    fclose(FID);
    fclose(FID2);
    fclose('all');
    delete(fname);
    
    movefile('temp.txt',fname,'f');
end