function miwrite(filename,MIs)

M=size(MIs,1);
N=size(MIs,2);
fid=fopen(filename,'w');
for i=1:M
    for j=1:N
        fprintf(fid,'%0.6f\t',MIs(i,j));
    end
    fprintf(fid,'\n');
end
fclose(fid);