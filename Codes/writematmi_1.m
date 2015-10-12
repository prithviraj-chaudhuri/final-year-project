function writematmi(filename,matrix)
% data format of new c based MI program.
% matrix data eq row==number of row /col in MI mat 

fd=fopen(filename,'w');
%fd1=fopen(file1,'w');
[row,col]=size(matrix);
buff=[];
 wrtbuff=zeros(1,col);
 
for ii=1:1:row
    wrtbuff(1:col)=matrix(ii,:);
    tbuff=sprintf('%d ',1);
    
    tbuff=sprintf(repmat('%.5f,',1,col),wrtbuff);
    tbuff=[tbuff,sprintf('\n')];
    %buff=[buff,tbuff];
    fwrite(fd,tbuff);
    clear tbuff;
end

   % fwrite(fd,buff);
   
clear buff;
fclose(fd);