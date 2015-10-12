function MIs  = FastPairMI(data,h)
% data : the input data, rows correspond to genes(variables)
%        columns correspond to arrays (samples)  
% h    : the std of the Gaussian kernel for density estimation 

MIs = zeros(size(data,1));
h_square = h^2;
L = size(data,2); % excluding the class label
for i=1:L
    
    tmp = data - repmat(data(:,i),1,L);
    tmp = exp(-(tmp.^2)/(2*h_square));
    tmp1 = sum(tmp,2);

    tmp2 = tmp*tmp';
    for j=1:size(tmp2,1)
        tmp2(j,:) = tmp2(j,:)./tmp1(j);
        tmp2(:,j) = tmp2(:,j)./tmp1(j);
    end
    MIs = MIs + log(tmp2);
    clear tmp2

% %   The following commented line does the same job as lines 16~22
%     MIs = MIs + log((tmp*tmp')./(tmp1*tmp1'));
end
MIs = MIs/L + log(L);

%save(filenam,'-ascii','-tabs','MIs')
% M=size(MIs,1);
% N=size(MIs,2);
% fid=fopen('MIs.txt','w');
% for i=1:M'
%     for j=1:N,
%         fprintf(fid,'%.6\t',MIs(i,j));
%     end
%     fprintf(fid,'\n');
% end
% fclose(fid);

return