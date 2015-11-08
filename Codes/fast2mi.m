function MIs=fast2mi(dat1,dat2,h)

[ro1,L]=size(dat1);
[ro2,L]=size(dat2);
MIs = zeros(ro1,ro2);
h_square = h^2;

for i=1:L
    tmp1=dat1-repmat(dat1(:,i),1,L);
    tmp1=exp(-(tmp1.^2)/(2*h_square));
    tmp1v = sum(tmp1,2);
    
    tmp2=dat2-repmat(dat2(:,i),1,L);
    tmp2=exp(-(tmp2.^2)/(2*h_square));
    tmp2v = sum(tmp2,2);
    
%     tmp1v=1./tmp1v;
%     tmp2v=1./tmp2v;
%    tmpv=tmp1v*tmp2v';
    tmp12=tmp1*tmp2';
 %   tmp12=tmp12.*tmpv;
    
    for j=1:ro1
         for k=1:ro2
            tmp12(j,k) = tmp12(j,k)/(tmp1v(j)*tmp2v(k));
            %tmp12(j,k) = tmp12(j,k)/tmp2v(k);
         end
    end
%      for k=1:ro2
%          tmp12(:,k) = tmp12(:,k)/tmp2v(k);
%      end
%     
    MIs = MIs + log(tmp12);
    
end
MIs = MIs/L + log(L);