function [mantel_z0,mantel_z,p_value] = mantel_test(mat1,mat2,nperm)

% test if mat1 and mat2 have same size
if(size(mat1) ~=  size(mat2))
    warning("matrices don't have same size");
    mantel_z0 = NaN;
    return 
end

vec1=asvect(mat1);
vec2=asvect(mat2);
mantel_z0 = ((vec1-mean(vec1))/std(vec1))*...
    ((vec2'-mean(vec1))/std(vec2))/(length(vec1)-1);
mantel_z = zeros(nperm, 1);
[m,~]=size(mat1);

n=((m-1)*m)/2;

for i = (1:nperm)
    ni=randperm(n);
    mantel_z(i) =(vec1-mean(vec1))*(vec2(ni)-mean(vec2(ni)))'...
        /(std(vec1)*std(vec2(ni))*(length(vec1)-1));
end

p_value = (sum(mantel_z >= mantel_z0))/(nperm + 1);


end
function vect=asvect(mat)

[m,~]=size(mat);
n=((m-1)*m)/2;
k=1;
vect=zeros(1,n);
for i = (2:m)
    for j = (1:(i-1))
        vect(k)=mat(i,j);
        k=k+1;
    end
end

end