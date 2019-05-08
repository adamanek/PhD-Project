%close all;
%clearvars -except headgroupangles20Adppc headgroupangles20A
A=[];
x=[];
y=[];
data=headgroupangles20A;
for i=0:0.5:20
    j=i+0.5;
    c=0;
    for g=0:3.5:140 
        count=0;
        h=g+3.5;
        for k=1:20316
            if data(k,2)>=g && data(k,2)<h && data(k,1)>=i && data(k,1)<j
                count=count+1;
            end
            if data(k,1)>=i && data(k,1)<j
                c=c+1;
                
            end
        end
        l=int64(i*2+1);
        f=int64((g/3.5)+1);
        A(f,l)=count;
        x(l)=i;
        y(f)=g;
 
    end
end
%A=flipud(A);
B=[]; %empty matrix for normalisation
for i=1:length(A)
    for j=1:length(A)
        B(j,i)=A(j,i)/sum(A(:,i));
        if sum(A(:,i))==0
            B(j,i)=0;
        end
    end
end
figure()
pcolor(x,y,A);
xlabel('Distance/Angstrom');
ylabel('Headgroup angle/Degrees');


