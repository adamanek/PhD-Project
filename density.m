%close all;
%clearvars -except headgroupangles20Adppc headgroupangles20A
clear A B C
A=[];
x=[];
y=[];
data = dppcdens;
data1 = choldens;

for i=0:20:460 %x-axis
    j=i+20;
    for g=0:20:460  %y-axis
        count=0;
        h=g+10;
        c=0;
        for k=1:length(data(:,1))
            if data(k,2)>=g && data(k,2)<h && data(k,1)>=i && data(k,1)<j
                count=count+1;
            end
            if data1(k,2)>=g && data1(k,2)<h && data1(k,1)>=i && data1(k,1)<j
                c=c+1;
                
            end
        end
        l=int64((i+20)/20);
        f=int64((g+20)/20);
        A(f,l)=count;
        B(f,l)=c;
        x(l)=i;
        y(f)=g;
 
    end
end

for k=1:length(A(:,1))
    for j=1:length(A(1,:))
        if A(k,j)==0 && B(k,j)==0 || B(k,j)==0
            C(k,j)=0;
        else
            C(k,j)=A(k,j)/B(k,j);
        end
    end
end
C=A-B
figure
pcolor(x,y,C);
colormap(jet);
xlabel('Distance/Angstrom');
ylabel('Headgroup angle/Degrees');

% den=[];
% den1=[];
% RA=[];
% for m=0:0.2:19.8
%     n=m+0.2;
%     areaM = pi()*m^2;
%     areaN = pi()*n^2;
%     area = areaN-areaM;
%     count=0;
%     for k=1:length(data(:,1))
%         if data(k,1)>=m && data(k,1)<n
%                 count=count+1;
%         end
%     end
%     l=int64(m*5+1);
%     den(l,1)=m;
%     den(l,2)=count/area;
%     count=0;
%     for k=1:length(data1(:,1))
%         if data1(k,1)>=m && data1(k,1)<n
%                 count=count+1;
%         end
%     end
%     den1(l,1)=m;
%     den1(l,2)=count/area;
%     RA(l,1)=den1(l,2)/den(l,2);
% end
% figure()
% plot(den(:,1),RA(:,1));
% xlabel('Distance from Prodan / Angstrom');
% ylabel('Density ratio')