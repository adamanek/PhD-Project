bin = zeros(8,4); %empty array to save ratios
plo = zeros(1,4); %empty array for area-normalised pops
S = []; %empty matrix to store snapshot data
F=[];
near_all = [];mid_all = [];snd_all =[];far_all=[];
%%This section defines data to be used through out the code (makes it more
%%simple than manually changing everything in the code)

data = hydrationdppclaurPfixed;  %Data to be analysed
snap = lenPfix;   %snapshot data (number of lipids per snapshot)
hydration = 15;  %the hydrations to loop over

%Defining region of interest thresholds (radii) in Angstrom    
 R=[0,10,14,18,25]; %[innermost radius, near radius, mid radius, far radius, outermost edge]
  
%area of regions calculations
area0=4/3*pi()*(R(1)^3);
areaX=4/3*pi()*(R(2)^3);
area1=4/3*pi()*(R(3)^3);
area2=4/3*pi()*(R(4)^3);
area3=4/3*pi()*(R(5)^3);
area_near=areaX-area0;  %near region area
area_X=area1-areaX;
area_mid=area2-area1;   %mid region area
area_far=area3-area2;   %far region area

for i=1:1:hydration    
    %initialising region counts
    near=0;
    X_c=0;
    mid=0;
    far=0;
    
    %finds out the number of lipids with hydration of 'i' per specified
    %region
    for k=1:length(data(:,1))
        if data(k,2)==i && data(k,1)>R(1) && data(k,1)<R(2)
            near=near+1;
            near_all=[near_all, data(k,2)];
        end
        if data(k,2)==i && data(k,1)>R(2) && data(k,1)<=R(3)
            X_c=X_c+1;
            snd_all=[snd_all, data(k,2)];
        end
        if data(k,2)==i && data(k,1)>R(3) && data(k,1)<=R(4)
            mid=mid+1;
            mid_all = [mid_all, data(k,2)];
        end
        if data(k,2)==i && data(k,1)>R(4) && data(k,1)<=R(5)
            far=far+1;
            far_all = [far_all, data(k,2)];
        end
    end
    l=int64(i);
    M(l,1)=near;
    M(l,2)=X_c;
    M(l,3)=mid;
    M(l,4)=far;
    %saves hydrations
    bin(l,1)=i;
    
    %calculates the area-normalised populations
    plo(l,1)=(near)/area_near;
    plo(l,2)=(X_c)/area_X;
    plo(l,3)=(mid)/area_mid;
    plo(l,4)=(far)/area_far;
    
    %%This section calculates area-normalised hydrations for each region and each snapshot 
    
    %sets up the intial loop values (so it only does the first snapshot)
    t=1;
    tf=snap(1,1);
    
    for m=1:1:(length(snap)-1)
        %initialisng region counts
        near=0;
        X_c=0;
        mid=0;
        far=0;
        
        %loops through data per snapshot (defined in 'data') and then saves
        %it into a 3D matrix (each page corresponds to a different
        %region). 
        for j=t:1:tf
            if data(j,2)==i && data(j,1)>R(1) && data(j,1)<R(2)
                near=near+1;
            end
            if data(j,2)==i && data(j,1)>R(2) && data(j,1)<=R(3)
                X_c=X_c+1;
            end
            if data(j,2)==i && data(j,1)>R(3) && data(j,1)<=R(4)
                mid=mid+1;
            end
            if data(j,2)==i && data(j,1)>R(4) && data(j,1)<=R(5)
                far=far+1;
            end
        end
        S(m,l,1)=(near/area_near);
        S(m,l,2)=(X_c/area_X);
        S(m,l,3)=(mid/area_mid);
        S(m,l,4)=(far/area_far);
        
        %changes the values so that the next snap is analysed
        t=t+snap(m,1);
        tf=tf+snap(m+1,1);
    end
end

%%This bit calculates the probabiliy of lipids having specified hydration
%%(divides the region normalised lipid counts for each hydration by the
%%total number of lipids in the region) and for st error

SN=[]; %empty matrix for normalised snapshot data
lik=[];  %empty for calculating likelihood of normalised lipids
for k=1:1:4
    for j=1:1:hydration
        lik(j,k)=plo(j,k)/sum(plo(:,k));
        for m=1:1:(length(snap)-1)
            SN(m,j,k)=S(m,j,k)/sum(S(m,:,k));
            if sum(S(m,:,k))==0 %to prevent division by 0
                SN(m,j,k)=0;
            end
        end
    end
end

%standard error calc
SE=[];
for i=1:1:hydration
    SE(i,1)=std(SN(:,i,1))/sqrt(length(SN(:,i,1)));
    SE(i,2)=std(SN(:,i,2))/sqrt(length(SN(:,i,2)));
    SE(i,3)=std(SN(:,i,3))/sqrt(length(SN(:,i,3)));
    SE(i,4)=std(SN(:,i,4))/sqrt(length(SN(:,i,4)));
end

% plotting
figure()
errorbar(lik,SE);
xlabel('Hydration');
ylabel('Probability');
legend('0-10 A','10-14 A','14-18 A','18-25 A');

    

