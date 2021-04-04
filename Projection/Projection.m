[modelPoints, number] = ReadObj();
originalImage = imread('foxiang.PNG');
[projectionMatrix, rawPicturePoints, rawModelPoints] = GetProjectionMatrix(1,6);

picturePoints = projectionMatrix * modelPoints;

picturePoints = NormalizeNoChangeZ(picturePoints);

maxx=ceil(max(picturePoints(1,:)))+2;
maxy=ceil(max(picturePoints(2,:)))+2;
used=zeros(maxx,maxy,2);
index=[];

for i = 1:length(picturePoints)
    x=picturePoints(1,i);
    y=picturePoints(2,i);
    z=picturePoints(3,i);
    x=ceil(x);
    y=ceil(y);
    if x<1
        x=1;
    else if x>maxx-1
            x=maxx-1;
        end
    end
    if y<1
        y=1;
    else if y>maxy-1
            y=maxy-1;
        end
    end
    
    if used(x,y,1)<z && used(x,y,1)~=0
        continue;
    end
    
    
    used(x,y,1)=z;
    used(x,y,2)=i;
end


filteredpicturePoint=[];rawusedIndex=used(:,:,2);
size1=size(rawusedIndex,1);
size2=size(rawusedIndex,2);
rawusedIndex=reshape(rawusedIndex,size1*size2,1);
rawusedIndex=rawusedIndex';
rawusedIndex(:,all(rawusedIndex==0,1))=[];
rawusedIndex=rawusedIndex';
for i=1:length(rawusedIndex)
    filteredpicturePoint(i,:)=picturePoints(:,rawusedIndex(i));
end
filteredpicturePoint=filteredpicturePoint';
DrawPicture(filteredpicturePoint,originalImage,projectionMatrix,rawModelPoints,rawPicturePoints);