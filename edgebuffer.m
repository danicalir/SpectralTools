function [edge] = edgebuffer(Z,edgeelements)
for iy=1:size(Z,1)
    znan=isnan(Z(iy,:));
    zadd1=[znan;znan];
    for iedgex=1:edgeelements
        zadd1(iedgex.*2+1,:)=[0 zadd1(iedgex*2-1,1:end-1)]; %add the next column with 1 below original NaN position
        zadd1(iedgex.*2+2,:)=[zadd1(iedgex*2,2:end) 0]; %add the next column with 1 above original NaN position
    end
    clipx(iy,:)=sum(zadd1);
    clipx(iy,logical(clipx(iy,:)))=NaN;
end
for ix=1:size(Z,2)
    znan=isnan(Z(:,ix));
    zadd1=[znan znan];
    for iedgey=1:edgeelements
        zadd1(:,iedgey.*2+1)=[0;zadd1(1:end-1,iedgey*2-1)]; %add the next column with 1 below original NaN position
        zadd1(:,iedgey.*2+2)=[zadd1(2:end,iedgey*2);0]; %add the next column with 1 above original NaN position
    end
    clipy(:,ix)=sum(zadd1,2);
    clipy(logical(clipy(:,ix)),ix)=NaN;
end

edge=clipy+clipx;

Zclipped=Z+edge;