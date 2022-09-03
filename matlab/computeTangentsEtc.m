function [lengths,tangents,tangents_cmplx,normals,normals_cmplx,allTangents,next,curveSize] = computeTangentsEtc(boundaryCurve,openCurve)
allTangents = [];
for k=1:numel(boundaryCurve)
    next{k} = circshift(boundaryCurve{k},-1,1);
    tangents{k} = next{k}-boundaryCurve{k};  
   
    if (openCurve(k))
        next{k}(end,:) = 0; %does not exist
        tangents{k}(end,:)=[];
        curveSize(k) = size(boundaryCurve{k},1)-1;
    else
        curveSize(k) = size(boundaryCurve{k},1);
    end
    %each finite element is definted as the edge between i and (i+1)
    
    lengths{k} = sqrt(sum(tangents{k}.^2,2));
    tangents{k} = tangents{k}./lengths{k};
    normals{k} = [tangents{k}(:,2), -tangents{k}(:,1)];
    tangents_cmplx{k} = tangents{k}(:,1)+1i*tangents{k}(:,2);
    normals_cmplx{k} = normals{k}(:,1)+1i*normals{k}(:,2);
    allTangents = [allTangents;  tangents_cmplx{k}];
end
end

