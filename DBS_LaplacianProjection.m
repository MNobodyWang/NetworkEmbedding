% Uses network Laplacian to generate a low-dimensional embedding of a
% network, within context of comparing DBS on and off conditions
% Inputs: roiAMats_DBS_ON: connectivity matrices for DBS ON cases (subjects
% x regions x regions)
% roiAMats_DBS_OFF: same for DBS OFF cases, patients ordered in same order
% as roiAMats_DBS_OFF
% testSubj: which subject is currently being tested in leave-one-out and
% should be left out of creating the low-dimensional embedding (should be
% between 1 and subjects count)
% Created on 20211123 by Max B Wang

function [onProj,offProj]=DBS_LaplacianProjection(roiAMats_DBS_ON,roiAMats_DBS_OFF,testSubj)

% Remove test subject and calculate low dimensional Laplacian embedding

pruned_OFF=roiAMats_DBS_OFF;
pruned_OFF(testSubj,:,:)=[];

ave_W=squeeze(mean(pruned_OFF,1));
aveD_OFF=sum(ave_W);
aveL_Off=diag(aveD_OFF)-ave_W;
[meanVecs_Basis,meanVals_Basis]=eig(aveL_Off);

numSubjs=size(roiAMats_DBS_ON,1);
numRegions=size(roiAMats_DBS_ON,2);

onProj=zeros(numSubjs,numRegions-1);
offProj=zeros(numSubjs,numRegions-1);

for subj=1:numSubjs
    W_on=squeeze(roiAMats_DBS_ON(subj,:,:));
    D_on=sum(W_on);
    L_on=diag(D_on)-W_on;
    
    for vec=2:numRegions
        projMVec=L_on*meanVecs_Basis(:,vec);
        onProj(subj,vec-1)=(projMVec.'*meanVecs_Basis(:,vec))/meanVals_Basis(vec,vec);
    end
    
    W_off=squeeze(roiAMats_DBS_OFF(subj,:,:));
    D_off=sum(W_off);
    L_off=diag(D_off)-W_off;
    
    for vec=2:numRegions
        projMVec=L_off*meanVecs_Basis(:,vec);
        offProj(subj,vec-1)=(projMVec.'*meanVecs_Basis(:,vec))/meanVals_Basis(vec,vec);
    end
end