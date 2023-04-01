function  []=writeCSVFile(TT,roiName, neighborMatrix)

columnNumber =max(sum(neighborMatrix,2))+42;
outputMatrix=zeros(size(neighborMatrix,1),columnNumber);
outputMatrix(:,1)=(1:size(neighborMatrix,1))';
outputMatrix(:,2:42)=TT;
outputMatrix(:,43)=sum(neighborMatrix,2);
for i=1:size(neighborMatrix,1)
   tempRow=neighborMatrix(i,:);
   neighborIndex=find(tempRow==1);
   outputMatrix(i,44:44+length(neighborIndex)-1)=neighborIndex;

end
title={'ROI','CellId','X_position','Y_position','115In_IL6','141Pr_CD14','142Nd_FOXP3',	'143Nd_CD16','144Nd_CD69','145Nd_CD4','146Nd_CD8','147Sm_Collagen1','148Nd_caspase3','149Sm_CD31','150Nd_E_cadherin','151Eu_B7H4','152Sm_VISTA','153Eu_CD7','154Sm_CD163','155Gd_CD103','156Gd_PDL1','158Gd_LAG3','159Tb_CD68','160Gd_CD11B','161Dy_CD20','162Dy_CD11C','163Dy_CD15','164Dy_Granzyme_B','165Ho_PD1','166Er_KI67','167Er_GATA3','168Er_HLADR','169Tm_CD45RA','170Er_CD3','171Yb_TNFa','172Yb_IL1b','173Yb_CD45RO','174Yb_CD57','175Lu_CD25','176Yb_PAN-Keratin','194Pt_aSMA','198Pt_VIMENTIN','89Y_CD45','Number_Neighbors'};
for i=1:max(sum(neighborMatrix,2))
    title{i+44}=strcat('neighbour_',num2str(i));
end
csvFileName=strcat('.\\DiscriminantGraph\\test\\',roiName,'.csv');
fid=fopen(csvFileName,'w');
for i=1:size(outputMatrix,1)+1
    i
    if i==1
      for j=1:length(title)
          if j<length(title)
            fprintf(fid,'%s,',title{j});
          else
            fprintf(fid,'%s\n',title{j});
          end
      end
      
    else
       for j=1:length(title)
          if j==1
              fprintf(fid,'%s,',roiName);
          elseif j>1 && j<length(title)
              fprintf(fid,'%.3f,',outputMatrix(i-1,j-1));
          else
              fprintf(fid,'%.3f\n',outputMatrix(i-1,j-1));
          end
       end
        
      
    end
end  
fclose(fid)
end