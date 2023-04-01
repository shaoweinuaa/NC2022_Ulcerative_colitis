clear all
clc
close all

stain_file_name ={'141Pr_CD14.ome.tiff','142Nd_FOXP3.ome.tiff','143Nd_CD16.ome.tiff','144Nd_CD69.ome.tiff','145Nd_CD4.ome.tiff','146Nd_CD8.ome.tiff',...
'147Sm_Collagen1.ome.tiff','148Nd_caspase3.ome.tiff','149Sm_CD31.ome.tiff','150Nd_E_cadherin.ome.tiff','151Eu_B7H4.ome.tiff','152Sm_VISTA.ome.tiff',...
  '153Eu_CD7.ome.tiff','154Sm_CD163.ome.tiff','155Gd_CD103.ome.tiff','156Gd_PDL1.ome.tiff', '158Gd_LAG3.ome.tiff','159Tb_CD68.ome.tiff','160Gd_CD11B.ome.tiff',...
  '161Dy_CD20.ome.tiff','162Dy_CD11C.ome.tiff','163Dy_CD15.ome.tiff','164Dy_Granzyme_B.ome.tiff','165Ho_PD1.ome.tiff','166Er_KI67.ome.tiff', '167Er_GATA3.ome.tiff',...
  '168Er_HLADR.ome.tiff','169Tm_CD45RA.ome.tiff','170Er_CD3.ome.tiff', '171Yb_TNFa.ome.tiff','172Yb_IL1b.ome.tiff', '173Yb_CD45RO.ome.tiff',...
  '174Yb_CD57.ome.tiff','175Lu_CD25.ome.tiff','176Yb_PAN-Keratin.ome.tiff','194Pt_aSMA.ome.tiff','198Pt_VIMENTIN.ome.tiff','89Y_CD45.ome.tiff','115In_IL6.ome.tiff'}



input='.\\IHC\\';  %%%%%% 16-bits .tiff image files
denoiseRoot='.\\Denoise\\';  %%%% denoised images

inputDir=dir(input);
roiNumber=0; %%%%%%%total ROI number
sampleNumber=0; %%%%%%%%%%%total sample number
sampleROI={};   %%%%%% the index of ROIs for each sample 
patient_name={};
feature={};
cellNumber=1;
cellROI={};
for i=3:length(inputDir)
    sampleName =inputDir(i).name; 
    sampleNumber=sampleNumber+1;
    sampleROI{sampleNumber}=[];
    saveDenoiseSampleRoot=strcat(denoiseRoot,sampleName);
    if ~exist(saveDenoiseSampleRoot)
        mkdir(saveDenoiseSampleRoot);
    end
    
    inputSampleRoot=strcat(input,'\\',sampleName);
    roiDir=dir(inputSampleRoot);
    for j=3:length(roiDir)
        dataOrigin=[];
        dataDenoise=[];
        roiNumber=roiNumber+1;
        sampleROI{sampleNumber}=[sampleROI{sampleNumber},roiNumber];
        roiName=roiDir(j).name
        saveDenoiseROIRoot=strcat(saveDenoiseSampleRoot,'\\',roiName,'\\');
        if ~exist(saveDenoiseROIRoot)
               mkdir(saveDenoiseROIRoot);
        end 
        inputROIRoot=strcat(inputSampleRoot,'\\',roiName);
        inputROIDir=dir(inputROIRoot);
        cell_position={};
        stain_id=0;  
        stain={};
        [neighborMatrix,nuclei_xy]=getNucleiCoordinate(inputROIRoot);
        output=zeros(size(nuclei_xy,1),39);  
        inputAllImages=[];
        for k=3:length(inputROIDir)
            flag=0;
            tiffFileName=inputROIDir(k).name
            if sum(contains(stain_file_name,tiffFileName))==1
               stain_id=stain_id+1;
               inputImageFileRoot=strcat(inputROIRoot,'\\',tiffFileName);
               im = tiffread2(inputImageFileRoot);
               dataOrigin=im.data;
               dataOrigin(dataOrigin>255)=255;
               data_thresh=prctile(dataOrigin(find(dataOrigin~=0)),90);
               dataOrigin=uint8(dataOrigin*ceil(255.0/data_thresh));
               inputAllImages(stain_id,:,:)=dataOrigin;
               index=str2num(tiffFileName(1:3))
               outputImage=spillImage(inputAllImages,dataOrigin,index);
               mag= medfilt2(dataOrigin,[3,3]);
               for step=0.01:0.02:1
                   temp_data=imadjust(mag,[0,step],[0,1]);
                   ratio_temp=length(find(temp_data==max(max(temp_data))))/(size(temp_data,1)*size(temp_data,2));
                   if ratio_temp<0.5
                      dataDenoise=temp_data;
                      break
                      end
                 
                      if step==0.99  && ratio_temp>0.5
                         dataDenoise=mag;
                      end
                  end
                  
                  t=tabulate(reshape(dataDenoise,[1,size(dataDenoise,1)*size(dataDenoise,2)]));
                  sum_t=0;
                  temp_c=0;
                  for c=size(t,1):-1:1
                      if c==size(t,1) || (sum_t+t(c,3)/100)<0.5
                         dataDenoise(find(dataDenoise==t(c,1)))=255;
                         sum_t=sum_t+t(c,3)/100;
                         temp_c=t(c,1);
                      end
                      if sum_t>0.5
                         break
                      end
                      
                  end
                  
                  dataDenoise(find(dataDenoise<temp_c))=0;
                  s_cell = regionprops(imbinarize(dataDenoise), 'Area','Centroid','PixelIdxList');
                  cell_number=1;
                  temp_cell=[];
                  for ii=1:length(s_cell)
                      if s_cell(ii).Area>=2
                        [idx,B]=knnsearch(nuclei_xy,s_cell(ii).Centroid,'k',1);
                         x=floor(s_cell(ii).PixelIdxList/size(dataDenoise,1))+1;
                         y=s_cell(ii).PixelIdxList-(x-1)*size(dataDenoise,1);
                         coordinate= [y,x];
                         if B>15 && length(intersect(coordinate,nuclei_xy))==0
                            dataDenoise(s_cell(ii).PixelIdxList)=0;  
                         else        
                           
                          [nuclei_id,distance]=knnsearch(coordinate,nuclei_xy);
                          index=find(distance<15);                   
                          for kkk=1:length(index)
                             tempZ=repmat(nuclei_xy(index(kkk),:),[size(coordinate,1),1])-coordinate;
                             newDistance=sqrt(sum(tempZ.*tempZ,2));
                             idx_new=find(newDistance<20);
                             PixelIdxList=s_cell(ii).PixelIdxList;
                             output(index(kkk),stain_id)= output(index(kkk),stain_id)+mean(dataOrigin(PixelIdxList(idx_new)));
                          end

                         end
                      else
                         dataDenoise(s_cell(ii).PixelIdxList)=0;
                      end
                    
                  end
                  cell_position{stain_id}=temp_cell;
                  stain{stain_id}=tiffFileName;  
                  imwrite(dataDenoise,strcat(saveDenoiseROIRoot,tiffFileName));
                  end

           end
           feature{roiNumber}=[nuclei_xy,output];
           
       end
              
          
          
          
end

     
      




     
   
     
     

     
   
   


     
     
     
     



    
    
    
    
   




