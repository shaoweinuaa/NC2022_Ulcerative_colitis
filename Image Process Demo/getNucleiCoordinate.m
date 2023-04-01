function [neighborMatrix,nuclei_xy,xSize,ySize]=getNucleiCoordinate(nucleiRoot)
 nuclei_xy=[];
 radius_xy=[];
 nucleiSampleTiffFile=strcat(nucleiRoot,'\\191Ir_DNA1.ome.tiff');
 im = tiffread2(nucleiSampleTiffFile);
 data=im.data;
 data(data>255)=255;
 data_thresh=prctile(data(find(data~=0)),90);
 data=uint8(data*ceil(255.0/data_thresh));
 mag= medfilt2(data,[3,3]);
 xSize=size(mag,1);
 ySize=size(mag,2);
%  figure
%  imshow(mag)
 for step=0.01:0.02:1
     temp_data=imadjust(mag,[0,step],[0,1]);
     ratio_temp=length(find(temp_data==max(max(temp_data))))/(size(temp_data,1)*size(temp_data,2));
     if ratio_temp<0.2
         data2=temp_data;
         break
     end
                 
     if step==0.99  && ratio_temp>0.075
        data2=mag;
        end
  end
                  
  t=tabulate(reshape(data2,[1,size(data2,1)*size(data2,2)]));
  sum_t=0;
  temp_c=0;
  for c=size(t,1):-1:1
      if c==size(t,1) || (sum_t+t(c,3)/100)<0.1
         data2(find(data2==t(c,1)))=255;
                   
         sum_t=sum_t+t(c,3)/100;
         temp_c=t(c,1);
      end
      if sum_t>0.05
         break
      end
                      
  end
  data2(find(data2<temp_c))=0;
%   figure
%   imshow(data2)
  s_cell = regionprops(imbinarize(data2), 'Area','Centroid','PixelIdxList','MajorAxisLength','MinorAxisLength');
  for ii=1:length(s_cell)
      if s_cell(ii).Area>5 && s_cell(ii).Area<1000 
        nuclei_xy=[nuclei_xy;s_cell(ii).Centroid];
        radius_xy=[radius_xy;s_cell(ii).MajorAxisLength/2];
      else  
        data2(s_cell(ii).PixelIdxList)=0; 
      end
  end
  A=zeros(size(nuclei_xy,1),size(nuclei_xy,1));
  k=1;
  for i=1:size(nuclei_xy,1)
      for j=i+1:size(nuclei_xy,1)
        dist(i,j)=norm(nuclei_xy(i,:)-nuclei_xy(j,:),2)-radius_xy(i)-radius_xy(j);
        if dist(i,j)<20
           small_distance_xy(k,:)=[i,j];
           A(i,j)=1;
           k=k+1;
        
        end
      end
  end
  neighborMatrix=A+A';


