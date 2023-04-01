function outputImage=spillImage(inputAllImages,originalImage,index)

metalNumber=[144,146,148,162,163,164,167,168,172,173,174];
emmitIndex=[3,5,7,20,21,22,25,26,30,31,32];
spillValue=[2.7,4.3,2.8,3.3,5.1,2.3,2.1,3.2,4.9,5.1,5.4]/100;

if length(index)>0 &&  length(find(metalNumber==index))>0 
   indicator = find(metalNumber==index);
   emmitImage=squeeze(inputAllImages(indicator,:,:));
   outputImage=uint8(double(ceil(originalImage))-emmitImage*spillValue(indicator));
   outputImage(outputImage<0)=0;
   


else
   outputImage=originalImage;
end