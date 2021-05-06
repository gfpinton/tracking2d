%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GIANMARCO PINTON
% WRITTEN: 2020-04-23
% LAST MODIFIED:  2020-04-26
% LAUNCH 2D QUALITY WEIGHTED ADAPTED TRACKING, EASY WRAPPER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%addpath /celerina/gfp/mfs/dumbmat/
%addpath /celerina/gfp/mfs/gel_2Dtracking_data/

%nohup matlab -nodisplay < YOUR_MATLAB_FILE.m > output.txt &
addpath Sandhya_func_library/


load l115test1_50_1200fact11_rf.mat
i=10;
figure(32);imagesc(dbzero(abs(hilbert(double(squeeze(rf2(:,i,:)))))),[-60 0]);
%rf2=rand(200,500,600);
searchX=15; searchY=20;
windowX=25; windowY=30;
%rf=single(rf2(500:900,:,450:800));
rf=single(rf2(400:900,:,:));
rf=permute(rf,[1 3 2]);
nX=size(rf,1), nY=size(rf,2), nZ=size(rf,3);
idc=zeros(2,nX,nY,nZ-1,'int32');
%for k=1:nZ
%    figure(32);imagesc(dbzero(abs(hilbert(double(squeeze(rf(:,:,k)))))),[-45 0]);
%    drawnow
%end

%% ITERATION 0 
tic
[cc dd] = mex_ncorrv2d_idc_omp2(rf,idc,int32(windowX), int32(windowY), int32(searchX), int32(searchY));
toc

clear rf2; save -V7.3 cc2d_idc_test_mex2_medfilt_rf_workspace2 
    
figure(1), imagesc(squeeze(dd(1,:,:,10))), colorbar
figure(2), imagesc(squeeze(dd(2,:,:,10))), colorbar
figure(3), imagesc(squeeze(cc(:,:,10))), colorbar
figure(4), plot(squeeze(mean(mean(dd(1,round(nX/2)-10:round(nX/2)+10,round(nY/2)-10:round(nY/2)+10,:),2),3))), hold on
figure(4), plot(squeeze(mean(mean(dd(2,round(nX/2)-10:round(nX/2)+10,round(nY/2)-10:round(nY/2)+10,:),2),3)))
figure(4), plot(squeeze(dd(1,round(nX/2),round(nY/2),:)))
figure(4), plot(squeeze(dd(2,round(nX/2),round(nY/2),:))), hold off

%% ITERATION 1 %%
medwin=[3 3 2]; medwin=int32(medwin);
thresh=0.5;

idc=dd; idct=find(cc<thresh); idc(2*idct-1)=0; idc(2*idct)=0;
idc(1,:,:,:)=mex_medfilt3_nonzero_float2(squeeze(idc(1,:,:,:)),medwin(1),medwin(2),medwin(3)); 
idc(2,:,:,:)=mex_medfilt3_nonzero_float2(squeeze(idc(2,:,:,:)),medwin(1),medwin(2),medwin(3)); 
idc2=zeros(2,nX,nY,nZ-1,'int32');
idc2(:,searchX+round(windowX/2):searchX+round(windowX/2)+nX-windowX-searchX*2+1-1,searchY+round(windowY/2):searchY+round(windowY/2)+nY-windowY-searchY*2+1-1,:)=round(idc);

figure(3), imagesc(squeeze(idc(1,:,:,10))), colorbar
figure(4), imagesc(squeeze(idc(2,:,:,10))), colorbar
figure(5), plot(squeeze(mean(mean(idc(1,round(nX/2)-10:round(nX/2)+10,round(nY/2)-10:round(nY/2)+10,:),2),3))), hold on
figure(5), plot(squeeze(mean(mean(idc(2,round(nX/2)-10:round(nX/2)+10,round(nY/2)-10:round(nY/2)+10,:),2),3)))
figure(5), plot(squeeze(idc(1,round(nX/2),round(nY/2),:)))
figure(5), plot(squeeze(idc(2,round(nX/2),round(nY/2),:))), hold off


figure(3), imagesc(squeeze(idc2(1,:,:,10))), colorbar
figure(4), imagesc(squeeze(idc2(2,:,:,10))), colorbar


searchX=3; searchY=4;
tic
[cc dd] = mex_ncorrv2d_idc_omp2(rf,idc2,int32(windowX), int32(windowY), int32(searchX), int32(searchY));
toc

figure(1), imagesc(squeeze(dd(1,:,:,10))), colorbar
figure(2), imagesc(squeeze(dd(2,:,:,10))), colorbar
figure(3), imagesc(squeeze(cc(:,:,10))), colorbar
figure(4), plot(squeeze(mean(mean(dd(1,round(nX/2)-10:round(nX/2)+10,round(nY/2)-10:round(nY/2)+10,:),2),3))), hold on
figure(4), plot(squeeze(mean(mean(dd(2,round(nX/2)-10:round(nX/2)+10,round(nY/2)-10:round(nY/2)+10,:),2),3)))
figure(4), plot(squeeze(dd(1,round(nX/2),round(nY/2),:)))
figure(4), plot(squeeze(dd(2,round(nX/2),round(nY/2),:))), hold off

%% ITERATION 2 %%
medwin=[3 3 2]; medwin=int32(medwin);
thresh=0.5;

idc=dd; idct=find(cc<thresh); idc(2*idct-1)=0; idc(2*idct)=0;
idc(1,:,:,:)=mex_medfilt3_nonzero_float2(squeeze(idc(1,:,:,:)),medwin(1),medwin(2),medwin(3)); 
idc(2,:,:,:)=mex_medfilt3_nonzero_float2(squeeze(idc(2,:,:,:)),medwin(1),medwin(2),medwin(3)); 
idc2=zeros(2,nX,nY,nZ-1,'int32');
idc2(:,searchX+round(windowX/2):searchX+round(windowX/2)+nX-windowX-searchX*2+1-1,searchY+round(windowY/2):searchY+round(windowY/2)+nY-windowY-searchY*2+1-1,:)=round(idc);

figure(3), imagesc(squeeze(idc(1,:,:,10))), colorbar
figure(4), imagesc(squeeze(idc(2,:,:,10))), colorbar
figure(5), plot(squeeze(mean(mean(idc(1,round(nX/2)-10:round(nX/2)+10,round(nY/2)-10:round(nY/2)+10,:),2),3))), hold on
figure(5), plot(squeeze(mean(mean(idc(2,round(nX/2)-10:round(nX/2)+10,round(nY/2)-10:round(nY/2)+10,:),2),3)))
figure(5), plot(squeeze(idc(1,round(nX/2),round(nY/2),:)))
figure(5), plot(squeeze(idc(2,round(nX/2),round(nY/2),:))), hold off
figure(3), imagesc(squeeze(idc2(1,:,:,10))), colorbar
figure(4), imagesc(squeeze(idc2(2,:,:,10))), colorbar


searchX=2; searchY=2;
tic
[cc dd] = mex_ncorrv2d_idc_omp2(rf,idc2,int32(windowX), int32(windowY), int32(searchX), int32(searchY));
toc

figure(1), imagesc(squeeze(dd(1,:,:,10))), colorbar
figure(2), imagesc(squeeze(dd(2,:,:,10))), colorbar
figure(3), imagesc(squeeze(cc(:,:,10))), colorbar
figure(4), plot(squeeze(mean(mean(dd(1,round(nX/2)-10:round(nX/2)+10,round(nY/2)-10:round(nY/2)+10,:),2),3))), hold on
figure(4), plot(squeeze(mean(mean(dd(2,round(nX/2)-10:round(nX/2)+10,round(nY/2)-10:round(nY/2)+10,:),2),3)))
figure(4), plot(squeeze(dd(1,round(nX/2),round(nY/2),:)))
figure(4), plot(squeeze(dd(2,round(nX/2),round(nY/2),:))), hold off


%% FINAL MEDIAN FILTER %%
medwin=[3 3 2]; medwin=int32(medwin);
thresh=0.5;
ddf=dd; ddft=find(cc<thresh); length(ddft)/length(cc(:)), ddf(2*ddft-1)=0; ddf(2*ddft)=0;
ddf(1,:,:,:)=mex_medfilt3_nonzero_float2(squeeze(ddf(1,:,:,:)),medwin(1),medwin(2),medwin(3)); 
ddf(2,:,:,:)=mex_medfilt3_nonzero_float2(squeeze(ddf(2,:,:,:)),medwin(1),medwin(2),medwin(3)); 
figure(4), plot(squeeze(mean(mean(ddf(1,round(nX/2)-10:round(nX/2)+10,round(nY/2)-10:round(nY/2)+10,:),2),3))), hold on
figure(4), plot(squeeze(mean(mean(ddf(2,round(nX/2)-10:round(nX/2)+10,round(nY/2)-10:round(nY/2)+10,:),2),3)))
figure(4), plot(squeeze(dd(1,round(nX/2),round(nY/2),:)))
figure(4), plot(squeeze(dd(2,round(nX/2),round(nY/2),:))), hold off


save -V7.3 cc2d_idc_test_mex2_medfilt_rf_workspace2 
    

figure(4), plot(squeeze(mean(mean(ddf(1,round(nX/2)-10:round(nX/2)+10,round(nY/2)-10:round(nY/2)+10,:),2),3))), hold on
figure(4), plot(squeeze(mean(mean(ddf(2,round(nX/2)-10:round(nX/2)+10,round(nY/2)-10:round(nY/2)+10,:),2),3)))
grid on
xlabel('Time (pixels)'), ylabel('Interframe displacement (pixels)')
legend('Vertical displacement','Horizontal displacement') 
saveFig(gcf,'figures/displacement',400)


for k=1:size(dd,4)
    k
    figure(1), imagesc(squeeze(ddf(1,:,:,k))), drawnow
    %figure(1), imagesc(dbzero(abs(hilbert(double(squeeze(rf(:,:,k)))))),[-45 0]); drawnow
end

figure(1)
for k=1:10:size(dd,4)
    k
    figure(1), imagesc(squeeze(ddf(1,:,:,k)),[-1 1]*15)
    axis equal, axis tight
    drawnow
    fname = ['print -r100 -dpng figures/dd1_movie/dd1_' sprintf('%0.4d',k)]
    eval(fname)
end

figure(1)
for k=1:10:size(dd,4)
    k
    figure(1), imagesc(squeeze(ddf(2,:,:,k)),[-1 1]*15)
    axis equal, axis tight
    drawnow
    fname = ['print -r100 -dpng figures/dd2_movie/dd2_' sprintf('%0.4d',k)]
    eval(fname)
end
