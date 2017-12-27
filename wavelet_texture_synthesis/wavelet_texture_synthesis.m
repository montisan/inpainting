
[inpaintedImg,origImg,fillImg,C,D,time,iter] = wavelet_inpaint('forest_s.PNG','forest_area.PNG',[0 255 0])
figure;
subplot(231);image(uint8(origImg)); title('Original image');
subplot(232);image(uint8(fillImg)); title('Fill region');
subplot(233);image(uint8(inpaintedImg)); title('Inpainted image');
subplot(234);imagesc(C); title('Confidence term');
subplot(235);imagesc(D); title('Data term');

figure;
subplot(121);imagesc(C); title('Confidence term');
subplot(122);imagesc(D); title('Data term');
imwrite(uint8(inpaintedImg),'f_inpainted.bmp');


