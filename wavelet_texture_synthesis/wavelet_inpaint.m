
function [inpaintedImg,origImg,fillImg,C,D,time,iter] = wavelet_inpaint(imgFilename,fillFilename,fillColor)

warning off MATLAB:divideByZero

[img,fillImg,fillRegion] = loadimgs(imgFilename,fillFilename,fillColor);
img = double(img);
s = size(fillRegion);
fillRegion = double(fillRegion);
origImg = img;

global fRegion;
[imgA,imgH,imgV,imgD,fRegion] = preprocessing(img,fillRegion);

origImgA = imgA;
origImgH = imgH;
origImgV = imgV;
origImgD = imgD;
ind = img2ind(imgA);
sz = [size(imgA,1) size(imgA,2)];


global sourceRegion;
sourceRegion = ~fRegion;

% Initialize isophote values
global Ix;
global Iy;
Ix = sum(imgH,3)/(3*255); Iy = sum(imgV,3)/(3*255);
temp = Ix; Ix = -Iy; Iy = temp;   % Rotate gradient 90 degrees

% Initialize confidence and data terms
C = double(sourceRegion);
D = repmat(-.1,sz);

% Visualization stuff


%Initialize details record
global detail;
detail = (sum(imgH,3)/(3*255)).^2+(sum(imgV,3)/(3*255)).^2+(sum(imgD,3)/(3*255)).^2;
detail = normr(detail);
detail(~isfinite(detail))=0;

iter = 0;
tb = clock;
% Seed 'rand' for reproducible results (good for testing)
rand('state',0);
% Loop until entire fill region has been covered
while any(fRegion(:))
    % Find contour & normalized gradients of fill region
    iter = iter + 1;
    dR = find(conv2(fRegion,[1,1,1;1,-8,1;1,1,1],'same')>0);
    [Nx,Ny] = gradient(double(~fRegion));
   % N = [Nx(dR(:)) Ny(dR(:))];
   % N = normr(N);  
   % N(~isfinite(N))=0; % handle NaN and Inf
      
    for k=dR'
        % Compute confidences along the fill front
        [Hp,q] = getadaptivepatch(sz,k);
        %q = Hp(~(fillRegion(Hp)));
        C(k) = sum(C(q))/numel(Hp);
       
        % Compute confidences along the fill front
       % Hpp = Hp';dRR = dR';
        Hpp = matrix2vector(Hp);
        bP = intersect(Hpp,dR);
        NbP = [Nx(bP(:)) Ny(bP(:))];
        NbP = normr(NbP);  
        NbP(~isfinite(NbP))=0;
        Nk = [sum(NbP(:,1)) sum(NbP(:,2))]/numel(bP);
        
        Ixy = [Ix(q(:)) Iy(q(:))];
        Ikxy = [sum(Ixy(:,1)) sum(Ixy(:,2))]/numel(q);
        
        D(k) = abs(Ikxy(1)*Nk(1)+Ikxy(2)*Nk(2)) + 0.001;
    end
    priorities = C(dR).* D(dR);
    
    % Find patch with maximum priority, Hp
    [unused,ndx] = max(priorities(:));
    p = dR(ndx(1));
    [Hp,q,rows,cols,texture] = getadaptivepatch(sz,p);
    toFill = fRegion(Hp);
    
    toFill = ~(~toFill);
    
     Ixy = [Ix(q(:)) Iy(q(:))];
     direct = [sum(Ixy(:,1)) sum(Ixy(:,2))]/numel(q);
   
    p_topleft = [rows(1) cols(1)];
   
   
   [Hq] = bestexemplar(imgA,imgA(rows,cols,:),toFill,p_topleft,direct,p,texture);
    
    % Find exemplar that minimizes error, Hq
   % Hq = bestexemplar(imgA,imgA(rows,cols,:),toFill',sourceRegion);
    %Hq = bestexemplar(imgA,imgA(rows,cols,:),toFill);
    %Hq = bestexemplar(imgA,imgA(rows,cols,:),toFill,pp,direct);
    
     % Update fill region
    fRegion(Hp(toFill)) = false;
    
    %sourceRegion(Hp(toFill)) = true;
  
    % Propagate confidence & isophote values
    C(Hp(toFill))  = C(p);
    Ix(Hp(toFill)) = Ix(Hq(toFill));
    Iy(Hp(toFill)) = Iy(Hq(toFill));
    detail(Hp(toFill)) = detail(Hq(toFill));
    % Copy image data from Hq to Hp
    ind(Hp(toFill)) = ind(Hq(toFill));
    imgA(rows,cols,:) = ind2img(ind(rows,cols),origImgA);  
    imgH(rows,cols,:) = ind2img(ind(rows,cols),origImgH);  
    imgV(rows,cols,:) = ind2img(ind(rows,cols),origImgV);  
    imgD(rows,cols,:) = ind2img(ind(rows,cols),origImgD);  
    
    
end

inpaintedImg = wavelet2img(imgA,imgH,imgV,imgD);
te = clock;
time = etime(te,tb);
img;

%---------------------------------------------------------------------
% Loads the an image and it's fill region, using 'fillColor' as a marker
% value for knowing which pixels are to be filled.
%---------------------------------------------------------------------
function [img,fillImg,fillRegion] = loadimgs(imgFilename,fillFilename,fillColor)
img = imread(imgFilename); fillImg = imread(fillFilename);
fillRegion = fillImg(:,:,1)==fillColor(1) & ...
    fillImg(:,:,2)==fillColor(2) & fillImg(:,:,3)==fillColor(3);

%---------------------------------------------------------------------
% Returns the wavelet datas by preprocessing
%---------------------------------------------------------------------
function [imgA,imgH,imgV,imgD,fRegion] = preprocessing(img,fillRegion)
for i=1:3
    [c,s] = wavedec2(img(:,:,i),1,'haar');
    imgA(:,:,i) = appcoef2(c,s,'haar',1);
    imgH(:,:,i) =detcoef2('h',c,s,1);%水平方向
    imgV(:,:,i) =detcoef2('v',c,s,1);%垂直方向
    imgD(:,:,i) =detcoef2('d',c,s,1);%斜线方向
end
[cm,sm] = wavedec2(fillRegion,1,'haar');
global fRegion;
fRegion = appcoef2(cm,sm,'haar',1);
fRegion(:) = fRegion(:)~=0


%---------------------------------------------------------------------
% Returns the indices for a 9x9 patch centered at pixel p.
%---------------------------------------------------------------------
function [Hp,q,rows,cols,texture] = getadaptivepatch(sz,p)
% [x,y] = ind2sub(sz,p);  % 2*w+1 == the patch size
w=4; p=p-1; y=floor(p/sz(1))+1; p=rem(p,sz(1)); x=floor(p)+1;
global detail;
global fRegion;
rows = max(x-w,1):min(x+w,sz(1));
cols = (max(y-w,1):min(y+w,sz(2)))';
Hp = sub2ndx(rows,cols,sz(1));
q = Hp(~(fRegion(Hp)));

localdetail = detail(Hp);
localfill = fRegion(Hp);
outdamagedetail = localdetail(~localfill);
hasdetail = find(outdamagedetail>0);
%infRegion = find(fRegion(hasdetail)>0);
texture = (numel(hasdetail)/numel(outdamagedetail));
while texture<0.5 & w < 7
    w = w+1;
    rows = max(x-w,1):min(x+w,sz(1));
    cols = (max(y-w,1):min(y+w,sz(2)))';
    Hp = sub2ndx(rows,cols,sz(1));
    q = Hp(~(fRegion(Hp)));
 %   hasdetail = find(detail(Hp)>0);
 %   infRegion = find(fRegion(hasdetail)>0);
 
    localdetail = detail(Hp);
    localfill = fRegion(Hp);
    outdamagedetail = localdetail(~localfill);
    hasdetail = find(outdamagedetail>0);
    texture = (numel(hasdetail)/numel(outdamagedetail));
end


%---------------------------------------------------------------------
% Returns the indices for a 9x9 patch centered at pixel p.
%---------------------------------------------------------------------
function [Hp,rows,cols] = getpatch(sz,p)
% [x,y] = ind2sub(sz,p);  % 2*w+1 == the patch size
w=4; p=p-1; y=floor(p/sz(1))+1; p=rem(p,sz(1)); x=floor(p)+1;
rows = max(x-w,1):min(x+w,sz(1));
cols = (max(y-w,1):min(y+w,sz(2)))';
Hp = sub2ndx(rows,cols,sz(1));

%---------------------------------------------------------------------
% Converts the (rows,cols) subscript-style indices to Matlab index-style
% indices.  Unforunately, 'sub2ind' cannot be used for this.
%---------------------------------------------------------------------
function N = sub2ndx(rows,cols,nTotalRows)
X = rows(ones(length(cols),1),:);
Y = cols(:,ones(1,length(rows)));
N = X+(Y-1)*nTotalRows;

function V = matrix2vector(matrix)
[n m] = size(matrix);
t=1;
for k=0:n-1
    V(t:t+m-1) = matrix(k+1,1:m);
    t = t + k*m;
end

%---------------------------------------------------------------------
% Converts an indexed image into an RGB image, using 'img' as a colormap
%---------------------------------------------------------------------
function img2 = ind2img(ind,img)
for i=3:-1:1, temp=img(:,:,i); img2(:,:,i)=temp(ind); end;


%---------------------------------------------------------------------
% Converts an RGB image into a indexed image, using the image itself as
% the colormap.
%---------------------------------------------------------------------
function ind = img2ind(img)
s=size(img); ind=reshape(1:s(1)*s(2),s(1),s(2));

%---------------------------------------------------------------------
% Scans over the entire image (with a sliding window)
% for the exemplar with the lowest error. Calls a MEX function.
%---------------------------------------------------------------------
function Hq = bestexemplar(img,Ip,toFill,p_topleft,direct,p,texture)
m=size(Ip,1); mm=size(img,1); n=size(Ip,2); nn=size(img,2);
imgsz = [mm nn]; Ipsz = [m n];
ind = img2ind(img);
global sourceRegion;
kown = ind(sourceRegion);


%[Sp,srows,scols] = getSearchArea(imgsz,p,Ipsz);

%lsp = matrix2vector(Sp);

%search = intersect(lsp,kown);

mindiff = +inf;

in = false;
%[Hq] = getadaptivepatch(imgsz,p);
for i = kown'
    [isPatchPoint distance rows cols angle] = forSearchBestPatch(i,p_topleft,direct,imgsz,Ipsz,toFill,ind);
    if isPatchPoint
        patch = img(rows,cols,:);
        mex = sum(sqrt((Ip-patch).^2),3)/3;
        kownarea_mex = mex(~toFill);
        %diff = (sum(kownarea_mex)/numel(kownarea_mex))*(1+(distance/50)^((1+texture)))*(1+(angle/(pi/2))^(4*(1/(1+texture))));
        diff = (sum(kownarea_mex)/numel(kownarea_mex))*(1+(distance/25)^(2*(1+texture)))*(1+(angle/(pi/2))^(3*(1/(1+texture))));
        %*(1+distance/200);
        if diff < mindiff
            mindiff = diff;
            Hq = sub2ndx(rows,cols,imgsz(1));
            %rowsq = rows;
            %colsq = cols;
            in = true;
        end
    end
end
if ~in
    [Hq] = getadaptivepatch(imgsz,p);
end

function angle = angleOfTowVector(v1,v2)
if norm(v1)==0 | norm(v2)==0
    angle=0;
else
     C = dot(v1,v2)/(norm(v1)*norm(v2));
    angle = acos(C);
    if(angle > pi/2.0)
        angle = pi - angle;
    end
end

function [isPatchPoint distance rows cols angle] = forSearchBestPatch(i,p_topleft,direct,imgsz,Ipsz,toFill,ind)
p = i;
m = imgsz(1);
i=i-1; y=floor(i/m)+1; i=rem(i,m); x=floor(i)+1; 
iDirect=[x y]-p_topleft;
%isOnDirect = angleOfTowVector(iDirect,direct)<=pi/2;
distance = sqrt(sum(iDirect.^2));
isSearchArea = distance < 25; 
isPatchPoint = isSearchArea;

rows=0;cols=0;angle=0;
global Ix;
global Iy;
if isPatchPoint
    [issample rows cols] = searchPatch(imgsz,Ipsz,p);
    isPatchPoint = issample;
    if isPatchPoint
        %ind = img2ind(img);
        patchind = ind(rows,cols);
        konwarea_pi = patchind(toFill);       
        patchIxy = [Ix(konwarea_pi(:)) Iy(konwarea_pi(:))];
        patch_direct = [sum(patchIxy(:,1)) sum(patchIxy(:,2))]/numel(konwarea_pi);
        angle = angleOfTowVector(patch_direct,direct);
    end
end




function [issample rows cols] = searchPatch(imgsz,Ipsz,p)
global sourceRegion;
p=p-1; y=floor(p/imgsz(1))+1; p=rem(p,imgsz(1)); x=floor(p)+1; 
rows = max(x,1):min(x+Ipsz(1)-1,imgsz(1));
cols = (max(y,1):min(y+Ipsz(2)-1,imgsz(2)))';
Hp = sub2ndx(rows,cols,imgsz(1));
issample = (numel(rows)==Ipsz(1) & numel(cols)==Ipsz(2));
%issample = size(Hp)==Ipsz;
if issample
    issample = sum((sum(sourceRegion(Hp)))')==numel(Hp); 
    % issample = ~any((any(~sourceRegion(rows,cols)))');
    %issample = any(~sourceRegion(rows,cols));
%    sourceRegion(rows,cols');
 %   for i=rows
  %      for j=col
end
    
    

function inpaintedImg = wavelet2img(imgA,imgH,imgV,imgD)
for i=1:3
    inpaintedImg(:,:,i) = idwt2(imgA(:,:,i),imgH(:,:,i),imgV(:,:,i),imgD(:,:,i),'haar');
end


