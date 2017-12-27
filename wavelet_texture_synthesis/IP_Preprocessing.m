function [CA,CH,CV,CD,WTMask] = IP_Preprocessing(unrepairedImage)

[M,N] = size(unrepairedImage);
for i=1:M
    for j=1:N
        if (unrepairedImage(i,j) == 255)
            Mask(i,j) = 1;
        else
            Mask(i,j) = 0;
        end
    end
end

[c,s] = wavedec2(unrepairedImage,1,'haar');
CA = appcoef2(c,s,'haar',1);
CH =detcoef2('h',c,s,1);%水平方向
CV =detcoef2('v',c,s,1);%垂直方向
CD =detcoef2('d',c,s,1);%斜线方向

[cm,sm] = wavedec2(Mask,1,'haar');
WTMask = appcoef2(cm,sm,'haar',1);
[WTM,WTN] = size(WTMask);
for i=1:WTM
    for j=1:WTN
        if(WTMask(i,j) ~= 0)
            WTMask(i,j) = 1;
        end
    end
end


