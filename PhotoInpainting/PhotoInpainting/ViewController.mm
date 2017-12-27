//
//  ViewController.m
//  PhotoInpainting
//
//  Created by hncoder on 2017/12/26.
//  Copyright © 2017年 com. All rights reserved.
//

#import "ViewController.h"
#import "MaskView.h"
#include "stdio.h"
#include "stdlib.h"
#include "waveDec2.h"
#include "waveRec2.h"

@interface ViewController ()

@property (nonatomic, strong) UIImageView *imageView;
//@property (nonatomic, strong) MaskView *maskView;

@end

@implementation ViewController

- (UIImageView *)imageView {
    if (!_imageView) {
        _imageView = [[UIImageView alloc] initWithFrame:self.view.bounds];
        _imageView.contentMode = UIViewContentModeScaleAspectFit;
    }
    return _imageView;
}

- (int)waveTest {
    int i, j, k;
    IntDataPanel2D ImageDataSrc;
    WTInfo * pImageDataDec = NULL;
    IntDataPanel2D * pImageDataRec = NULL;
    double * pTracePointer = NULL;
    int nImageRow, nImageCol;
    FILE * fpImageDataFile = NULL;
    
    int nCLen;
    int nTest;
    double dbTest;
    
    /* ∂¡»° ˝æ› */
    if ((fpImageDataFile = fopen("LENA.txt", "r")) == NULL) {
        return 1;
    }
    
    fscanf(fpImageDataFile, "%d", &nImageRow);
    fscanf(fpImageDataFile, "%d", &nImageCol);
    
    ImageDataSrc.m_nSizeRow = nImageRow;
    ImageDataSrc.m_nSizeCol = nImageCol;
    ImageDataSrc.m_pData2D = (int *)calloc(nImageRow*nImageCol, sizeof(int));
    if (ImageDataSrc.m_pData2D == NULL)
        return 1;
    
    for (i = 0; i < nImageRow*nImageCol; i++)
        fscanf(fpImageDataFile, "%d", &(ImageDataSrc.m_pData2D[i]));
    fclose(fpImageDataFile);
    
    /* –°≤®±‰ªª */
    pImageDataDec = waveDec2Int(&ImageDataSrc, 3, WAVE_DB3, ExtT_SYM);
    
    pTracePointer = pImageDataDec->m_pC;
    printf("\n\n------------%d--A-------\n", 0);
    for (i = 0; i < DATA2D(pImageDataDec->m_pS, 0, 0, 2); i++) {
        for (j = 0; j < DATA2D(pImageDataDec->m_pS, 0, 1, 2); j++) {
            printf("%10f\t", *pTracePointer);
            pTracePointer++;
        }
        printf("\n");
    }
    for (k = 1; k < pImageDataDec->m_nWTLevel; k++) {
        printf("\n\n------------%d--H-------\n", k-1);
        for (i = 0; i < DATA2D(pImageDataDec->m_pS, k, 0, 2); i++) {
            for (j = 0; j < DATA2D(pImageDataDec->m_pS, k, 1, 2); j++) {
                printf("%10f\t", *pTracePointer);
                pTracePointer++;
            }
            printf("\n");
        }
        printf("\n\n------------%d--V-------\n", k-1);
        for (i = 0; i < DATA2D(pImageDataDec->m_pS, k, 0, 2); i++) {
            for (j = 0; j < DATA2D(pImageDataDec->m_pS, k, 1, 2); j++) {
                printf("%10f\t", *pTracePointer);
                pTracePointer++;
            }
            printf("\n");
        }
        printf("\n\n------------%d--D-------\n", k-1);
        for (i = 0; i < DATA2D(pImageDataDec->m_pS, k, 0, 2); i++) {
            for (j = 0; j < DATA2D(pImageDataDec->m_pS, k, 1, 2); j++) {
                printf("%10f\t", *pTracePointer);
                pTracePointer++;
            }
            printf("\n");
        }
    }
    //nCLen = DATA2D(pImageDataDec->m_pS, 0, 0, 2)*DATA2D(pImageDataDec->m_pS, 0, 1, 2);;
    //for (i = 1; i < pImageDataDec->m_nWTLevel+1; i++) {
    //    j = DATA2D(pImageDataDec->m_pS, i, 0, 2);
    //    j = DATA2D(pImageDataDec->m_pS, i, 1, 2);
    //    nCLen = DATA2D(pImageDataDec->m_pS, i, 0, 2)*DATA2D(pImageDataDec->m_pS, i, 1, 2)*3;
    //}
    //for (i = 0; i < nCLen; i++) {
    //    // dbTest = pImageDataDec->m_pC[i];
    //    printf("%10f", pImageDataDec->m_pC[i]);
    //}
    
    /* –°≤®ƒÊ±‰ªª */
    printf("\n--------Inverse Result--------\n");
    pImageDataRec = waveRec2Int(pImageDataDec);
    for (i = 0; i < pImageDataRec->m_nSizeRow; i++) {
        for (j = 0; j < pImageDataRec->m_nSizeCol; j++)
            // nTest = DATA2D(pImageDataRec->m_pData2D, i, j, pImageDataRec->m_nSizeCol);
            printf("%5d", DATA2D(pImageDataRec->m_pData2D, i, j, pImageDataRec->m_nSizeCol));
        printf("\n");
    }
    
    WI_FREE(pImageDataDec);
    IDP_FREE(pImageDataRec);
    
    return 0;
}

//- (MaskView *)maskView {
//    if (!_maskView) {
//        _maskView = [[MaskView alloc] init];
//        _maskView.backgroundColor = [UIColor clearColor];
//    }
//
//    return _maskView;
//}

- (void)drawImage
{
    CGImageRef cgimage = [self.imageView.image CGImage];
    
    size_t width = CGImageGetWidth(cgimage); // 图片宽度
    size_t height = CGImageGetHeight(cgimage); // 图片高度
    unsigned char *data = (unsigned char *)calloc(width * height * 4, sizeof(unsigned char)); // 取图片首地址
    size_t bitsPerComponent = 8; // r g b a 每个component bits数目
    size_t bytesPerRow = width * 4; // 一张图片每行字节数目 (每个像素点包含r g b a 四个字节)
    CGColorSpaceRef space = CGColorSpaceCreateDeviceRGB(); // 创建rgb颜色空间
    
    CGContextRef context =
    CGBitmapContextCreate(data,
                          width,
                          height,
                          bitsPerComponent,
                          bytesPerRow,
                          space,
                          kCGImageAlphaPremultipliedLast | kCGBitmapByteOrder32Big);
    CGContextDrawImage(context, CGRectMake(0, 0, width, height), cgimage);
    
    for (size_t i = 0; i < height; i++)
    {
        for (size_t j = 0; j < width; j++)
        {
            size_t pixelIndex = i * width * 4 + j * 4;
            
            unsigned char red = data[pixelIndex];
            unsigned char green = data[pixelIndex + 1];
            unsigned char blue = data[pixelIndex + 2];
            
            // 修改颜色
            red = 0;
            data[pixelIndex] = red;
            data[pixelIndex+1] = green;
            data[pixelIndex+2] = blue;
        }
    }
    
    cgimage = CGBitmapContextCreateImage(context);
    self.imageView.image = [UIImage imageWithCGImage:cgimage];
    
    CGContextRelease(context);
    CGImageRelease(cgimage);
    free(data);
}

- (void)viewDidLoad {
    [super viewDidLoad];
    // Do any additional setup after loading the view, typically from a nib.
    
    [self waveTest];
    
    UIImage *image = [UIImage imageNamed:@"inpainting_2.jpg"];
    
    CGFloat width = [UIScreen mainScreen].bounds.size.width;
    CGFloat height = width * image.size.height / image.size.width;
    if (height > [UIScreen mainScreen].bounds.size.height) {
        height = [UIScreen mainScreen].bounds.size.height;
        width = height * image.size.width / image.size.height;
    }
    
    self.imageView.image = image;
    self.imageView.frame = CGRectMake(0, 0, width, height);
    
    [self.view addSubview:self.imageView];
    self.imageView.center = CGPointMake(self.view.frame.size.width/2.0, self.view.frame.size.height/2.0);
    
    [self drawImage];
    
//    self.maskView.frame = self.imageView.frame;
//    [self.view addSubview:self.maskView];
}


- (void)didReceiveMemoryWarning {
    [super didReceiveMemoryWarning];
    // Dispose of any resources that can be recreated.
}


@end
