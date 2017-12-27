#include "utilityFunWT.h"

DBDataPanel2D * SubWaveRec2(WTCOEFSet * pSrcWTCoef, WaveletBASE * pWaveBase);

DBDataPanel2D * waveRec2(WTInfo * pSrcWTInfo)
{
	int i, j;
	DBDataPanel2D * pTempCoefA_LL = NULL;
	WaveletBASE * pWaveletBase = NULL;
	DBDataPanel2D * pRecResult = NULL;

	WTCOEFSet WTCoefSet;
	DBDataPanel2D pCoefPanel[4];

	if (!(pRecResult = (DBDataPanel2D*)calloc(1, sizeof(DBDataPanel2D))))
		exit(1);
	/* 根据nWaveType获取小波基的滤波器系数 */
	if (!(pWaveletBase = SetWaveletBase(pSrcWTInfo->m_nWaveType))) {
		DDP_FREE(pRecResult);
		exit(1);
	}

	/* 小波系数集合与四个2D数据关联 */
	WTCoefSet.A_LL = &(pCoefPanel[0]);	WTCoefSet.H_LH = &(pCoefPanel[1]);
	WTCoefSet.V_HL = &(pCoefPanel[2]);	WTCoefSet.D_HH = &(pCoefPanel[3]);

	/* 获取最高层的小波系数，0-低通分量，1-水平分量，2-垂直分量，3-对角分量 */
	pCoefPanel[0].m_pData2D = pSrcWTInfo->m_pC;
	pCoefPanel[0].m_nSizeRow = DATA2D(pSrcWTInfo->m_pS, 0, 0, 2);
	pCoefPanel[0].m_nSizeCol = DATA2D(pSrcWTInfo->m_pS, 0, 1, 2);
	for (i = 1; i < 4; i++) {
		pCoefPanel[i].m_pData2D = pCoefPanel[i-1].m_pData2D+pCoefPanel[i-1].m_nSizeRow*pCoefPanel[i-1].m_nSizeCol;
		pCoefPanel[i].m_nSizeRow = DATA2D(pSrcWTInfo->m_pS, 1, 0, 2);
		pCoefPanel[i].m_nSizeCol = DATA2D(pSrcWTInfo->m_pS, 1, 1, 2);
	}

	for (i = 1; i<= pSrcWTInfo->m_nWTLevel; i++) {
		pTempCoefA_LL = SubWaveRec2(&WTCoefSet, pWaveletBase);
		WTCoefSet.A_LL = GetSubMatrix(pTempCoefA_LL, DATA2D(pSrcWTInfo->m_pS, i+1, 0, 2), DATA2D(pSrcWTInfo->m_pS, i+1, 1, 2));
		DDP_FREE(pTempCoefA_LL);
		if (i == pSrcWTInfo->m_nWTLevel)
			break;
		for (j = 1; j < 4; j++) {
			if (j == 1)
				pCoefPanel[j].m_pData2D = pCoefPanel[3].m_pData2D+pCoefPanel[3].m_nSizeRow*pCoefPanel[3].m_nSizeCol;
			else
				pCoefPanel[j].m_pData2D = pCoefPanel[j-1].m_pData2D+pCoefPanel[j-1].m_nSizeRow*pCoefPanel[j-1].m_nSizeCol;
			pCoefPanel[j].m_nSizeRow = DATA2D(pSrcWTInfo->m_pS, i+1, 0, 2);
			pCoefPanel[j].m_nSizeCol = DATA2D(pSrcWTInfo->m_pS, i+1, 1, 2);
		}
	}

	pRecResult = WTCoefSet.A_LL;
	WB_FREE(pWaveletBase);
	return pRecResult;
}


DBDataPanel2D * SubWaveRec2(WTCOEFSet * pSrcWTCoef, WaveletBASE * pWaveBase)
{
	DBDataPanel2D * pResultA_LL = NULL;

	DBDataPanel2D * pTempDataA_LL = NULL, * pTempDataH_LH = NULL, * pTempDataV_HL = NULL, * pTempDataD_HH = NULL;
	DBDataPanel2D * pTempDataA_LL2 = NULL, * pTempDataH_LH2 = NULL, * pTempDataV_HL2 = NULL, * pTempDataD_HH2 = NULL;
	DBDataPanel2D * pTempDataAH = NULL, * pTempDataVD = NULL;
	DBDataPanel2D * pTempDataAH2 = NULL, * pTempDataVD2 = NULL;
	
	/* 行插值 */
	pTempDataA_LL = DYADUP(pSrcWTCoef->A_LL, DYAD_ODD, DIR_ROW);
	pTempDataH_LH = DYADUP(pSrcWTCoef->H_LH, DYAD_ODD, DIR_ROW);
	pTempDataV_HL = DYADUP(pSrcWTCoef->V_HL, DYAD_ODD, DIR_ROW);
	pTempDataD_HH = DYADUP(pSrcWTCoef->D_HH, DYAD_ODD, DIR_ROW);

	/* 列滤波 */
	pTempDataA_LL2 = CONV2D(pTempDataA_LL, pWaveBase->LO_R, pWaveBase->nFilterLen, DIR_COL, CONVT_FULL);
	pTempDataH_LH2 = CONV2D(pTempDataH_LH, pWaveBase->HI_R, pWaveBase->nFilterLen, DIR_COL, CONVT_FULL);
	pTempDataV_HL2 = CONV2D(pTempDataV_HL, pWaveBase->LO_R, pWaveBase->nFilterLen, DIR_COL, CONVT_FULL);
	pTempDataD_HH2 = CONV2D(pTempDataD_HH, pWaveBase->HI_R, pWaveBase->nFilterLen, DIR_COL, CONVT_FULL);
	DDP_FREE(pTempDataA_LL);	DDP_FREE(pTempDataH_LH);	DDP_FREE(pTempDataV_HL);	DDP_FREE(pTempDataD_HH);

	/* 组合后列差值 */
	pTempDataAH = SumMatrix(pTempDataA_LL2, pTempDataH_LH2);
	DDP_FREE(pTempDataA_LL2);	DDP_FREE(pTempDataH_LH2);
	pTempDataAH2 = DYADUP(pTempDataAH, DYAD_ODD, DIR_COL);
	DDP_FREE(pTempDataAH);

	pTempDataVD = SumMatrix(pTempDataV_HL2, pTempDataD_HH2);
	DDP_FREE(pTempDataV_HL2);	DDP_FREE(pTempDataD_HH2);
	pTempDataVD2 = DYADUP(pTempDataVD, DYAD_ODD, DIR_COL);
	DDP_FREE(pTempDataVD);

	/* 行滤波 */
	pTempDataAH = CONV2D(pTempDataAH2, pWaveBase->LO_R, pWaveBase->nFilterLen, DIR_ROW, CONVT_FULL);
	DDP_FREE(pTempDataAH2);
	pTempDataVD = CONV2D(pTempDataVD2, pWaveBase->HI_R, pWaveBase->nFilterLen, DIR_ROW, CONVT_FULL);
	DDP_FREE(pTempDataVD2);

	/* 组合 */
	pResultA_LL = SumMatrix(pTempDataAH, pTempDataVD);

	DDP_FREE(pTempDataAH);
	DDP_FREE(pTempDataVD);
	return pResultA_LL;
}


IntDataPanel2D * waveRec2Int(WTInfo * pSrcWTInfo)
{
	int i;
	IntDataPanel2D * pRecResult = NULL;
	DBDataPanel2D * pDBRecResult = NULL;

	if (!pSrcWTInfo)
		exit(1);

	if (!(pRecResult = (IntDataPanel2D *)calloc(1, sizeof(IntDataPanel2D))))
		exit(1);

	if (!(pDBRecResult = waveRec2(pSrcWTInfo))) {
		IDP_FREE(pRecResult);
		exit(1);
	}
	pRecResult->m_nSizeRow = pDBRecResult->m_nSizeRow;
	pRecResult->m_nSizeCol = pDBRecResult->m_nSizeCol;

	if (!(pRecResult->m_pData2D = (int *)calloc(pRecResult->m_nSizeRow*pRecResult->m_nSizeCol, sizeof(int)))) {
		IDP_FREE(pRecResult);
		DDP_FREE(pDBRecResult);
	};

	for (i = 0; i < (pRecResult->m_nSizeRow)*(pRecResult->m_nSizeCol); i++)
		pRecResult->m_pData2D[i] = (unsigned int)pDBRecResult->m_pData2D[i];

	DDP_FREE(pDBRecResult);
	return pRecResult;
}
