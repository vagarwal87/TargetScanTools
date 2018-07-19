// ComputeContextScore.cpp : Defines the entry point for the console application.
// This script was written by Kyle Kai-How Farh
// It computes the context score and 3' pairing score as described in Grimson et al, 2007

//#include "stdafx.h"
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>




#define MAX_3UTR_SIZE	1500



int Min( int nData1, int nData2 ) { return nData1 > nData2 ? nData2 : nData1; }
int Max( int nData1, int nData2 ) { return nData1 < nData2 ? nData2 : nData1; }


double fMin( double fData1, double fData2 ) { return fData1 > fData2 ? fData2 : fData1; }
double fMax( double fData1, double fData2 ) { return fData1 < fData2 ? fData2 : fData1; }

void compute_context_score(	int nSeedTypeID, int nSeedStartPos, int nSeedEndPos, char *s3UTRSeq, int n3UTRSize, char *sCompMiRNASeq, 
							double *fRetContextScore, double *fRetContextScore1, double *fRetContextScore2, double *fRetContextScore3,
							double *fRetParamValue1, double *fRetParamValue2, double *fRetParamValue3 );

void end_program();
void complement_seq( char *sSeq, int nSeqLen );
char complement( char cBase )
{
  return cBase=='A' ? 'T' : cBase=='C' ? 'G' : cBase=='G' ? 'C' : cBase=='T' ? 'A' : cBase;
}



int main(int argc, char* argv[])
{
	int i, nSiteTypeID, nSiteStartPos, nSiteEndPos, n3UTRSize, nMiRNASize;
	
	char s3UTRSeq[20000], smiRCompMatureSeq[10000];
	
	double fScore, fScore1, fScore2, fScore3, fParam1, fParam2, fParam3;


	//argc = 5;
	//char sTemp[5][10000];
	//strcpy( sTemp[1], "8mer" );
	//strcpy( sTemp[2], "1" );
	//strcpy( sTemp[3], "TAAGGGCCAAAAAAAAAAAAAAAAAAA" );
	//strcpy( sTemp[4], "GGCCCTTAGGAAAAGGGCCCTTTTTTTTTTTTTTTTTC" );


	if ( argc != 5 ) { end_program(); }

	// nSeedTypeID =	1 : 8mer
	//					2 : 7m8
	//					3 : 7A1
	//					4 : 6mer
	//
	nSiteTypeID =	strcmp( argv[1], "8mer" ) == 0 ?	1 :
					strcmp( argv[1], "7m8" ) == 0 ?		2 :
					strcmp( argv[1], "7A1" ) == 0 ?		3 :
					strcmp( argv[1], "6mer" ) == 0 ?	4 : -1;

	nMiRNASize =	strlen( argv[3] );
	n3UTRSize =		strlen( argv[4] );
	if ( nSiteTypeID < 1 || nSiteTypeID > 4 || nMiRNASize < 15 || n3UTRSize < nMiRNASize ) { end_program(); }
	
	for(i=0;i<nMiRNASize;i++)
	{
		if ( argv[3][i]!='A' && argv[3][i]!='C' && argv[3][i]!='G' && argv[3][i]!='T' ) { printf("\nERROR. SEQUENCE OTHER THAN A,C,G,T FOUND: %s\n", argv[3] ); exit(0); }
	}
	for(i=0;i<n3UTRSize;i++)
	{
		if ( argv[4][i]!='A' && argv[4][i]!='C' && argv[4][i]!='G' && argv[4][i]!='T' ) { printf("\nERROR. SEQUENCE OTHER THAN A,C,G,T FOUND: %s\n", argv[4] ); exit(0); }
	}

	nSiteStartPos = atoi( argv[2] );
	nSiteEndPos =	nSiteStartPos + ( nSiteTypeID==1 ? 7 : nSiteTypeID==2 ? 6 : nSiteTypeID==3 ? 6 : nSiteTypeID==4 ? 5 : 0 );

	if ( nSiteStartPos < 0 || nSiteStartPos >= n3UTRSize || nSiteEndPos < 0 || nSiteEndPos >= n3UTRSize ) 
	{ 
		printf("\nERROR. nSiteStartPos=%d  nSiteEndPos=%d  n3UTRSize=%d\n", nSiteStartPos, nSiteEndPos, n3UTRSize ); exit(0); 
	}

	strcpy( s3UTRSeq, argv[4] );
	strcpy( smiRCompMatureSeq, argv[3] );
	complement_seq( smiRCompMatureSeq, nMiRNASize );

	// printf("\nsmiRCompMatureSeq = %s", smiRCompMatureSeq );

	compute_context_score(	nSiteTypeID, nSiteStartPos, nSiteEndPos, s3UTRSeq, n3UTRSize, smiRCompMatureSeq,
							&fScore, &fScore1, &fScore2, &fScore3, &fParam1, &fParam2, &fParam3 );

	printf("\n");
	printf("\nTotalContext: %f", fScore );
	printf("\nLocalAU:      %f  %f", fScore1, fParam1 );
	printf("\n3'Pairing:    %f  %f", fScore2, fParam2 );
	printf("\nSitePosition: %f  %f", fScore3, fParam3 );
	printf("\n");
} // end of main()



void end_program()
{
	printf("\n\nERROR. USAGE = ComputeContextScore  SiteType[8mer,7m8,7A1,6mer]  StartPos[0-BASED POS]  MatureMicroRNASeq  3'UTRSeq\n");
	exit(0);
}




// proofread on 04/25/07
void complement_seq( char *sSeq, int nSeqLen )
{
	int i;
	char *sBuffer;

	if ( nSeqLen <= 0 )
	{
		printf("\nERROR. nSeqLen=%d\n", nSeqLen );
		exit(0);
	}

	sBuffer = (char *)calloc( nSeqLen+1, sizeof(char) );
	if ( sBuffer==NULL )
	{
		printf("\n\nMEMORY ALLOCATION ERROR. sBuffer\n\n");
		exit(0);
	}

	for(i=0;i<nSeqLen;i++)
	{
		sBuffer[i] = complement( sSeq[nSeqLen-1-i] );
	}

	for(i=0;i<nSeqLen;i++)
	{
		sSeq[i] = sBuffer[i];
	}

	free( sBuffer );

} // end of void complement_seq()





//
// sCompMiRNASeq : REVERSE-COMPLEMENTARY TO miRNA SEQ
//
void compute_context_score(	int nSeedTypeID, int nSeedStartPos, int nSeedEndPos, char *s3UTRSeq, int n3UTRSize, char *sCompMiRNASeq, 
							double *fRetContextScore, double *fRetContextScore1, double *fRetContextScore2, double *fRetContextScore3,
							double *fRetParamValue1, double *fRetParamValue2, double *fRetParamValue3 )
{
	int i, j, k, nMotifSize, nSeedSize, nDist, nMatureMiRNASize;
	int nStartPos, nEndPos, nSearch3EndPos, nSearch5EndPos, nDist1, nDist2, nOffset, bCoreMatchFound;
	int nMaxUTRPos, nMaxMiRNAPos, nMaxMotifSize, nMaxOffset, nAlreadyPairedSize, nPenaltyNucCount;
	int nUpWindowSize, nDnWindowSize;

	double fScore, fMeanRepression, fContextScoreSum, fParamValue1, fParamValue2, fParamValue3;
	double fSlope1, fSlope2, fSlope3, fIntercept1, fIntercept2, fIntercept3, fContextScore1, fContextScore2, fContextScore3;
	double fMaxScore, fDenominator, fDenominator2;
	double fUpNumerator, fUpDenominator, fDnNumerator, fDnDenominator;

	
	// nSeedTypeID =	1 : 8mer
	//					2 : 7m8
	//					3 : 7A1
	//					4 : 6mer
	if ( nSeedTypeID < 1 || nSeedTypeID > 4 )						{ printf("\nERROR. nSeedTypeID=%d\n", nSeedTypeID ); exit(0); }
	if ( (int)strlen( s3UTRSeq ) != n3UTRSize || n3UTRSize <= 0 )	{ printf("\nERROR. (int)strlen( s3UTRSeq ) = %d  n3UTRSize = %d\n", (int)strlen( s3UTRSeq ), n3UTRSize ); exit(0); }


	char sSiteSeq[1000];
	int nSiteMatchLen = nSeedTypeID==1 ? 7 : nSeedTypeID==2 ? 7 : nSeedTypeID==3 ? 6 : nSeedTypeID==4 ? 6 : 0;
	int nMiRNASeqSize = strlen( sCompMiRNASeq );
	
	if ( nSeedTypeID==1 )		{ strncpy( sSiteSeq, sCompMiRNASeq+nMiRNASeqSize-8, nSiteMatchLen ); }
	else if ( nSeedTypeID==2 )	{ strncpy( sSiteSeq, sCompMiRNASeq+nMiRNASeqSize-8, nSiteMatchLen ); }
	else if ( nSeedTypeID==3 )	{ strncpy( sSiteSeq, sCompMiRNASeq+nMiRNASeqSize-7, nSiteMatchLen ); }
	else if ( nSeedTypeID==4 )	{ strncpy( sSiteSeq, sCompMiRNASeq+nMiRNASeqSize-7, nSiteMatchLen ); }
	else { printf("\nERROR. nSeedTypeID=%d", nSeedTypeID ); exit(0); }

	sSiteSeq[nSiteMatchLen] = '\0'; 
	if ( strncmp( sSiteSeq, s3UTRSeq+nSeedStartPos, nSiteMatchLen ) != 0 )
	{
		printf("\nERROR. SEED SITE NOT FOUND.\n"); exit(0);
	}
	
	
	(*fRetContextScore) =	0.0;
	(*fRetContextScore1) =	0.0;
	(*fRetContextScore2) =	0.0;
	(*fRetContextScore3) =	0.0;

	(*fRetParamValue1) =	0.0;
	(*fRetParamValue2) =	0.0;
	(*fRetParamValue3) =	0.0;



	//***********************************************
	// AU CONTENT
	nSeedSize = nSeedEndPos - nSeedStartPos + 1;

	if (	( nSeedTypeID==1 && nSeedSize != 8 ) || ( nSeedTypeID==2 && nSeedSize != 7 ) ||
			( nSeedTypeID==3 && nSeedSize != 7 ) || ( nSeedTypeID==4 && nSeedSize != 6 ) || 
			nSeedStartPos < 0 || nSeedStartPos >= n3UTRSize || nSeedEndPos < 0 || nSeedEndPos >= n3UTRSize || nSeedStartPos >= nSeedEndPos )
	{
		printf("\nERROR. nSeedTypeID=%d  nSeedStartPos=%d  nSeedEndPos=%d  n3UTRSize=%d\n", nSeedTypeID, nSeedStartPos, nSeedEndPos, n3UTRSize ); exit(0);
	}
	nStartPos =	nSeedTypeID==1 ? (nSeedStartPos- 1) : nSeedTypeID==2 ? (nSeedStartPos- 1) : nSeedTypeID==3 ? (nSeedStartPos- 2) : nSeedTypeID==4 ? (nSeedStartPos- 2) : -1;
	nEndPos =	nSeedTypeID==1 ? (nSeedStartPos-30) : nSeedTypeID==2 ? (nSeedStartPos-30) : nSeedTypeID==3 ? (nSeedStartPos-31) : nSeedTypeID==4 ? (nSeedStartPos-31) : -1;
	nStartPos =	Max( nStartPos, 0 );
	nEndPos =	Max( nEndPos, 0 );

	fUpNumerator =		0.0;
	fUpDenominator =	0.0;
	fDnNumerator =		0.0;
	fDnDenominator =	0.0;
	nUpWindowSize =		0;
	nDnWindowSize =		0;
	int nWindowSize =	nStartPos - nEndPos + 1;

	fParamValue1 =	0.0; 
	fDenominator =	0.0;
	fDenominator2 =	0.0; 
	
	for(i=nStartPos,j=nWindowSize;i>=nEndPos;i--,j--)
	{
		fScore =	(	( nSeedTypeID==1 || nSeedTypeID==2 ) ? (1.0/(double)(nStartPos-i+1)) : 
						( nSeedTypeID==3 || nSeedTypeID==4 ) ? (1.0/(double)(nStartPos-i+2)) : 0.0 );
		
		if ( fScore < 0.0 || fScore > 1.0 ) { printf("\nERROR. fScore=%f  nStartPos=%d  nEndPos=%d  nWindowSize=%d  nSeedStartPos=%d\n", fScore, nStartPos, nEndPos, nWindowSize, nSeedStartPos ); exit(0); }
		
		if ( s3UTRSeq[i]=='A' || s3UTRSeq[i]=='T' || s3UTRSeq[i]=='W' || s3UTRSeq[i]=='a' || s3UTRSeq[i]=='t' || s3UTRSeq[i]=='w' ) { fParamValue1 += fScore; fUpNumerator += fScore; }
		fDenominator +=		fScore;
		fUpDenominator +=	fScore;
		nUpWindowSize++;
	} // end of for i
	

	nStartPos =	Min( nSeedEndPos+ 1, n3UTRSize-1 ); 
	nEndPos =	Min( nSeedEndPos+30, n3UTRSize-1 );
	
	for(i=nStartPos;i<=nEndPos;i++)
	{
		fScore =	(	( nSeedTypeID==1 || nSeedTypeID==3 ) ? (1.0/(double)(i-nStartPos+2)) : 
						( nSeedTypeID==2 || nSeedTypeID==4 ) ? ( i==nStartPos ? 0.5 : (1.0/(double)(i-nStartPos+1)) ) : 0.0 );
		if ( fScore < 0.0 || fScore > 1.0 ) { printf("\nERROR. fScore=%f\n", fScore ); exit(0); }

		if ( s3UTRSeq[i]=='A' || s3UTRSeq[i]=='T' || s3UTRSeq[i]=='W' || s3UTRSeq[i]=='a' || s3UTRSeq[i]=='t' || s3UTRSeq[i]=='w' ) { fParamValue1 += fScore; fDnNumerator += fScore; }
		fDenominator +=		fScore;
		fDnDenominator +=	fScore;
		nDnWindowSize++;
	} // end of for i

	if ( fDenominator > 0.0 ) { fParamValue1 /= fDenominator; }
	if ( fParamValue1 < 0.0 || fParamValue1 > 1.0 ) { printf("\nERROR. fParamValue1 = %f\n", fParamValue1 ); }

	fMeanRepression =	nSeedTypeID==1 ? -0.31 : nSeedTypeID==2 ? -0.161 : nSeedTypeID==3 ? -0.099 : nSeedTypeID==4 ? -0.015 : 0.0;
	
	fSlope1 =			nSeedTypeID==1 ? -0.64 : nSeedTypeID==2 ? -0.50 : nSeedTypeID==3 ? -0.42 : nSeedTypeID==4 ? -0.241 : 0.0;
	fIntercept1 =		nSeedTypeID==1 ? 0.055 : nSeedTypeID==2 ? 0.108 : nSeedTypeID==3 ? 0.137 : nSeedTypeID==4 ? 0.115  : 0.0;
	fContextScore1 =	fSlope1 * fParamValue1 + fIntercept1 - fMeanRepression;


	//***********************************************
	// 3' PAIRING
	nMatureMiRNASize = (int)strlen( sCompMiRNASeq );
	if ( nMatureMiRNASize < 17 ) { printf("\nERROR. nMatureMiRNASize=%d\n", nMatureMiRNASize ); exit(0); }

	nSearch3EndPos =	Max(0, nSeedStartPos-1 );
	nSearch5EndPos =	Max(0, nSearch3EndPos-30);
	if ( nSearch3EndPos < -2 || nSearch3EndPos >= n3UTRSize || nSearch5EndPos > nSearch3EndPos )
	{
		printf("\nERROR. nSearch3EndPos=%d  nSearch5EndPos=%d  n3UTRSize=%d\n", nSearch3EndPos, nSearch5EndPos, n3UTRSize );
		exit(0);
	}
	nAlreadyPairedSize = ( nSeedTypeID==1 || nSeedTypeID==2 ) ? 8 : 7;

	fMaxScore = 0.0; nMaxUTRPos = -1; nMaxMiRNAPos = -1; nMaxMotifSize = -1; nMaxOffset = 0;
	if ( nSearch3EndPos == -1 || nSearch3EndPos == -2 )
	{
		fMaxScore = 0.0;
	}
	else
	{
		for(i=nSearch5EndPos;i<=nSearch3EndPos;i++)
		{
			for(j=0;j<=(nMatureMiRNASize-nAlreadyPairedSize-1);j++)
			{
				for(nMotifSize=(nMatureMiRNASize-nAlreadyPairedSize);nMotifSize>=1;nMotifSize--)
				{
					if ( (i+nMotifSize-1) <= nSearch3EndPos && (j+nMotifSize-1) <= (nMatureMiRNASize-nAlreadyPairedSize-1) )
					{
						if ( strncmp( s3UTRSeq+i, sCompMiRNASeq+j, nMotifSize ) == 0 ) 
						{
							fScore = 0.0; bCoreMatchFound = 0;
							for(k=j;k<=(j+nMotifSize-1);k++)
							{
								if ( k >= (nMatureMiRNASize-16) && k <= (nMatureMiRNASize-13) )			{ fScore += 1.0; bCoreMatchFound = 1; }
								else																	{ fScore += 0.5; }
							}
							nDist1 = nSearch3EndPos+1 -							(i+nMotifSize-1);
							nDist2 = (nMatureMiRNASize-nAlreadyPairedSize) -	(j+nMotifSize-1);
							if ( nDist1 < 0 || nDist1 > (nSearch3EndPos+1) || nDist2 < 0 || nDist2 > (nMatureMiRNASize-nAlreadyPairedSize) )
							{
								printf(	"\nERROR. nMotifSize=%d  nDist1=%d  nSearch3EndPos=%d  nDist2=%d  nMatureMiRNASize-nAlreadyPairedSize)=%d\n", 
										nMotifSize, nDist1, nSearch3EndPos, nDist2, nMatureMiRNASize-nAlreadyPairedSize );
								exit(0);
							}
							
							nOffset = abs( nDist1 - nDist2 );
							if ( nOffset <= 2 ) { nPenaltyNucCount = 0;				}
							else				{ nPenaltyNucCount = nOffset - 2;	}
							
							fScore -= ( (double)nPenaltyNucCount * 0.5 ); // OFFSET PENALTY

							if ( fScore > fMaxScore )
							{
								fMaxScore =		fScore;
								nMaxUTRPos =	i;
								nMaxMiRNAPos =	j;
								nMaxMotifSize =	nMotifSize;
								nMaxOffset =	nOffset;
							} // end of if fMaxScore

						} // end of if strncmp()

					} // end of if (i+nMotifSize-1) ...
				} // end of for k
			} // end of for j
		} // end of for i
	} // end of if nSearch3EndPos else

	//exit(0);

	if ( fMaxScore < 0.0 )
	{
		printf("\nERROR. fMaxScore=%f  nSearch3EndPos=%d\n\ns3UTRSeq = %s\n\nsCompMiRNASeq = %s\n", fMaxScore, nSearch3EndPos, s3UTRSeq, sCompMiRNASeq );
		exit(0);
	}

	fParamValue2 =		fMaxScore;
	fSlope2 =			nSeedTypeID==1 ? -0.0041 : nSeedTypeID==2 ? -0.031 : nSeedTypeID==3 ? -0.0211 : nSeedTypeID==4 ? -0.00278 : 0.0;
	fIntercept2 =		nSeedTypeID==1 ? -0.299  : nSeedTypeID==2 ? -0.094 : nSeedTypeID==3 ? -0.053  : nSeedTypeID==4 ? -0.0091  : 0.0;
	fContextScore2 =	fSlope2 * fParamValue2 + fIntercept2 - fMeanRepression;


	//***********************************************
	// SEED POSITION
	nDist = Min( Min( abs( 1 - ((nSeedStartPos+1+nSeedEndPos+1)/2) ), abs( n3UTRSize - ((nSeedStartPos+1+nSeedEndPos+1)/2) ) ), MAX_3UTR_SIZE ); // ALLOW DIST UP TO ONLY 1500bp
	if ( nDist <= 0 || nDist >= n3UTRSize )
	{
		printf("\nERROR. nSeedStartPos=%d  nSeedEndPos=%d  nDist=%d  n3UTRSize=%d\n", nSeedStartPos, nSeedEndPos, nDist, n3UTRSize );
		exit(0);
	}
	fParamValue3 =		(double)nDist;
	fSlope3 =			nSeedTypeID==1 ? 0.000172 : nSeedTypeID==2 ? 0.000091 : nSeedTypeID==3 ? 0.000072 : nSeedTypeID==4 ? 0.000049 : 0.0;
	fIntercept3 =		nSeedTypeID==1 ? -0.38    : nSeedTypeID==2 ? -0.198   : nSeedTypeID==3 ? -0.131   : nSeedTypeID==4 ? -0.033   : 0.0;
	fContextScore3 =	fSlope3 * fParamValue3 + fIntercept3 - fMeanRepression;

	fContextScoreSum =		fMeanRepression + fContextScore1 + fContextScore2 + fContextScore3;
	fContextScoreSum =		fMin( fContextScoreSum, 0.0 ); // ACTIVATION IS NOT PERMITTED.

	(*fRetContextScore) =	fContextScoreSum;
	(*fRetContextScore1) =	fContextScore1;
	(*fRetContextScore2) =	fContextScore2;
	(*fRetContextScore3) =	fContextScore3;

	(*fRetParamValue1) =	fParamValue1;
	(*fRetParamValue2) =	fParamValue2;
	(*fRetParamValue3) =	fParamValue3;

	
	if ( nSeedStartPos <= 14 ) // RIBOSOMAL SHADOWING
	{
		printf("\nWARNING: RIBOSOMAL SHADOWING. SEED START POS <= 15. CONTEXT SCORE NOT COMPUTED.");
		(*fRetContextScore) =	0.0;
		(*fRetContextScore1) =	0.0;
		(*fRetContextScore2) =	0.0;
		(*fRetContextScore3) =	0.0;
	}
	
} // end of void compute_context_score()


