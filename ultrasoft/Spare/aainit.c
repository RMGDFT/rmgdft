/************************** SVN Revision Information **************************
 **    $Id$    **
******************************************************************************/

#include <stdio.h>
#include <float.h>
#include <stdlib.h>
#include <math.h>
#include "main.h"

void aainit (int lli, int mix, int lx, int mx, int nlx, REAL ap[][9][9], int lpx[][9],
             int lpl[][9][9])
{
    int ll, lp, mp, np, l, m, k, n, p, il, ik, ii, mpp, li, lj, mi, mj, ilp;
    REAL y = 0.0, f, t, s, fs, ss, fts, ts;
    REAL cc[3][5][3][5][5];
    REAL aaR[5][9][3][5][3][5], aaI[5][9][3][5][3][5];
    REAL uR[5][9][9], uI[5][9][9];
    REAL sumR, sumI;


    ll = 2 * lli - 1;
    if (ll > lx)
        error_handler ("lx too small");

    /*initialize cc */
    for (l = 0; l < 3; l++)
    {
        for (m = 0; m < 5; m++)
        {
            for (k = 0; k < 3; k++)
            {
                for (n = 0; n < 5; n++)
                {
                    for (p = 0; p < 5; p++)
                        cc[l][m][k][n][p] = 0.0;
                }
            }
        }
    }


    cc[0][0][0][0][0] = 1.0;

    if (lli > 1)
    {
        cc[0][0][1][0][1] = 1.0;
        cc[0][0][1][1][1] = 1.0;
        cc[0][0][1][2][1] = 1.0;

        cc[1][0][1][0][2] = sqrt (6.0 / 5.0);
        cc[1][0][1][1][2] = sqrt (3.0 / 5.0);
        cc[1][0][1][2][0] = -1.0;
        cc[1][0][1][2][2] = sqrt (1.0 / 5.0);
        cc[1][1][1][1][0] = 1.0;
        cc[1][1][1][1][2] = sqrt (4.0 / 5.0);
        cc[1][1][1][2][2] = sqrt (3.0 / 5.0);
        cc[1][2][1][2][2] = sqrt (6.0 / 5.0);
    }

    if (lli > 2)
    {
        cc[0][0][2][0][2] = 1.0;
        cc[0][0][2][1][2] = 1.0;
        cc[0][0][2][2][2] = 1.0;
        cc[0][0][2][3][2] = 1.0;
        cc[0][0][2][4][2] = 1.0;

        cc[1][0][2][0][3] = sqrt (45.0 / 35.0);
        cc[1][0][2][1][3] = sqrt (30.0 / 35.0);
        cc[1][0][2][2][1] = -sqrt (1.0 / 5.0);
        cc[1][0][2][2][3] = sqrt (18.0 / 35.0);
        cc[1][0][2][3][1] = -sqrt (3.0 / 5.0);
        cc[1][0][2][3][3] = sqrt (9.0 / 35.0);
        cc[1][0][2][4][1] = -sqrt (6.0 / 5.0);
        cc[1][0][2][4][3] = sqrt (3.0 / 35.0);
        cc[1][1][2][0][3] = sqrt (15.0 / 35.0);
        cc[1][1][2][1][1] = sqrt (3.0 / 5.0);
        cc[1][1][2][1][3] = sqrt (24.0 / 35.0);
        cc[1][1][2][2][1] = sqrt (4.0 / 5.0);

        cc[1][1][2][2][3] = sqrt (27.0 / 35.0);
        cc[1][1][2][3][1] = sqrt (3.0 / 5.0);
        cc[1][1][2][3][3] = sqrt (24.0 / 35.0);
        cc[1][1][2][4][3] = sqrt (15.0 / 35.0);
        cc[1][2][2][0][1] = -sqrt (6.0 / 5.0);
        cc[1][2][2][0][3] = sqrt (3.0 / 35.0);
        cc[1][2][2][1][1] = -sqrt (3.0 / 5.0);
        cc[1][2][2][1][3] = sqrt (9.0 / 35.0);
        cc[1][2][2][2][1] = -sqrt (1.0 / 5.0);
        cc[1][2][2][2][3] = sqrt (18.0 / 35.0);
        cc[1][2][2][3][3] = sqrt (30.0 / 35.0);
        cc[1][2][2][4][3] = sqrt (45.0 / 35.0);

        cc[2][0][2][0][4] = sqrt (70.0 / 49.0);
        cc[2][0][2][1][4] = sqrt (35.0 / 49.0);
        cc[2][0][2][2][2] = -sqrt (20.0 / 49.0);
        cc[2][0][2][2][4] = sqrt (15.0 / 49.0);
        cc[2][0][2][3][2] = -sqrt (30.0 / 49.0);
        cc[2][0][2][3][4] = sqrt (5.0 / 49.0);
        cc[2][0][2][4][0] = 1.0;
        cc[2][0][2][4][2] = -sqrt (20.0 / 49.0);
        cc[2][0][2][4][4] = sqrt (1.0 / 49.0);

        cc[2][1][2][1][2] = sqrt (30.0 / 49.0);
        cc[2][1][2][1][4] = sqrt (40.0 / 49.0);
        cc[2][1][2][2][2] = sqrt (5.0 / 49.0);
        cc[2][1][2][2][4] = sqrt (30.0 / 49.0);
        cc[2][1][2][3][0] = -1.0;
        cc[2][1][2][3][2] = -sqrt (5.0 / 49.0);
        cc[2][1][2][3][4] = sqrt (16.0 / 49.0);
        cc[2][1][2][4][2] = -sqrt (30.0 / 49.0);
        cc[2][1][2][4][4] = sqrt (5.0 / 49.0);

        cc[2][2][2][2][0] = 1.0;
        cc[2][2][2][2][2] = sqrt (20.0 / 49.0);
        cc[2][2][2][2][4] = sqrt (36.0 / 49.0);
        cc[2][2][2][3][2] = sqrt (5.0 / 49.0);
        cc[2][2][2][3][4] = sqrt (30.0 / 49.0);
        cc[2][2][2][4][2] = -sqrt (20.0 / 49.0);
        cc[2][2][2][4][4] = sqrt (15.0 / 49.0);

        cc[2][3][2][3][2] = sqrt (30.0 / 49.0);
        cc[2][3][2][3][4] = sqrt (40.0 / 49.0);
        cc[2][3][2][4][4] = sqrt (35.0 / 49.0);

        cc[2][4][2][4][4] = sqrt (70.0 / 49.0);
    }

    for (l = 0; l < ll; l++)
    {
        for (li = 0; li < lli; li++)
        {
            for (mi = 0; mi < 2 * li + 1; mi++)
            {
                for (lj = 0; lj < li; lj++)
                {
                    for (mj = 0; mj < 2 * lj + 1; mj++)
                    {
                        cc[li][mi][lj][mj][l] = cc[lj][mj][li][mi][l];
                    }           /*end for mj */
                }               /*end for lj */
                for (mj = 0; mj < mi; mj++)
                {
                    cc[li][mi][li][mj][l] = cc[li][mj][li][mi][l];
                }               /*end for mj */
            }                   /*end for mi */
        }                       /*end for li */
    }                           /*end for l */

    for (l = 0; l < 3; l++)
    {
        for (m = 0; m < 5; m++)
        {
            for (k = 0; k < 3; k++)
            {
                for (n = 0; n < 5; n++)
                {
                    for (p = 0; p < 5; p++)
                        cc[l][m][k][n][p] *= sqrt (1.0 / fourPI);
                }
            }
        }
    }

    /*initialize uR and uI */
    for (l = 0; l < 5; l++)
    {
        for (m = 0; m < 9; m++)
        {
            for (k = 0; k < 9; k++)
            {
                uR[l][m][k] = 0.0;
                uI[l][m][k] = 0.0;
            }
        }
    }

    uR[0][0][0] = 1.0;
    uI[0][0][0] = 0.0;

    if (lli > 1)
    {
        y = 1.0 / sqrt (2.0);
        uR[1][0][0] = y;
        uR[1][0][2] = -y;
        uR[1][1][1] = 1.0;
        uI[1][2][0] = y;
        uI[1][2][2] = y;

        uI[2][0][0] = y;
        uI[2][0][4] = -y;
        uR[2][1][1] = y;
        uR[2][1][3] = -y;
        uR[2][2][2] = 1.0;
        uI[2][3][1] = y;
        uI[2][3][3] = y;
        uR[2][4][0] = y;
        uR[2][4][4] = y;
    }
    if (lli > 2)
    {
        f = sqrt (5.0) / 4.0;
        t = sqrt (3.0) / 4.0;
        uR[3][0][0] = f;
        uR[3][0][2] = -t;
        uR[3][0][4] = t;
        uR[3][0][6] = -f;
        uI[3][1][0] = -f;
        uI[3][1][2] = -t;
        uI[3][1][4] = -t;
        uI[3][1][6] = -f;
        uI[3][2][1] = y;
        uI[3][2][5] = -y;
        uR[3][3][3] = 1.0;
        uR[3][4][1] = y;
        uR[3][4][5] = y;
        uI[3][5][0] = -t;
        uI[3][5][2] = f;
        uI[3][5][4] = f;
        uI[3][5][6] = -t;
        uR[3][6][0] = -t;
        uR[3][6][2] = -f;
        uR[3][6][4] = f;
        uR[3][6][6] = t;

        s = sqrt (7.0) / 4.0;
        fs = sqrt (5.0 / 6.0) / 2.0;
        ss = sqrt (7.0 / 6.0) / 2.0;
        fts = sqrt (14.0 / 6.0) / 2.0;
        ts = sqrt (10.0 / 6.0) / 2.0;
        uR[4][0][0] = fs;
        uR[4][0][4] = fts;
        uR[4][0][8] = fs;
        uI[4][1][1] = -0.25;
        uI[4][1][3] = -s;
        uI[4][1][5] = -s;
        uI[4][1][7] = -0.25;
        uR[4][2][1] = -0.25;
        uR[4][2][3] = s;
        uR[4][2][5] = -s;
        uR[4][2][7] = 0.25;
        uR[4][3][2] = -y;
        uR[4][3][6] = -y;
        uI[4][4][0] = y;
        uI[4][4][8] = -y;
        uI[4][5][2] = y;
        uI[4][5][6] = -y;
        uR[4][6][1] = -s;
        uR[4][6][3] = -0.25;
        uR[4][6][5] = 0.25;
        uR[4][6][7] = s;
        uI[4][7][1] = s;
        uI[4][7][3] = -0.25;
        uI[4][7][5] = -0.25;
        uI[4][7][7] = s;
        uR[4][8][0] = -ss;
        uR[4][8][4] = ts;
        uR[4][8][8] = -ss;
    }

    /*initialize aaR and aaI */
    for (l = 0; l < 5; l++)
    {
        for (m = 0; m < 9; m++)
        {
            for (k = 0; k < 3; k++)
            {
                for (n = 0; n < 5; n++)
                {
                    for (p = 0; p < 3; p++)
                    {
                        for (lp = 0; lp < 5; lp++)
                        {
                            aaR[l][m][k][n][p][lp] = 0.0;
                            aaI[l][m][k][n][p][lp] = 0.0;
                        }
                    }
                }
            }
        }
    }


    for (lp = 0; lp < ll; lp++)
    {
        for (l = 0; l < lli; l++)
        {
            for (m = 0; m < 2 * l + 1; m++)
            {
                for (k = 0; k < lli; k++)
                {
                    for (n = 0; n < 2 * k + 1; n++)
                    {
                        for (np = 0; np < 2 * k + 1; np++)
                        {
                            for (p = 0; p < 2 * l + 1; p++)
                            {
                                mp = (p - l) + (np - k) + lp;
                                if ((mp >= 0) && (mp <= 2 * lp))
                                {
                                    aaR[lp][mp][l][m][k][n] =
                                        aaR[lp][mp][l][m][k][n] + (uR[l][m][p] * uR[k][n][np] -
                                                                   uI[l][m][p] * uI[k][n][np]) *
                                        cc[l][p][k][np][lp];
                                    aaI[lp][mp][l][m][k][n] =
                                        aaI[lp][mp][l][m][k][n] + (uR[l][m][p] * uI[k][n][np] +
                                                                   uI[l][m][p] * uR[k][n][np]) *
                                        cc[l][p][k][np][lp];
                                }
                            }
                        }
                    }
                }
            }
        }
    }

/*	for(ik=0;ik<nlx;ik++) {
		for(il=0;il<nlx;il++){
			for(ii=0;ii<mx;ii++) {
				lpl[il][ik][ii]=0;
			}
			lpx[il][ik]=0;
		}
	}
*/
    for (lp = 0; lp < ll; lp++)
    {
        for (mp = 0; mp < 2 * lp + 1; mp++)
        {
            for (l = 0; l < lli; l++)
            {
                for (m = 0; m < 2 * l + 1; m++)
                {
                    for (k = 0; k < lli; k++)
                    {
                        for (n = 0; n < 2 * k + 1; n++)
                        {
                            sumR = 0.0;
                            sumI = 0.0;
                            for (mpp = 0; mpp < 2 * lp + 1; mpp++)
                            {
                                sumR =
                                    sumR + uR[lp][mp][mpp] * aaR[lp][mpp][l][m][k][n] +
                                    uI[lp][mp][mpp] * aaI[lp][mpp][l][m][k][n];
                                sumI =
                                    sumI + uR[lp][mp][mpp] * aaI[lp][mpp][l][m][k][n] -
                                    uI[lp][mp][mpp] * aaR[lp][mpp][l][m][k][n];
                            }
                            if (sqrt (sumR * sumR + sumI * sumI) > 0.001)
                            {
                                il = l * l + m;
                                ik = k * k + n;
                                ilp = lp * lp + mp;
                                if (fabs (sumI) > 0.001)
                                {
                                    printf
                                        ("error!! l=%d m=%d k=%d n=%d lp=%d mp=%d il=%d ik=%d ilp=%d sumI=%f\n",
                                         l, m, k, n, lp, mp, il, ik, ilp, sumI);
                                    error_handler ("there is something wrong about the sumI");
                                }
                                ap[ilp][il][ik] = sumR;
                                lpl[il][ik][lpx[il][ik]] = ilp;
/*if(pct.gridpe==0) printf("ap[%d][%d][%d]=%f lpx[%d][%d]=%d lpl[%d][%d][%d]=%d\n",ilp,il,ik,ap[ilp][il][ik],il,ik,lpx[il][ik],il,ik,lpx[il][ik],ilp);*/
                                lpx[il][ik] = lpx[il][ik] + 1;
                            }   /*end for if */
                        }       /*end for n */
                    }
                }
            }
        }
    }
}
