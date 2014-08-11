#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <sys/time.h>
#include <assert.h>
#include <pthread.h>
#include "string.h"


void gridart_preprocessing(int num_grid, int num_cols, float center, float *mov, float *gridx, float *gridy){

    *mov = (float)num_cols/2. - center;
    if(*mov - ceil(*mov) < 1e-6)
        *mov += 1e-6;

    int i;
    for(i=0; i<= num_grid; i++){
        gridx[i] = -num_grid/2. + i;
        gridy[i] = -num_grid/2. + i;
    }

}

bool calc_quadrant(float theta_q){
    const double PI = 3.141592653589793238462;
    bool quadrant;
    if ((theta_q >= 0 && theta_q < PI/2) ||
            (theta_q >= PI && theta_q < 3*PI/2)) {
        quadrant = true;
    } else {
        quadrant = false;
    }
    return quadrant;
}

void calc_coordinates(
        int num_grid,
        float xi, float yi,
        float sinq, float cosq,
        float *gridx, float *gridy,
        float *coordx, float *coordy)
{
    float srcx, srcy, detx, dety;
    float slope, islope;
    int n;

    // Find the corresponding source and
    // detector locations for a given line
    // trajectory of a projection (Projection
    // is specified by sinp and cosp).   
    srcx = xi*cosq-yi*sinq;
    srcy = xi*sinq+yi*cosq;
    detx = -xi*cosq-yi*sinq;
    dety = -xi*sinq+yi*cosq;

    // Find the intersection points of the 
    // line connecting the source and the detector
    // points with the reconstruction grid. The 
    // intersection points are then defined as: 
    // (coordx, gridy) and (gridx, coordy)
    slope = (srcy-dety)/(srcx-detx);
    islope = 1/slope;
    for (n = 0; n <= num_grid; n++) {
        coordx[n] = islope*(gridy[n]-srcy)+srcx;
        coordy[n] = slope*(gridx[n]-srcx)+srcy;
    }
}

void merge_trim_coord(
        int num_grid,
        float *coordx, float *coordy, 
        float *gridx, float* gridy, 
        int *alen, int *blen, 
        float *ax, float *ay, 
        float *bx, float *by)
{
    int n;
    // Merge the (coordx, gridy) and (gridx, coordy)
    // on a single array of points (ax, ay) and trim
    // the coordinates that are outside the
    // reconstruction grid. 
    *alen = 0;
    *blen = 0;
    for (n = 0; n <= num_grid; n++) {
        if (coordx[n] > gridx[0]) {
            if (coordx[n] < gridx[num_grid]) {
                ax[*alen] = coordx[n];
                ay[*alen] = gridy[n];
                (*alen)++;
            }
        }
        if (coordy[n] > gridy[0]) {
            if (coordy[n] < gridy[num_grid]) {
                bx[*blen] = gridx[n];
                by[*blen] = coordy[n];
                (*blen)++;
            }
        }
    }
}

void sort_inters(
        int ind_cond, 
        int alen, int blen, 
        float *ax, float *ay, 
        float *bx, float *by, 
        float *coorx, float *coory)
{
    int i=0, j=0, k=0;
    int a_ind;
    while (i < alen && j < blen)
    {
        a_ind = (ind_cond) ? i : (alen-1-i);
        if (ax[a_ind] < bx[j]) {
            coorx[k] = ax[a_ind];
            coory[k] = ay[a_ind];
            i++;
            k++;
        } else {
            coorx[k] = bx[j];
            coory[k] = by[j];
            j++;
            k++;
        }
    }
    while (i < alen) {
        a_ind = (ind_cond) ? i : (alen-1-i);
        coorx[k] = ax[a_ind];
        coory[k] = ay[a_ind];
        i++;
        k++;
    }
    while (j < blen) {
        coorx[k] = bx[j];
        coory[k] = by[j];
        j++;
        k++;
    }
}

void calculate_dist_len(int len, int num_grid, float *coorx, float *coory, float *leng, float *leng2, int *indi){
    int n, x1, x2, i1, i2;
    float diffx, diffy, midx, midy;
    int indx, indy;

    for (n = 0; n < len-1; n++) {
        diffx = coorx[n+1] - coorx[n];
        diffy = coory[n+1] - coory[n];
        leng2[n] = diffx * diffx + diffy * diffy;
        leng[n] = sqrt(leng2[n]);
        midx = (coorx[n+1] + coorx[n])/2;
        midy = (coory[n+1] + coory[n])/2;
        x1 = midx + num_grid/2.;
        x2 = midy + num_grid/2.;
        i1 = (int)(midx + num_grid/2.);
        i2 = (int)(midy + num_grid/2.);
        indx = i1 - (i1 > x1);
        indy = i2 - (i2 > x2);
        indi[n] = (indx+(indy*num_grid));
    }
}

float calc_simdata(int sliceID, int num_grids, int len, int *indi, float *leng, float *recon){
    float simdata = 0.0;
    int n;

    int index = sliceID*num_grids*num_grids;
    for (n = 0; n < len-1; n++) {
        //simdata += recon_data[indi[n]] * leng[n];
        simdata += recon[indi[n]+index]*leng[n];
    }

    return simdata;
}

void update_recon(float simdata, float *data, int *indi, int num_grids, float *leng2, float *leng, int io, int sliceID, int len, float *recon){
    float upd;
    int n;

    // Note: The indices (indi) and the corresponding
    // weights (leng) are the same for all slices. So,
    // there is no need to calculate them for each slice.
    float a2 = 0.0;
    for (n = 0; n < len-1; n++) {
        a2 += leng2[n];
    }

    int index = sliceID*num_grids*num_grids;
    upd = (data[io]-simdata)/a2;
    for (n = 0; n < len-1; n++) {
        //accumulate_float(sliceID, indi[n], upd*leng[n]);
        recon[indi[n]+index] += upd*leng[n];
    }

}

void myart(float* data, float* theta, float center,
        int num_projs, int tot_num_slices,
        int num_cols, int num_grid, int num_iter, 
        float *recon) // Recon is reductiom object
{
    float* gridx = (float *)malloc((num_grid+1)*sizeof(float));
    float* gridy = (float *)malloc((num_grid+1) * sizeof(float));
    float* coordx = (float *)malloc((num_grid+1) * sizeof(float));
    float* coordy = (float *)malloc((num_grid+1) * sizeof(float));
    float* ax = (float *)malloc((num_grid+1) * sizeof(float));
    float* ay = (float *)malloc((num_grid+1) * sizeof(float));
    float* bx = (float *)malloc((num_grid+1) * sizeof(float));
    float* by = (float *)malloc((num_grid+1) * sizeof(float));
    float* coorx = (float *)malloc((2*num_grid) * sizeof(float));
    float* coory = (float *)malloc((2*num_grid) * sizeof(float));
    float* leng = (float *)malloc((2*num_grid) * sizeof(float));
    float* leng2 = (float *)malloc((2*num_grid) * sizeof(float));
    int* indi = (int *)malloc((2*num_grid) * sizeof(int));

    assert( coordx != NULL && coordy != NULL &&
            ax != NULL && ay != NULL && by != NULL && bx != NULL &&
            coorx != NULL && coory != NULL &&
            leng != NULL && leng2 != NULL && indi != NULL);

    int m, q, kk;
    int quadrant;
    float sinq, cosq;

    float xi, yi;
    int alen, blen, len;
    int io;
    float simdata;
    float mov;

    int iter;

    gridart_preprocessing(num_grid, num_cols, center, &mov , gridx, gridy);

    for(iter=0; iter<num_iter; iter++){

        // For each slice
        for (kk = 0; kk < tot_num_slices; kk++) {

            // For each projection angle
            for (q = 0; q < num_projs; q++) {
                //printf("Iter=%d; Slice=%d; Projection=%d; Cols=%d; Center=%.10f; Mov=%.10f\n", iter, kk, q, num_cols, center, mov);

                // Calculate the sin and cos values 
                // of the projection angle and find
                // at which quadrant on the cartesian grid.
                // (Can be calculated independently)
                float theta_q = theta[q];
                quadrant = calc_quadrant(theta_q);

                //calculate_sincos(theta_q, &sinq, &cosq);
                // sincosf(theta_q, &sinq, &cosq);
                sinq = sinf(theta_q);
                cosq = cosf(theta_q);

                //printf("kk=%d; q=%d; theta_q=%.10f; quadrant=%d; sinq=%.10f; cosq=%.10f; center=%.10f; mov=%.10f\n", 
                //        kk, q, theta_q, quadrant, sinq, cosq, center, mov);
                for (m = 0; m < num_cols; m++) {
                    // Calculate coordinates
                    xi = -1e6;
                    yi = -(num_cols-1)/2. + m + mov;
                    calc_coordinates(num_grid, xi, yi, sinq, cosq, gridx, gridy, 
                            coordx, coordy); // Outputs coordx and coordy

                    // Merge the (coordx, gridy) and (gridx, coordy)
                    // Output alen and after
                    merge_trim_coord(num_grid, coordx, coordy, gridx, gridy, &alen, &blen, ax, ay, bx, by);

                    // Sort the array of intersection points (ax, ay)
                    // The new sorted intersection points are 
                    // stored in (coorx, coory).
                    // if quadrant=1 then a_ind = i; if 0 then a_ind = (alen-1-i)
                    sort_inters(quadrant, alen, blen, ax, ay, bx, by, coorx, coory);

                    // Calculate the distances (leng) between the 
                    // intersection points (coorx, coory). Find 
                    // the indices of the pixels on the  
                    // reconstruction grid (ind_recon).
                    len = alen+blen;
                    calculate_dist_len(len, num_grid, coorx, coory, leng, leng2, indi);

                    //*******************************************************
                    // Below is for updating the reconstruction grid and 
                    // is algorithm specific part. 

                    // Calculate simdata 
                    simdata = calc_simdata(kk, num_grid, len, indi, leng, recon);

                    // Update recon
                    io = m+(kk*num_cols)+q*(tot_num_slices*num_cols);
                    update_recon(simdata, data, indi, num_grid, leng2, leng, io, kk, len, recon);
                    //*******************************************************
                }
            }
        }
    }

    free(coordx);
    free(coordy);
    free(ax);
    free(ay);
    free(bx);
    free(by);
    free(coorx);
    free(coory);
    free(leng);
    free(leng2);
    free(indi);
}


