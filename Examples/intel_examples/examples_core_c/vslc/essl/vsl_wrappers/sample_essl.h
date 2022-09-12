/*******************************************************************************
* Copyright 2003-2020 Intel Corporation.
*
* This software and the related documents are Intel copyrighted  materials,  and
* your use of  them is  governed by the  express license  under which  they were
* provided to you (License).  Unless the License provides otherwise, you may not
* use, modify, copy, publish, distribute,  disclose or transmit this software or
* the related documents without Intel's prior written permission.
*
* This software and the related documents  are provided as  is,  with no express
* or implied  warranties,  other  than those  that are  expressly stated  in the
* License.
*******************************************************************************/

#ifndef _SAMPLE_ESSL_H_
#define _SAMPLE_ESSL_H_

#ifdef __cplusplus
extern "C" {
#endif

/*
// (A) Convolution/correlation via Fourier transfrom
*/

void sconf(
    int init,
    float h[], int inc1h,
    float x[], int inc1x, int inc2x,
    float y[], int inc1y, int inc2y,
    int nh, int nx, int m, int iy0, int ny,
    void* aux1, int naux1, void* aux2, int naux2);

void scorf(
    int init,
    float h[], int inc1h,
    float x[], int inc1x, int inc2x,
    float y[], int inc1y, int inc2y,
    int nh, int nx, int m, int iy0, int ny,
    void* aux1, int naux1, void* aux2, int naux2);

/*
// (B) Convolution/correlation via direct method
*/

void scond(
    float h[], int inch,
    float x[], int incx,
    float y[], int incy,
    int nh, int nx, int iy0, int ny);

void scord(
    float h[], int inch,
    float x[], int incx,
    float y[], int incy,
    int nh, int nx, int iy0, int ny);

/*
// (C) Convolution/correlation via direct method
//     with decimation of the output sequence
*/

/* (C.1) single precision */

void sdcon(
    float h[], int inch,
    float x[], int incx,
    float y[], int incy,
    int nh, int nx, int iy0, int ny, int id);

void sdcor(
    float h[], int inch,
    float x[], int incx,
    float y[], int incy,
    int nh, int nx, int iy0, int ny, int id);

/* (C.2) single precision data, but double precision calculations */

void sddcon(
    float h[], int inch,
    float x[], int incx,
    float y[], int incy,
    int nh, int nx, int iy0, int ny, int id);

void sddcor(
    float h[], int inch,
    float x[], int incx,
    float y[], int incy,
    int nh, int nx, int iy0, int ny, int id);

/* (C.3) double precision data and calculations */

void ddcon(
    double h[], int inch,
    double x[], int incx,
    double y[], int incy,
    int nh, int nx, int iy0, int ny, int id);

void ddcor(
    double h[], int inch,
    double x[], int incx,
    double y[], int incy,
    int nh, int nx, int iy0, int ny, int id);

#ifdef __cplusplus
}
#endif

#endif/*_SAMPLE_ESSL_H_*/
