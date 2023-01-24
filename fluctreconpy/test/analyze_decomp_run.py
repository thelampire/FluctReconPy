#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 18:24:11 2023

@author: mlampert
"""
import numpy as np
import pickle

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

def analyze_decomp_run(filename='tmp/test_fluct_decomp_all_BS_0.5_5.0_NL_0.0_0.5_SR_10.0_10.0.pickle',
                       pdf=True):

    if pdf:
        pdf_pages=PdfPages(filename+'.pdf')
    ymax=2
    imin=1
    number=20
    fnum=20

    with open(filename,'rb') as f: results=pickle.load(f)

    spatial_resolution=results['spatial_resolution']
    noise_level=results['noise_level']
    blob_size=results['blob_size']
    fwhm_calc=results['fwhm_calc']/(2*np.sqrt(2*np.log(2)))
    fwhm_orig=results['fwhm_orig']/(2*np.sqrt(2*np.log(2)))
    pos_calc=results['pos_calc']
    pos_orig=results['pos_orig']
    density_calc=results['density_calc']
    density_samp=results['density_samp']
    nwin_time=results['nwin_time']

    fig,ax=plt.subplots(figsize=(8.5/2.54,8.5/2.54))
    for j in range(len(spatial_resolution)):
        for i in range(imin,len(blob_size)):
            y=(np.abs(fwhm_calc-fwhm_orig)/(spatial_resolution[j]))[:,i,j,:,0]

            ax.errorbar(noise_level,
                        np.sum(y,axis=1)/nwin_time,
                        yerr=np.sqrt(np.var(y,axis=1)/(nwin_time-1)),
                        label=str(blob_size[i]*spatial_resolution[j])+'mm',
                        fmt='-o',
                        )
            # ax.set_ylim([0,2])
            ax.set_xlabel('Relative noise amplitude')
            ax.set_ylabel('Radial FWHM error')
            ax.set_title('Spat res:'+str(spatial_resolution[j])+'mm')

    plt.legend()
    plt.tight_layout(pad=0.1)
    if pdf:
        pdf_pages.savefig()

    fig,ax=plt.subplots(figsize=(8.5/2.54,8.5/2.54))
    for j in range(len(spatial_resolution)):
        for i in range(imin,len(blob_size)):
            y=(np.abs(fwhm_calc-fwhm_orig)/(spatial_resolution[j]))[:,i,j,:,1]

            ax.errorbar(noise_level,
                     np.sum(y,axis=1)/nwin_time,
                     yerr=np.sqrt(np.var(y,axis=1)/(nwin_time-1)),
                     label=str(blob_size[i]*spatial_resolution[j])+'mm',
                     fmt='-o',)
            # ax.set_ylim([0,2])
            ax.set_xlabel('Relativenoise amplitude')
            ax.set_ylabel('Vertical FWHM error')
            ax.set_title('Spat res:'+str(spatial_resolution[j])+'mm')
    plt.legend()
    plt.tight_layout(pad=0.1)
    if pdf:
        pdf_pages.savefig()

    fig,ax=plt.subplots(figsize=(8.5/2.54,8.5/2.54))
    for j in range(len(spatial_resolution)):
        for i in range(imin,len(blob_size)):
            y=(np.abs(pos_calc-pos_orig)/(spatial_resolution[j]))[:,i,j,:,0]
            ax.errorbar(noise_level,
                        np.sum(y,axis=1)/nwin_time,
                        yerr=np.sqrt(np.var(y,axis=1)/(nwin_time-1)),
                        label=str(blob_size[i]*spatial_resolution[j])+'mm',
                        fmt='-o')
            # ax.set_ylim([0,2])
            ax.set_xlabel('Relative noise amplitude')
            ax.set_ylabel('Radial position error')
            ax.set_title('Spat res:'+str(spatial_resolution[j])+'mm')
    plt.legend()
    plt.tight_layout(pad=0.1)
    if pdf:
        pdf_pages.savefig()

    fig,ax=plt.subplots(figsize=(8.5/2.54,8.5/2.54))
    for j in range(len(spatial_resolution)):
        for i in range(imin,len(blob_size)):
            y=(np.abs(pos_calc-pos_orig)/(spatial_resolution[j]))[:,i,j,:,1]

            ax.errorbar(noise_level,
                        np.sum(y,axis=1)/nwin_time,
                        yerr=np.sqrt(np.var(y,axis=1)/(nwin_time-1)),
                        label=str(blob_size[i]*spatial_resolution[j])+'mm',
                        fmt='-o',
                        )
            # ax.set_ylim([0,2])
            ax.set_xlabel('Relative noise amplitude')
            ax.set_ylabel('Vertical position error')
            ax.set_title('Spat res:'+str(spatial_resolution[j])+'mm')
    plt.legend()
    plt.tight_layout(pad=0.1)

    if pdf:
        pdf_pages.savefig()

    fig,ax=plt.subplots(figsize=(8.5/2.54,8.5/2.54))
    for j in range(len(spatial_resolution)):
        for i in range(imin,len(blob_size)):
            y=(np.abs(1-density_calc*fwhm_calc[:,:,:,:,0]*fwhm_calc[:,:,:,:,1]/(density_samp*fwhm_orig[:,:,:,:,0]*fwhm_orig[:,:,:,:,1])))[:,i,j,:]
            ax.errorbar(noise_level,
                        np.sum(y,axis=1)/nwin_time,
                        yerr=np.sqrt(np.var(y,axis=1)/(nwin_time-1)),
                        label=str(blob_size[i]*spatial_resolution[j])+'mm',
                        fmt='-o')
            # ax.set_ylim([0,2])
            ax.set_xlabel('Relative noise amplitude')
            ax.set_ylabel('Density reconstruction error')
            ax.set_title('Spat res:'+str(spatial_resolution[j])+'mm')
    plt.legend()
    plt.tight_layout(pad=0.1)

    if pdf:
        pdf_pages.savefig()
        pdf_pages.close()


