#!/usr/bin/python -tt

from __future__ import division
import sys
import re
import numpy as np
from scipy.optimize import curve_fit
from scipy.optimize import fmin
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
import struct as stc

# this is a random change!! KNUT KNUT KNUT!!

"""Positron track analysis
author: Shih-ying Clare Huang
date: April 22nd, 2014

    - version 2.0
    - eliminate unnecessary functions for a clean python script for the task of analyzing positron range effect from MR
"""

def ReadFile(fname):
    """
    read the colume data from fname and return the data array
    (x,y,z) are in unit of "cm"
    """
    print 'Read the file: ',fname
    data = np.loadtxt(fname)
    #print data
    return data

def Plot3DScatter(xx,yy,zz):
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    ax.scatter(xx,yy,zz,c='r',marker='o')

def GetHisto(bin,data):
    hist,hbin=np.histogram(data,bin)
    #print data.shape,bin.shape,hist.shape,hbin.shape
    #print hist
    #print hbin
    center = (hbin[:-1] + hbin[1:])/2
    fig = plt.figure()
    ax = fig.add_subplot(111)
    width = 0.05
    ax.bar(center,hist,width,color='red')

def IsInBound(val,low,high):
    """
    check to make sure VAL in within the range (low,high]
    """
    if val < low:
        #print 'less than low'
        return low
    elif val >= high:
        #print 'more than high'
        return high-1
    else:
        return val

def FindCofMass(vol):
    """
    determine the center of mass for a given 3D volume
    i.e.X_center = sum(ix[...]*vol[...])/sum(vol)
        same for Y_center and Z_center
    """
    xlen,ylen,zlen = vol.shape
    xx = np.arange(0,xlen)
    yy = np.arange(0,ylen)
    zz = np.arange(0,zlen)
    ix,iy,iz = np.where(vol > 0)
    wvol = vol[vol>0]
    
    sumwvol = sum(wvol)
    X_center = round(sum(ix*wvol)/sumwvol)
    Y_center = round(sum(iy*wvol)/sumwvol)
    Z_center = round(sum(iz*wvol)/sumwvol)
    
    return X_center,Y_center,Z_center

def Position2Vol(fdir,maxEvent,run1,run2,dxyz,xyz0,xyz1):
    """"
    (1) read the (x,y,z)'s from PositronPositon.dat for a given positron state
    (2) digitize (x,y,z) into voxel dimention specified by dxyz (voxel dimension), xyz0 (starting position of an image plane), xyz1 (ending position of an image plane)
    (3) accumulate the existing position in a 3d array, vol
    (4) return the 3D arrays, vol0 and vol1, and its corresponding list of physical dimension in x,y,z direction
    """
    dx = dxyz[0]
    dy = dxyz[1]
    dz = dxyz[2]
    x0 = xyz0[0]
    y0 = xyz0[1]
    z0 = xyz0[2]
    x1 = xyz1[0]
    y1 = xyz1[1]
    z1 = xyz1[2]
    xx = np.arange(x0,x1+dx,dx)
    yy = np.arange(y0,y1+dy,dy)
    zz = np.arange(z0,z1+dz,dz)
    xxlen =  len(xx)
    yylen = len(yy)
    zzlen = len(zz)
    vol0 = np.zeros((xxlen,yylen,zzlen),dtype=int)
    vol1 = np.zeros((xxlen,yylen,zzlen),dtype=int)
    
    for k in range(run1,run2+1):
        fname = fdir+'/Run'+str(k)+'/PositronPosition.dat'
        #print fname
        data = ReadFile(fname)
        #print data.shape[0]
        #print data.shape[1]
        if maxEvent > data.shape[0]:
            maxNevt = data.shape[0]
        else:
            maxNevt = maxEvent
        xpos0 = data[0:(maxNevt-1),0]
        ypos0 = data[0:(maxNevt-1),1]
        zpos0 = data[0:(maxNevt-1),2]
        xpos1 = data[0:(maxNevt-1),3]
        ypos1 = data[0:(maxNevt-1),4]
        zpos1 = data[0:(maxNevt-1),5]
        #print np.min(xpos1)
        #print np.max(xpos1)
        #print np.min(ypos1)
        #print np.max(ypos1)
        #print np.min(zpos1)
        #print np.max(zpos1)
        
        # Get histogram of the input position data
        #GetHisto(xx,xpos1)
        #GetHisto(yy,ypos1)
        #GetHisto(zz,zpos1)
        
        # 3D scatter plot of input positions
        #Plot3DScatter(xpos1,ypos1,zpos1)
        
        for i in range(len(xpos0)):
            ix0 = IsInBound(round((xpos0[i]-x0)/dx),0,xxlen)
            iy0 = IsInBound(round((ypos0[i]-y0)/dy),0,yylen)
            iz0 = IsInBound(round((zpos0[i]-z0)/dz),0,zzlen)
            ix1 = IsInBound(round((xpos1[i]-x0)/dx),0,xxlen)
            iy1 = IsInBound(round((ypos1[i]-y0)/dy),0,yylen)
            iz1 = IsInBound(round((zpos1[i]-z0)/dz),0,zzlen)
            vol0[iy0,ix0,iz0] += 1
            vol1[iy1,ix1,iz1] += 1
    #print 'Sum of vol0: ',str(np.sum(vol0))
    #print 'Sum of vol1: ',str(np.sum(vol1))
    return vol0,vol1,xx,yy,zz

def Get2DPlane(vol,idir,indx,fname):
    """
    vol: 3D array
    idir: the direction of the image plane
    index: the index of the image plane
    im: return the 2D image in a desired plane
    """
    if idir == 0:
        im = vol[indx,:,:]
    elif idir == 1:
        im = vol[:,indx,:]
    else:
        im = vol[:,:,indx]
    # object oriented interface of matplotlib
    fig = plt.figure()
    ax = fig.add_subplot(111)
    # specificy the min and max of image display with vmin, vmax
    ax.imshow(im,cmap=plt.cm.gray,vmin=0,vmax=10,interpolation='nearest')
    fig.savefig(fname)
    return im

def ShowLineProfile(vol,idir,indx,xx,yy,zz,isplot):
    """
    Plot a line profile through a given image dimension at indx1, indx2
    """
    if idir == 0:
        xdata = yy
        ydata = np.sum(vol[indx,:,:],axis=0)
    elif idir == 1:
        xdata = xx
        ydata = np.sum(vol[:,indx,:],axis=0)
    else:
        xdata = zz
        ydata = np.sum(vol[:,:,indx],axis=0)
    if isplot == 1:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(xdata,ydata,'ro')
        #plt.plot(xdata,ydata,'ro')
    return xdata, ydata

def GetHalfLineProfile(im,xdata,iaxis,isplot,issavefile,fname):
    """
    Get the line profile by summing 'im' in axis=0 direction
    """
    ydata = np.sum(im,axis=iaxis)
    
    # only get the part with x > 0
    xdata_new = xdata[xdata >= 0.0]
    ydata_new = ydata[xdata >=0.0]
    
    if isplot == 1:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(xdata,ydata,'ro')
        
    if issavefile == 1:
        np.savetxt(fname,(xdata_new,ydata_new))
    return xdata_new,ydata_new

def GetR2(data,fitdata):
    fres = np.sum((data-fitdata)**2)
    ftot = np.sum((data-data.mean())**2)
    R2 = 1 - fres/ftot
    return R2
    
def gauss(x,*p):
    A,mu,sigma = p
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def PlotGaussFit(xdata,ydata,f0):
    coeff,var_matrix = curve_fit(gauss,xdata,ydata,p0=f0,maxfev=1000000)
    ydata_fit = gauss(xdata,*coeff)
    R2 = GetR2(ydata,ydata_fit)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.plot(xdata,ydata,'ro',label='Original data')
    ax.plot(xdata,ydata_fit,label='Fitted data')
    
    # load the output data
    outdata = np.zeros(len(coeff)+1)
    outdata[:-1] = coeff
    outdata[-1] = R2
    return outdata

def BiGaussian(x,a0,mu0,sig0,a1,mu1):
    return a0*np.exp(-((x-mu0)/sig0)**2) + a1*np.exp(-((x-mu1)/sig0)**2)

def PlotBiGaussFit(xdata,ydata,f0,isplot,fname):
    ymax = np.max(ydata)
    ydata_norm = ydata/ymax
    print 'initial values: ', f0
    coeff,cov = curve_fit(BiGaussian,xdata,ydata_norm,f0,maxfev=1000000)
    ydata_fit = ymax*BiGaussian(xdata,*coeff)
    R2 = GetR2(ydata,ydata_fit)
    if isplot == 1:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(xdata,ydata,'ro',label='Original data')
        ax.plot(xdata,ydata_fit,label='Fitted data')
        fig.savefig(fname)
    
    #compile all fit data
    #outdata = np.zeros(len(coeff)+1)
    #outdata[:-1] = coeff
    #outdata[-1] = R2
    #print outdata
    return coeff,R2
    
def BiExpo(x,a0,k0,a1,k1):
    return a0*np.exp(-k0*x) + a0*np.exp(-k1*x)
    
def PlotBiExpoFit(xdata,ydata,f0,isplot,fname):
    ymax = np.max(ydata)
    ydata_norm = ydata/ymax
    print 'initial values: ', f0
    # the below statement didn't help with the fit...
    #coeff,cov = curve_fit(BiExpo,xdata,ydata,maxfev=1000000,diag=(1./np.mean(xdata),1./np.max(ydata)))
    coeff,cov = curve_fit(BiExpo,xdata,ydata_norm,f0,maxfev=1000000)
    ydata_fit = ymax*BiExpo(xdata,*coeff)
    R2 = GetR2(ydata,ydata_fit)
    if isplot == 1:
        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(xdata,ydata,'ro',label='Original data')
        ax.plot(xdata,ydata_fit,label='Fitted data')
        fig.savefig(fname)

    return coeff,R2

def flatten3DArray( inputArray ):
    """
    returns a flattened 3D array, suitable for writing as a binary file
    source: CT_data_reader_2012.06.27.py from 
    """
    dimensions = inputArray.shape
    numberOfElements = dimensions[ 0 ] * dimensions[ 1 ] * dimensions[ 2 ]
    flattenedArray = np.zeros( numberOfElements ) 
    currentPixel = 0
    for iz in range( dimensions[ 2 ] ):
        for iy in range( dimensions[ 0 ] ):
            for ix in range( dimensions[ 1 ] ):
                flattenedArray[ currentPixel ] = inputArray[ iy, ix, iz]
                currentPixel += 1
    return( flattenedArray )

def SaveBinaryFile(intputArray,fname):
    # save into binary file
    outputType = 'float64'
    outputExtension = '.f64le'
    lefile = fname + outputExtension
    flatteninputArray = flatten3DArray(intputArray.astype( outputType ))
    flatteninputArray.tofile(lefile)
    
def Objfunc_BiGauss(xinput,coeff,pscale):
    a0,mu0,sig0,a1,mu1 = coeff
    fmax = BiGaussian(0.1,a0,mu0,sig0,a1,mu1)
    f = BiGaussian(xinput,a0,mu0,sig0,a1,mu1)
    err = abs(f - pscale*fmax)
    return err

def GetBiGaussFWM(coeff,pscale):
    """
    Determine the full-width [] maximum of a Bi-Gaussian function with coeff
    full-width HALF maximum: pscale = 0.5
    full-width TENTH maximum: pscale = 0.1
    """
    x0 = 0.0
    x_sol = fmin(Objfunc_BiGauss,x0,args=(coeff,pscale),xtol=0.00000001)
    print x_sol
    return x_sol

def Objfunc_BiExpo(xinput,coeff,pscale):
    a0,k0,a1,k1= coeff
    fmax = BiExpo(0.0001,a0,k0,a1,k1)
    f = BiExpo(xinput,a0,k0,a1,k1)
    err = abs(f - pscale*fmax)
    return err
    
def GetBiExpoFWM(coeff,pscale):
    """
    Determine the full-width [] maximum of a Bi-Gaussian function with coeff
    full-width HALF maximum: pscale = 0.5
    full-width TENTH maximum: pscale = 0.1
    """
    x0 = 0.0
    x_sol = fmin(Objfunc_BiExpo,x0,args=(coeff,pscale),xtol=0.00000001)
    print x_sol
    return x_sol


def main():
    # This command-line parsing code is provided.
    # Make a list of command line arguments, omitting the [0] element
    # which is the script itself.
    
    args = sys.argv[1:]
    if not args:
        print 'Usage: fildir run1 run1 maxEvent'
        sys.exit(1)
    fdir = args[0]
    run1 = int(args[1])
    run2 = int(args[2])
    maxEvent = int(args[3])
    
    dxyz = [0.05,0.05,0.05]
    xyz0 = [-5.0,-5.0,-5.0]
    xyz1 = [5.0,5.0,5.0]
    
    analysisdir = fdir + '/Analysis'
    
    # get the 3D volume
    vol0,vol1,xx,yy,zz = Position2Vol(fdir,maxEvent,run1,run2,dxyz,xyz0,xyz1)
    
    # for testing purpose
    #fname7 = analysisdir + '/Vol_PostPosition_Run' + str(run1) + '_' + str(run2) + '.npy'
    #vol1 = np.load(fname7)
    
    dx = dxyz[0]
    dy = dxyz[1]
    dz = dxyz[2]
    x0 = xyz0[0]
    y0 = xyz0[1]
    z0 = xyz0[2]
    x1 = xyz1[0]
    y1 = xyz1[1]
    z1 = xyz1[2]
    xx = np.arange(x0,x1+dx,dx)
    yy = np.arange(y0,y1+dy,dy)
    zz = np.arange(z0,z1+dz,dz)
    
    
    # find center of mass in the volume
    xc,yc,zc = FindCofMass(vol1)
    #print xc,yc,zc
    
    fname1 = analysisdir + '/CurveFit_Run' + str(run1) + '_' + str(run2) + '.txt'
    combodata = np.zeros((3,8))
    leidir = [2,2,1]
    leiaxis = [0,1,0]
    for idx in range(len(leidir)):
        """
        # plot the 2D plane
        if leidir[idx] == 0:
            # im shape: ix-by-iz
            indx = xc
            fname2 = analysisdir + '/XPositronPostPosition_Run' + str(run1) + '_' + str(run2) + '.pdf'
            fname3 = analysisdir + '/XLineProfileFit_PostPosition_Run' + str(run1) + '_' + str(run2) + '.pdf'
            iaxis = 1
            # sum over the 1st dim (= iz) ==> profile axis = xx
            prfHordata = xx
            #print zz
        elif leidir[idx] == 1:
            # im shape: iy-by-iz
            indx = yc
            fname2 = analysisdir + '/ZPositronPostPosition_Run' + str(run1) + '_' + str(run2) + '.pdf'
            fname3 = analysisdir + '/ZLineProfileFit_PostPosition_Run' + str(run1) + '_' + str(run2) + '.pdf'
            # sum over 0th dim (=iy) ==> profile axis = zz, it's different from the case of idir == 0
            # since it's summing over y-dir instead of x-dir
            iaxis = 0
            prfHordata = zz
        else:
            # im shape: iy-by-ix 
            indx = zc  
            fname2 = analysisdir + '/YPositronPostPosition_Run' + str(run1) + '_' + str(run2) +'.pdf'
            fname3 = analysisdir + '/YLineProfileFit_PostPosition_Run' + str(run1) + '_' + str(run2) + '.pdf'
            # sum over 1st dim (=ix) ==> profile axis = yy
            iaxis = 1
            prfHordata = yy
        """
        # plot the 2D plane
        if leidir[idx] == 2 and leiaxis[idx] == 0:
            # im shape: iy-by-ix
            indx = zc
            fname2 = analysisdir + '/XPositronPostPosition_Run' + str(run1) + '_' + str(run2) + '.pdf'
            fname3 = analysisdir + '/XLineProfileFit_PostPosition_Run' + str(run1) + '_' + str(run2) + '.pdf'
            fname4 = analysisdir + '/XLSFdata_Run' + str(run1) + '_' + str(run2)+ '.txt'
            # sum over the 1st dim (= iz) ==> profile axis = xx
            prfHordata = xx
        elif leidir[idx] == 2 and leiaxis[idx] == 1:
            # im shape: iy-by-ix
            indx = zc
            fname2 = analysisdir + '/YPositronPostPosition_Run' + str(run1) + '_' + str(run2) + '.pdf'
            fname3 = analysisdir + '/YLineProfileFit_PostPosition_Run' + str(run1) + '_' + str(run2) + '.pdf'
            fname4 = analysisdir + '/YLSFdata_Run' + str(run1) + '_' + str(run2)+ '.txt'
            prfHordata = yy
        elif leidir[idx] == 1 and leiaxis[idx] == 0:
            # im shape: iy-by-iz 
            indx = yc  
            fname2 = analysisdir + '/ZPositronPostPosition_Run' + str(run1) + '_' + str(run2) +'.pdf'
            fname3 = analysisdir + '/ZLineProfileFit_PostPosition_Run' + str(run1) + '_' + str(run2) + '.pdf'
            fname4 = analysisdir + '/ZLSFdata_Run' + str(run1) + '_' + str(run2)+ '.txt'
            # sum over 1st dim (=ix) ==> profile axis = yy
            prfHordata = zz
            
        # get 2D plane of interest
        im1 = Get2DPlane(vol1,leidir[idx],indx,fname2)
        
        # plot 1D plane
        xdata1,ydata1 = GetHalfLineProfile(im1,prfHordata,leiaxis[idx],0,1,fname4)
        
        # fit the half profile with BiGaussian fit
        """
        #f1 = [600.0,-3.0,0.1,0.8,1.0]
        f1 = [100.0,-10.0,0.05,0.1,0.1]
        fitcoeff,r2 = PlotBiGaussFit(xdata1,ydata1,f1,1,fname3)
        
        # Calculate the FWHM of the fitted Bi-Gaussian function
        fwhm = GetBiGaussFWM(fitcoeff,0.5)
        fwtm = GetBiGaussFWM(fitcoeff,0.1)
        """
        
        f1 = [0.0,0.0,0.0,0.0]
        fitcoeff,r2 = PlotBiExpoFit(xdata1,ydata1,f1,1,fname3)
        
        # Calculate the FWHM of the fitted Bi-Gaussian function
        fwhm = GetBiExpoFWM(fitcoeff,0.5)
        fwtm = GetBiExpoFWM(fitcoeff,0.1)
        
        
        # save the fit info, fwhm and fwtm, and R2 in an array
        print idx
        combodata[idx,:len(fitcoeff)] = fitcoeff
        combodata[idx,-3]=r2
        combodata[idx,-2]=fwhm
        combodata[idx,-1]=fwtm

    #save the combo curve fit results
    np.savetxt(fname1,combodata,fmt='%6.5e')

    # save into binary file
    fname5 = analysisdir + '/Vol_PostPosition_Run' + str(run1) + '_' + str(run2)
    SaveBinaryFile(vol1,fname5)
    
    #plt.show()
    #plt.close() 
    

if __name__ == '__main__':
    main()


