import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table
from spectral_cube import SpectralCube
import numpy as np
import subprocess

def generate_fakes(xoff=0.0, yoff=0.0, voff=0.0, datadir='../data/'):
    cube = SpectralCube.read(datadir+'m83.co10.K_correct.fits')
    dl = cube.header['CDELT2']*np.pi/180*4.6e6
    dv = cube.header['CDELT3']

    sub = cube[55:,:,0:600]

    mass = np.logspace(6.0,5.0,51)
    Lco = mass/4.35
    radius = (mass/np.pi/300)**0.5
    sigmar = (((radius/1.91)**2+(23.5**2))**0.5)/dl
    sigmav = ((0.7*(radius**0.5)+(2.6/2.355)**2)**(0.5))/dv
    amp = Lco/(sigmar**2*dl**2*dv*sigmav*(2*np.pi)**1.5)

    xloc = np.linspace(100,500,9)
    yloc = np.linspace(100,1500,29)
    vloc = np.linspace(5,45,5)

    v,y,x = np.meshgrid(vloc,yloc,xloc)

    sources = np.zeros((20,100,100))
    vcoord,ycoord,xcoord = np.meshgrid(np.arange(20),
                                       np.arange(100),
                                       np.arange(100),indexing = 'ij')
    noise = sub.filled_data[:]
    for thismass,thissigr,thissigv,thisamp in zip(mass,sigmar,sigmav,amp): 
        sources = np.exp(-((xcoord-49.5-xoff)**2/(2*thissigr**2)+
                           (ycoord-49.5-yoff)**2/(2*thissigr**2)+
                           (vcoord-10.0-voff)**2/(2*thissigv**2)))
        sources /= sources.sum() 
        sources *= thismass / 4.35 / dv / dl**2
        tt = np.tile(sources,(2,16,6))

        fakesources = np.copy(noise)
        fakesources[2:42,:,:]+=tt
        fakecube = SpectralCube(fakesources,sub.wcs,header=sub.header)
        fakecube.write(datadir+'m83.fakes_{:.2f}.fits'.format(np.log10(thismass)),overwrite=True)
    #    import pdb; pdb.set_trace()

def driveidl():
    command = 'idl -e "completeness"'
    subprocess.call(command,shell=True)
#call the IDL program here

def analyze_run(xoff=0.0, yoff=0.0, voff=0.0, datadir='../data/',
                outfile='completeness.fits'):
    logmass = np.linspace(6.0, 5.0, 51)
    yloc = np.linspace(50,1550,16)+xoff
    xloc = np.linspace(150,550,5)+yoff
    vloc = np.array([12.0,32.0])+voff

    xx,yy,vv = np.meshgrid(xloc,yloc,vloc,indexing='ij')

    xx = xx.ravel()
    yy = yy.ravel()
    vv = vv.ravel()
    idx = ~((xx==150)*(yy==1550))

    xx = xx[idx]
    yy = yy[idx]
    vv = vv[idx]
    threshold = 5
    origmass = np.array([])
    newmass = np.array([])
    compfrac = np.zeros_like(logmass)
    for idx2, lm in enumerate(logmass):
        xmatch = np.zeros(len(xx))-1
        t = Table.read('../measurements/m83.fakes_{0}_props_clfind.fits'.format(lm))
        for idx,(x, y, v) in enumerate(zip(xx, yy, vv)):
            dist = (np.abs(t['MOM1X']-x)+np.abs(t['MOM1Y']-y)+np.abs(t['MOM1V']-v))
            match = np.argmin(dist)
            if dist[match] < threshold:
                xmatch[idx] = match
        print lm, np.sum((xmatch > -1))/158.0
        compfrac[idx2] =  np.sum((xmatch > -1))/158.0
        xmatch = xmatch[xmatch!=-1]
        meanmass = np.median(t[xmatch.astype(np.bool)]['MASS_GCORR'].data)
#        newmass = np.r_[newmass,t[xmatch.astype(np.bool)]['MASS_GCORR'].data]
#        origmass = np.r_[origmass,lm*np.ones(len(xmatch))]
    results = Table()
    results['logmass'] = logmass
    results['FracRecover'] = compfrac
    results['MeanMass'] = meanmass
    results.write(outfile,overwrite=True)

def mcp(nrun=10):
    for i in np.arange(nrun):
        xoff = np.random.rand()*20-10
        yoff = np.random.rand()*20-10
        voff = np.random.rand()*5-2.5
        generate_fakes(xoff=xoff, yoff=yoff, voff=voff)
        driveidl()
        analyze_run(xoff=xoff, yoff=yoff, voff=voff,
                    outfile='comprun_{:05.3}_{:05.3}_{:05.3}.fits'.format(xoff,yoff,voff))
        
