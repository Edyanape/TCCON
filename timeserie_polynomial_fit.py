
## Imported modules ##

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import datetime

## Definition method inverse funtion ##

def subset(dset,limits):
	for k,limite in enumerate(limits):
		condition=(np.logical_and((min([dset[limits[k][0]]]) > limits[k][1]),(max([dset[limits[k][0]]]) < limits[k][2])))
		if k==0:
			conditions=condition

		else:
			conditions=(np.logical_and(conditions,condition))

	indices=np.where(conditions)[0]
	subdset=dset[indices]
	return subdset

## Definition of a inverse function ##

def fitcicle(tarr,y,ncoef,npoly):
    tcenter=np.average(tarr)
    tarrfull=24*3600*np.arange((int(np.max(tarr)-np.min(tarr)))/(24*3600))+np.min(tarr)
    tepoch=tarr-tcenter
    tepochfull=tarrfull-tcenter
    ntimes=len(tepoch)
    ntimesfull=len(tepochfull)
    kmat=np.zeros((ntimes,npoly+2*ncoef),dtype=float)
    kmatfull=np.zeros((ntimesfull,npoly+2*ncoef),dtype=float)
    omega=2*np.pi/(365*24*3600)
    for i,ttt in enumerate(tepoch):
        for k in range(npoly):
            kmat[i,k]=ttt**k
        for j in range(ncoef):
            kmat[i,2*j+npoly]=np.cos((j+1)*omega*ttt) 
            kmat[i,2*j+npoly+1]=np.sin((j+1)*omega*ttt) 
    for i,ttt in enumerate(tepochfull):
        for k in range(npoly):
            kmatfull[i,k]=ttt**k
        for j in range(ncoef):
            kmatfull[i,2*j+npoly]=np.cos((j+1)*omega*ttt) 
            kmatfull[i,2*j+npoly+1]=np.sin((j+1)*omega*ttt) 

    gainmatrix=KTKinvKT=np.dot(np.linalg.inv(np.dot(kmat.T,kmat)),kmat.T)
    x=np.dot(KTKinvKT,y)
    yfit=np.dot(kmatfull,x)
    yy=np.dot(kmat,x)
    rms=np.std(y-yy)
    Serr=rms**2*np.dot(gainmatrix,gainmatrix.T)
    errores=np.sqrt(np.diag(Serr))
    t0=datetime.datetime.fromtimestamp(0.0)
    tfit=[t0+datetime.timedelta(seconds=tepoch) for tepoch in tarrfull]
    return x,yy,tfit,yfit,errores

## Polynomial fit function ##
def fitcicleold(tarr,y,ncoef,npoly):
    tcenter=np.average(tarr)
    tarrfull=24*3600*np.arange((int(np.max(tarr)-np.min(tarr)))/(24*3600))+np.min(tarr)
    tepoch=tarr-tcenter
    tepochfull=tarrfull-tcenter
    ntimes=len(tepoch)
    ntimesfull=len(tepochfull)
    kmat=np.zeros((ntimes,npoly+2*ncoef),dtype=float)
    kmatfull=np.zeros((ntimesfull,npoly+2*ncoef),dtype=float)
    omega=2*np.pi/(365*24*3600)
    for i,ttt in enumerate(tepoch):
        for k in range(npoly):
            kmat[i,k]=ttt**k
        for j in range(ncoef):
            kmat[i,2*j+npoly]=np.cos((j+1)*omega*ttt) 
            kmat[i,2*j+npoly+1]=np.sin((j+1)*omega*ttt) 
    for i,ttt in enumerate(tepochfull):
        for k in range(npoly):
            kmatfull[i,k]=ttt**k
        for j in range(ncoef):
            kmatfull[i,2*j+npoly]=np.cos((j+1)*omega*ttt) 
            kmatfull[i,2*j+npoly+1]=np.sin((j+1)*omega*ttt) 

    gainmatrix=KTKinvKT=np.dot(np.linalg.inv(np.dot(kmat.T,kmat)),kmat.T)
    x=np.dot(KTKinvKT,y)
    yfit=np.dot(kmatfull,x)
    yy=np.dot(kmat,x)
    t0=datetime.datetime.fromtimestamp(0.0)
    tfit=[t0+datetime.timedelta(seconds=tepoch) for tepoch in tarrfull]
    return x,yy,tfit,yfit

## Extraction of data from spectrum files ##

## Document with spectrum data ##
filename = 'datos.txt'

## Import and transformation ##
dataset = pd.read_csv(filename,delim_whitespace=True)
df = pd.DataFrame(dataset)
df['Columv']=df['OVC_co2']*df['VSF_co2']
datetimes = np.array(pd.to_datetime(df['time']))
ndata=len(datetimes)

cols=np.array(df['Columv'])
print (ndata)

colavg=np.average(cols)
colstd=np.std(cols)

print (colavg,colstd)

## Identification of outliers ##

plt.plot(datetimes,cols,'.b')

index=np.where(np.abs(cols-colavg) < 3*colstd)[0]

plt.plot(datetimes[index],cols[index],'.r')
plt.title('ALTZOMONI-TCCON')
plt.ylabel(u'CO\u2082' + ' ' + '[Moléculas/cm\u00b2]' )

plt.show()

## Time serie without outliers and polynomial fit ##
datetimes=datetimes[index]
cols=cols[index]
plt.plot(datetimes,cols,'.b')
to=datetimes[0]
print (datetimes[0])
print (to)
tarr=np.array([(ttt-to)/np.timedelta64(1,'s') for ttt in datetimes])
x,yy,tfit,yfit,errores=fitcicle(tarr,cols,4,2)
print( x)
print (x[1]/x[0]*100*3600*24*365 ,'%/yr')
print (400 *x[1]/x[0]*3600*24*365 ,'ppm/yr')

plt.plot(datetimes,yy,'r-')
plt.title('ALTZOMONI-TCCON')
plt.ylabel(u'CO\u2082' + ' ' + '[Moléculas/cm\u00b2]' )
plt.show()

## Histogram of range of ppm in CO2 ##
dif=cols-yy
std=np.std(dif)
plt.hist(dif,100)
plt.show()

index=np.where(np.abs(dif) < 1.0E20)[0]


datetimes=datetimes[index]
cols=cols[index]
plt.plot(datetimes,cols,'.b')
#to=datetimes[0]
to=np.datetime64('1970-01-01T00:00:00Z')
print (datetimes[0])
print (to)
tarr=np.array([(ttt-to)/np.timedelta64(1,'s') for ttt in datetimes])
x,yy,tfit,yfit,errores=fitcicle(tarr,cols,2,2)
print( x)
print (errores)
print (errores[1]/x[1]*100.0 , '%')
print (x[1]/x[0]*100*3600*24*365 ,'%/yr')
print (errores[1]/x[0]*100*3600*24*365 ,'\pm %/yr')
print (400 *x[1]/x[0]*3600*24*365 ,'ppm/yr')
print (400 *errores[1]/x[0]*3600*24*365 ,'\pm ppm/yr')
promedio = 400 *x[1]/x[0]*3600*24*365
err = 400 *errores[1]/x[0]*3600*24*365 
plt.plot(datetimes,yy,'k-')
plt.show()

## Final time serie with smoothed fit ##

plt.plot(datetimes,cols,'.b')
plt.plot(tfit,yfit,'r-')
#plt.legend(['Columnas de'+ ' ' + u'O\u2082', 'Ajuste suavizado'  ])
plt.title("Tendencia ="+' '+ str(round(promedio,3)) + ' ± '+ str(round(err,3))+' ppm/año ≈' + ' '+str(round(promedio/400*100,3))+ ' ± '+str(round(err/400*100,3)),
         position=(0.5, 0.9),
         fontdict={'size': 10})
#plt.title('ALTZOMONI-TCCON')
plt.ylabel(u'CO\u2082' + ' ' + '[Moléculas/cm\u00b2]' )
plt.show()
