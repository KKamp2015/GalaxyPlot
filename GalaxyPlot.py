import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import imread
from astropy.io import ascii


def arms():

	narms=5
	arms=ascii.read("arms.csv",format='csv')#arms.inp rand through excel to make it a csv
	aarm=arms['a'][:]
	rmin=arms['rmin'][:]
	thmin=arms['thmin'][:]
	extent=arms['Extent'][:]

	nth=1000
	r1=np.zeros((nth,narms))
	th1=np.zeros((nth,narms))
	rad=180/np.pi

	for j in range(0,narms):
		th1[:,j]=np.arange(0,nth,1)/(nth-1)*extent[j]+thmin[j]
		th1[:,j]=np.arange(0,nth,1)/(nth-1)*extent[j]+thmin[j]
		r1[:,j]=rmin[j]*np.exp((th1[:,j]-thmin[j])/aarm[j])
		th1[:,j]=th1[:,j]*180/np.pi
		if j==1:
			g=np.where((th1[:,j]>370.0) & (th1[:,j]<= 410.0))
			r1[g,j]=r1[g,j]*(1. + 0.04*np.cos((th1[g,j]-390.)*180./(40.*rad)))
			g=np.where((th1[:,j]>315) &( th1[:,j]<=370))
			r1[g,j]=r1[g,j]*(1. - 0.07*np.cos((th1[g,j]-345.)*180./(55.*rad)))
			g=np.where((th1[:,j]>180)&(th1[:,j]<+315))
			r1[g,j]=r1[g,j]*(1 + 0.16*np.cos((th1[g,j]-260.)*180./(135.*rad)))
		if j==3:
			g=np.where((th1[:,j]>290)&(th1[:,1]<+395))
			r1[g,j]=r1[g,j]*(1. - 0.11*np.cos((th1[g,j]-350.)*180./(105.*rad)))
	yy=r1*np.cos(np.radians(th1))
	xx=-r1*np.sin(np.radians(th1))
	return [th1,r1,xx,yy]

def plotarms(r1,th1,xx,yy):
	r1=r1*8.0/8.5
	xx=xx*8.0/8.5
	yy=yy*8.0/8.5
	plt.scatter(xx,yy,c='white')
	plt.xlim([-16,16])
	plt.ylim([-16,16])

	dr=0.375
	for i in range(0,4):
		rout=r1[:][i]+dr
		rin=r1[:][i]-dr
		yyout=rout*np.cos(np.radians(th1[:][i]))
		xxout=-rout*np.sin(np.radians(th1[:][i]))
		yyin=rin*np.cos(np.radians(th1[:][i]))
		xxin=-rin*np.sin(np.radians(th1[:][i]))
		xxa=[xxin,xxout[::-1],xxin[0]]
		yya=[yyin,yyout[::-1],yyin[0]]
		#pch.Patch(xxa,yya,c='b')
	plt.scatter(0,8,marker="*",c="yellow",zorder=100,label="Sun") #this is the sun
	plt.text(-4.2,8.5,"Local",color='White')
	plt.text(-8.99662,-14.3011,'SAGITTARIUS-CARINA',color='White')
	plt.text(1.02473,-12.1064,'SCUTUM-CRUX',color='White')
	plt.text(7.0,11.3,'NORMA-CYGNUS',color='White')
	plt.text(-10.4,10.3,'PERSEUS',color='White')

	rbar=3.4
	ibar=20
	xb=np.sin(np.radians(ibar))*rbar
	yb=np.cos(np.radians(ibar))*rbar
	plt.plot([-xb,xb],[-yb,yb],linewidth=5,c='white')
	plt.scatter(0,0,marker='*',s=10,zorder=1000)


#Original file has to be read in excel to make it csv
infile=ascii.read("CoordBe.csv",data_start=2)
r=infile['avg'][:674]
print(r[673])
l=infile['l'][:674]
b=infile['b'][:674]
ax=plt.gca()
d=r*np.cos(np.radians(r))

x=d*np.sin(np.radians(l))/1000

y=(-d*np.cos(np.radians(d))+8000)/1000
'''
fig=plt.figure()
fig.set_size_inches(11, 8.5)
'''
plt.figure(figsize=(15,15))
plt.xlabel("x [kpc]")
plt.ylabel("y [kpc]")
th1,r1,xx,yy=arms()
plotarms(r1,th1,xx,yy)

#making image Backround
img=imread('MilkyWay.jpg')
p=20.5
plt.imshow(img,zorder=0,extent=[-p,p,-p,p])


plt.scatter(x[:],y[:],marker='o',c='C0',s=2,label="Be Stars")
thcirc=np.linspace(0,360,360)
rcirc=0.5
xcirc=rcirc*np.cos(np.radians(thcirc))
ycirc=rcirc*np.sin(np.radians(thcirc))+8
plt.plot(xcirc,ycirc,zorder=99,linewidth=1,c='r')
#plt.text(0,8.51,'500 pc',color='r',horizontalalignment='center')
#plt.legend(loc=4)
plt.title('Be Star Distribution over Galactic Spiral Arms')
plt.savefig("TopDownKeefeMilky2.jpg",format='jpg')
plt.show()


