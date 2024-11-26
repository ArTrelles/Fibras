
# convert gaia error data to fits file
from astropy.io import fits
#import pyfits
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from scipy import constants as const
#from FortranFile import FortranFile
import struct
#from galpy.util import bovy_coords
#from fortranfile import FortranFile

# input file name
inputfileint='GeneratedStarsle03.bin'
# output file
outfile='file.fits'

# Definir formato de registro: 33 doubles (d), 1 int (i), 7 doubles (d)
record_format = '<33d i 7d'
record_size = struct.calcsize(record_format)  # Tamaño del registro sin padding
print(f"Tamaño del registro (sin padding): {record_size} bytes")

# Padding: 4 bytes al inicio y al final
padded_record_size = record_size + 8  # Agregar padding inicial y final
print(f"Tamaño del registro con padding: {padded_record_size} bytes")

# Archivo de entrada
#input_file = "GeneratedStarsle03.bin"

# Lista para almacenar los datos
rdata = []

# Leer archivo binario con padding
# Convertir datos a un arreglo NumPy
rdata = np.array(rdata)
print(f"Datos cargados correctamente con forma: {rdata.shape}")
print(f"Primeros valores: {rdata[:5, :]}")


#print("Primeros valores de ao1:", rdata[:5, 0])  # ao1
#print("Primeros valores de ae1:", rdata[:5, 6])  # ae1
#print("Primeros valores de iden:", rdata[:5, 33])  # iden


# Verificar los datos
print(f"Datos cargados correctamente con forma: {rdata.shape}")


# Galactic coordinates
# True l and b
ao1=rdata[:,0]
ao2=rdata[:,1]
ao3=rdata[:,2]
ao4=rdata[:,3]
ao5=rdata[:,4] 
ao6=rdata[:,5] 

ae1=rdata[:,6]
ae2=rdata[:,7]
ae3=rdata[:,8]
ae4=rdata[:,9]
ae5=rdata[:,10] 
ae6=rdata[:,11] 

po1=rdata[:,12] 
po2=rdata[:,13] 
pe1=rdata[:,14] 
pe2=rdata[:,15] 
apo1=rdata[:,16] 
apo2=rdata[:,17] 
ape1=rdata[:,18] 
ape2=rdata[:,19]

xo_p= rdata[:,20]
yo_p= rdata[:,21]
zo_p= rdata[:,22]

vx_pn = rdata[:,23]
vy_pn = rdata[:,24]
vz_pn = rdata[:,25]

vxo_p = rdata[:,26]
vyo_p=rdata[:,27]
vzo_p=rdata[:,28]

mgen = rdata[:,29]
z = rdata[:,30]
PAge = rdata[:,31]

priot = rdata[:,32]
iden = rdata[:,33]

r7 = rdata[:,34]
gr7=rdata[:,35]
Z7 =rdata[:,36]

lo_p=rdata[:,37]
bo_p=rdata[:,38]
r_n=rdata[:,39]
gr_n=rdata[:,40]



tbhdu = fits.BinTableHDU.from_columns([\
  fits.Column(name='ao1',unit='(radians)',format='D',array=rdata[:,0]),\
  fits.Column(name='ao2',unit='(radians)',format='D',array=rdata[:,1]),\
  fits.Column(name='ao3',unit='(as)',format='D',array=rdata[:,2]),\
  fits.Column(name='ao4',unit='(as/yr)',format='D',array=rdata[:,3]),\
  fits.Column(name='ao5',unit='(as/yr)',format='D',array=rdata[:,4]),\
  fits.Column(name='ao6',unit='(km/s)',format='D',array=rdata[:,5]),\
  fits.Column(name='ae1',unit='(radians)',format='D',array=rdata[:,6]),\
  fits.Column(name='ae2',unit='(radians)',format='D',array=rdata[:,7]),\
  fits.Column(name='ae3',unit='(as)',format='D',array=rdata[:,8]),\
  fits.Column(name='ae4',unit='(as/yr)',format='D',array=rdata[:,9]),\
  fits.Column(name='ae5',unit='(as/yr)',format='D',array=rdata[:,10]),\
  fits.Column(name='ae6',unit='(km/s)',format='D',array=rdata[:,11]),\
  fits.Column(name='po1',unit='(magnitudes)',format='D',array=rdata[:,12]),\
  fits.Column(name='po2',unit='(magnitudes)',format='D',array=rdata[:,13]),\
  fits.Column(name='pe1',unit='(magnitudes)',format='D',array=rdata[:,14]),\
  fits.Column(name='pe2',unit='(magnitudes)',format='D',array=rdata[:,15]),\
  fits.Column(name='apo1',unit='(kelvin)',format='D',array=rdata[:,16]),\
  fits.Column(name='apo2',unit='(dex)',format='D',array=rdata[:,17]),\
  fits.Column(name='ape1',unit='(kelvin)',format='D',array=rdata[:,18]),\
  fits.Column(name='ape2',unit='(dex)',format='D',array=rdata[:,19]),\
  fits.Column(name='xo_p',unit='(kpc)',format='D',array=rdata[:,20]),\
  fits.Column(name='yo_p',unit='(kpc)',format='D',array=rdata[:,21]),\
  fits.Column(name='zo_p',unit='(kpc)',format='D',array=rdata[:,22]),\
  fits.Column(name='vx_pn',unit='(km/s)',format='D',array=rdata[:,23]),\
  fits.Column(name='vy_pn',unit='(km/s)',format='D',array=rdata[:,24]),\
  fits.Column(name='vz_pn',unit='(km/s)',format='D',array=rdata[:,25]),\
  fits.Column(name='vxo_pn',unit='(mag)',format='D',array=rdata[:,26]),\
  fits.Column(name='vyo_pn',unit='(km/s)',format='D',array=rdata[:,27]),\
  fits.Column(name='vzo_pn',unit='(km/s)',format='D',array=rdata[:,28]),\
  fits.Column(name='mgen',unit='(msun)',format='D',array=rdata[:,29]),\
  fits.Column(name='z',unit='(dex)',format='D',array=rdata[:,30]),\
  fits.Column(name='PAge',unit='(Gyr)',format='D',array=rdata[:,31]),\
  fits.Column(name='priot',format='D',array=rdata[:,32]),\
  fits.Column(name='iden',format='K',array=rdata[:,33]),\
  fits.Column(name='r7',unit='(mag)',format='D',array=rdata[:,34]),\
  fits.Column(name='gr7',unit='(mag)',format='D',array=rdata[:,35]),\
  fits.Column(name='z7',unit='(mag)',format='D',array=rdata[:,36]),\
  fits.Column(name='lo_p',unit='(radian)',format='D',array=rdata[:,37]),\
  fits.Column(name='bo_p',unit='(radian)',format='D',array=rdata[:,38]),\
  fits.Column(name='r_n',unit='(mag)',format='D',array=rdata[:,39]),\
  fits.Column(name='gr_n',unit='(dex)',format='D',array=rdata[:,40]) ])
tbhdu.writeto(outfile,overwrite=True)
