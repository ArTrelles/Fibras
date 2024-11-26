from astropy.io import fits

# Archivo FITS
fits_file = 'file.fits'

# Abrir el archivo FITS
hdul = fits.open(fits_file)

# Ver información general del archivo
print(hdul.info())


# Ver encabezado (header) de la tabla
header = hdul[1].header
print(repr(header))

# Obtener los datos de la tabla (BinTableHDU)
data = hdul[1].data  # Índice 1 porque la tabla está en la segunda extensión

# Mostrar las primeras filas
print(data[:5])


# Acceder a una columna por nombre
ao1 = data['ao1']
print(ao1[:10])  # Mostrar los primeros 10 valores de 'ao1'

# Otra columna
priot = data['priot']
print(priot[:10])  # Mostrar los primeros 10 valores de 'priot'

hdul.close()
