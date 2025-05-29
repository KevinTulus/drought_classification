import os
import shutil
from glob import glob
from datetime import datetime
import calendar
import ee
import geemap
import rasterio
from rasterio.mask import mask
from rasterio.warp import calculate_default_transform, reproject, Resampling
import fiona
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')  # Atur backend ke Agg sebelum memuat pyplot
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from sklearn.metrics import confusion_matrix, accuracy_score, classification_report
import klasifikasi

# Inisialisasi Google Earth Engine
ee.Initialize(project='spi-landsat')

# Function
def delete_files_in_directory(path):
    try:
        files = glob(os.path.join(path, '*'))
        for file in files:
            if os.path.isfile(file):
                os.remove(file)
        return("All files deleted successfully.")
    except OSError:
        return("Error occurred while deleting files.")

def crop_raster_with_shapefile(raster_path, shapefile_path, output_path):
    """
    Memotong raster menggunakan shapefile dan menyimpan hasilnya.
    """
    with fiona.open(shapefile_path, "r") as shapefile:
        # Membaca geometri dari shapefile
        shapes = [feature["geometry"] for feature in shapefile]

    with rasterio.open(raster_path) as src:
        # Crop raster dengan geometri shapefile
        out_image, out_transform = mask(src, shapes, crop=True)
        out_meta = src.meta.copy()

        # Update metadata untuk raster hasil crop
        out_meta.update({
            "driver": "GTiff",
            "height": out_image.shape[1],
            "width": out_image.shape[2],
            "transform": out_transform
        })

        # Simpan raster hasil crop
        with rasterio.open(output_path, "w", **out_meta) as dest:
            dest.write(out_image)

def replace_nodata(input_path, output_path, new_nodata_value):
    """
    Mengganti nilai NoData raster menjadi nilai baru.
    """
    with rasterio.open(input_path) as src:
        data = src.read(1)  # Membaca band pertama
        meta = src.meta

        # Ganti nilai NoData dari 0 ke -9999
        data[data == 0] = new_nodata_value

        # Perbarui metadata untuk NoData
        meta.update(nodata=new_nodata_value)

        # Simpan raster hasil perubahan
        with rasterio.open(output_path, "w", **meta) as dst:
            dst.write(data, 1)

def resample_raster(src_path, dst_path, ref_transform, ref_crs, ref_width, ref_height):
    with rasterio.open(src_path) as src:
        transform, width, height = calculate_default_transform(
            src.crs, ref_crs, ref_width, ref_height, *src.bounds
        )
        kwargs = src.meta.copy()
        kwargs.update({
            "transform": ref_transform,
            "width": ref_width,
            "height": ref_height,
            "crs": ref_crs
        })

        with rasterio.open(dst_path, "w", **kwargs) as dst:
            reproject(
                source=rasterio.band(src, 1),
                destination=rasterio.band(dst, 1),
                src_transform=src.transform,
                src_crs=src.crs,
                dst_transform=ref_transform,
                dst_crs=ref_crs,
                resampling=Resampling.bilinear  # Pilih metode resampling
            )

def copy_file(path_file, path_folder):
    # Salin file
    shutil.copy2(path_file, path_folder)  # copy2 mempertahankan metadata asli

def apply_nodata_to_raster(raster_path, dem=False):
    reference_raster_path = "static/output/VHI.tiff"
    # Membuka raster input
    with rasterio.open(raster_path) as input_raster:
        input_data = input_raster.read()  # Membaca semua data raster
        input_meta = input_raster.meta  # Metadata raster input

    # Membuka raster referensi
    with rasterio.open(reference_raster_path) as reference_raster:
        reference_data = reference_raster.read(1)  # Membaca band pertama sebagai referensi
        nodata_value = reference_raster.nodata  # Mendapatkan nilai NoData dari referensi

    if dem:
        nodata_value = -32767  # Default nilai NoData DEM

    # Membuat mask berdasarkan nilai NoData pada referensi
    nodata_mask = reference_data == -9999

    # Terapkan mask NoData ke raster input
    for band in range(input_data.shape[0]):
        input_data[band][nodata_mask] = nodata_value

    # Update metadata untuk menyimpan raster output dengan NoData
    input_meta.update(nodata=nodata_value, dtype='float32')

    # Menyimpan hasil raster dengan NoData diterapkan
    with rasterio.open(raster_path, 'w', **input_meta) as output_raster:
        output_raster.write(input_data.astype('float32'))

    return(f"NoData dari referensi telah diterapkan dan hasil disimpan di {raster_path}")

# Download process
def download_raster(year, month, wrs_path, wrs_row):
    """
    Mengunduh citra Landsat 8 Collection 2 Level-2 untuk satu bulan tertentu.
    """
    output_dir = "static/download" # Direktori output
    bands = ['SR_B4', 'SR_B5', 'ST_B10'] # Daftar band yang akan diunduh
    # bands_qa = ['SR_B4', 'SR_B5', 'ST_B10', 'QA_PIXEL'] # Daftar band yang akan diunduh

    # Tentukan bulan awal dan akhir
    start_date = datetime(year, month, 1)
    end_date = datetime(year, month, calendar.monthrange(year, month)[1])

    # Filter koleksi Landsat berdasarkan Path/Row yang ditentukan
    path_row = ee.Filter.And(
        ee.Filter.eq('WRS_PATH', wrs_path),
        ee.Filter.eq('WRS_ROW', wrs_row)
    )
    
    # Filter koleksi citra
    collection = (ee.ImageCollection("LANDSAT/LC08/C02/T1_L2")
                  .filter(path_row)
                  .filterDate(start_date.strftime('%Y-%m-%d'), end_date.strftime('%Y-%m-%d'))
                  .filter(ee.Filter.lt('CLOUD_COVER', 20))
                  .select(bands))
    
    # # Menambahkan cloud masking menggunakan QA band
    # def apply_cloud_mask(image):
    #     # Mendapatkan band 'pixel_qa' dari citra
    #     qa60 = image.select('QA_PIXEL')
        
    #     # Membuat mask untuk awan (bit 3 == 1 berarti ada awan)
    #     cloud_mask = qa60.bitwiseAnd(1 << 3).eq(0)  # 0 berarti tidak ada awan
        
    #     # Terapkan mask ke citra
    #     return image.updateMask(cloud_mask)
    
    # # Terapkan cloud mask pada semua citra dalam koleksi
    # collection = collection.map(apply_cloud_mask)
    
    # Periksa apakah koleksi memiliki citra
    if collection.size().getInfo() > 0:
        # Gabungkan koleksi menjadi median
        image = collection.median()
        year_month = f"{year}_{month:02d}"
        
        # Unduh tiap band secara terpisah
        for band in bands:
            band_image = image.select(band)
            band_name = band
            filename = os.path.join(output_dir, f"landsat_{band_name}.tif")
            
            # Ekspor citra menggunakan geemap
            geemap.ee_export_image(
                ee_object=band_image,
                filename=filename,
                scale=1000,
                region=collection.geometry(),
                crs="EPSG:4326"
            )
        
        # return False, "Berhasil mengunduh citra"
        return False
    else:
        # return True, f"Tidak ada citra memenuhi syarat untuk bulan {month}-{year}."
        return True

def crop_raster():
    shapefile_path = 'static/database/Jatim_shp/RBI_50K_2023_Jawa Timur.shp'
    input_folder = glob('static/download/*.tif')
    output_folder = 'static/output/Crop'
    for raster_file in input_folder:
        output_filename = os.path.basename(raster_file).replace(".tif", ".tiff")
        output_path = os.path.join(output_folder, output_filename)
        crop_raster_with_shapefile(raster_file, shapefile_path, output_path)

def change_nodata():
    # Proses semua raster di folder input
    new_nodata = -9999
    raster_files = glob('static/output/Crop/*.tiff')
    output_folder = f"static/output/Nodata"

    for raster_file in raster_files:
        # Nama file output
        output_filename = os.path.basename(raster_file)
        output_path = os.path.join(output_folder, output_filename)

        # Ganti NoData
        replace_nodata(raster_file, output_path, new_nodata)

def resample_download(wrs_path, wrs_row):
    # Baca metadata raster referensi
    with rasterio.open(glob(f'static/database/NDVI/{wrs_path}_{wrs_row}/*.tiff')[0]) as ref:
        ref_transform = ref.transform
        ref_crs = ref.crs
        ref_width = ref.width
        ref_height = ref.height
        ref_bounds = ref.bounds

    raster_files = glob('static/output/Nodata/*.tiff')
    output_folder = f"static/output/Resample"

    for raster_file in raster_files:
        # Nama file output
        output_filename = os.path.basename(raster_file)
        output_path = os.path.join(output_folder, output_filename)

        # Ganti NoData
        resample_raster(raster_file, output_path, ref_transform, ref_crs, ref_width, ref_height)

    return("DEM berhasil dibuat.\n")

def start_download(year, month, wrs_path, wrs_row):
    error = False
    error = download_raster(int(year), int(month), int(wrs_path), int(wrs_row))
    crop_raster()
    change_nodata()
    resample_download(wrs_path, wrs_row)
    return error

# Check and get metadata
def check_rasters_alignment():
    file = glob('static/output/Resample/*.tiff')
    error = []

    with rasterio.open(file[1]) as src1, rasterio.open(file[2]) as src2, rasterio.open(file[0]) as src3:
        # Memeriksa dimensi (size) ketiga raster
        size_match = (src1.width == src2.width == src3.width) and (src1.height == src2.height == src3.height)
        
        # Memeriksa transformasi spasial (affine transform) ketiga raster
        transform_match = (src1.transform == src2.transform == src3.transform)
        
        # Memeriksa CRS (Coordinate Reference System)
        crs_match = (src1.crs == src2.crs == src3.crs)

        # Menampilkan hasil pemeriksaan
        if size_match and transform_match and crs_match:
            return False, "Validasi berhasil. Semua raster memiliki ukuran, posisi, dan CRS yang sama.\n"
        else:
            if not size_match:
                error.append(f"Ukuran berbeda: Raster1 ({src1.width}x{src1.height}), Raster2 ({src2.width}x{src2.height}), Raster3 ({src3.width}x{src3.height})")
                # return{"error": f"Ukuran berbeda: Raster1 ({src1.width}x{src1.height}), Raster2 ({src2.width}x{src2.height}), Raster3 ({src3.width}x{src3.height})"}
            if not transform_match:
                error.append(f"Transformasi berbeda: Raster1 = {src1.transform}, Raster2 = {src2.transform}, Raster3 = {src3.transform}")
                # return{"error": f"Transformasi berbeda: Raster1 {src1.transform}, Raster2 {src2.transform}, Raster3 {src3.transform}"}
            if not crs_match:
                error.append(f"CRS berbeda: Raster1 = {src1.crs}, Raster2 = {src2.crs}, Raster3 = {src3.crs}")
                # return{"error": f"CRS berbeda: Raster1 {src1.crs}, Raster2 {src2.crs}, Raster3 {src3.crs}"}
            return True, error

def get_raster_metadata():
    input = glob('static/output/Resample/*.tiff')[0]
    with rasterio.open(input) as src:
        metadata = {
            "width": src.width,
            "height": src.height,
            "crs": src.crs.to_string() if src.crs else None,
            # "transform": src.transform,
            # "bounds": src.bounds,
            # "driver": src.driver,
            # "count": src.count,
            "dtypes": src.dtypes[0],
            "nodata": src.nodata,
            # "indexes": src.indexes[0],
            # "res": src.res,
        }
    return metadata

# Preprocessing and classification Process
def calculate_ndvi():
    with rasterio.open(glob('static/output/Resample/*B4.tiff')[0]) as red_src:
        red = red_src.read(1, masked=True).astype('float32') * 0.0000275 - 0.2
        red_meta = red_src.meta
    
    with rasterio.open(glob('static/output/Resample/*B5.tiff')[0]) as nir_src:
        nir = nir_src.read(1, masked=True).astype('float32') * 0.0000275 - 0.2
    
    numerator = nir - red
    denominator = nir + red
    ndvi = np.where(denominator == 0, 0, numerator / denominator)

    # ndvi.fill_value = -9999
    output_path = f'static/output/NDVI.tiff'
    red_meta.update(dtype='float32', count=1, nodata=-9999)

    with rasterio.open(output_path, 'w', **red_meta) as dst:
        dst.write(ndvi, 1)
    return("NDVI berhasil dibuat.\n")

def calculate_lst():
    with rasterio.open(glob('static/output/Resample/*B10.tiff')[0]) as tirs_src:
        tirs = tirs_src.read(1, masked=True).astype('float32')
        tirs_meta = tirs_src.meta

    lst = (tirs * 0.00341802 + 149) - 273.15

    # tirs.fill_value = -9999
    output_path = f'static/output/LST.tiff'
    tirs_meta.update(dtype='float32', count=1, nodata=-9999)
    
    with rasterio.open(output_path, 'w', **tirs_meta) as dst:
        dst.write(lst, 1)
    return("LST berhasil dibuat.\n")

def calculate_savi():
    with rasterio.open(glob('static/output/Resample/*B4.tiff')[0]) as red_src:
        red = red_src.read(1, masked=True).astype('float32') * 0.0000275 - 0.2
        red_meta = red_src.meta
    
    with rasterio.open(glob('static/output/Resample/*B5.tiff')[0]) as nir_src:
        nir = nir_src.read(1, masked=True).astype('float32') * 0.0000275 - 0.2
    
    savi = ((nir - red) / (nir + red + 0.5)) * (1.5)

    # savi.fill_value = -9999
    output_path = "static/output/SAVI.tiff"
    red_meta.update(dtype='float32', count=1, nodata=-9999)
    
    with rasterio.open(output_path, 'w', **red_meta) as dst:
        dst.write(savi, 1)
    return("SAVI berhasil dibuat.\n")

def calculate_vci(path, row):
    inputs = glob(f'static/database/NDVI/{path}_{row}/*.tiff')

    # ndvi_stack = []
    # for input in inputs:
    #     with rasterio.open(input) as src:
    #         ndvi_stack.append(src.read(1, masked=True).astype("float32"))

    # ndvi_stack = np.ma.stack(ndvi_stack)
    # ndvi_min = ndvi_stack.min(axis=0)
    # ndvi_max = ndvi_stack.max(axis=0)

    ndvi_min = float('inf')
    ndvi_max = float('-inf')

    for raster_file in inputs:
        with rasterio.open(raster_file) as src:
            data = src.read(1, masked=True).astype("float32")
            ndvi_min = min(ndvi_min, data.min())
            ndvi_max = max(ndvi_max, data.max())

    with rasterio.open(f'static/output/NDVI.tiff') as ndvi_src:
        ndvi = ndvi_src.read(1, masked=True).astype('float32')
        ndvi_meta = ndvi_src.meta
    
    vci = ((ndvi - ndvi_min) / (ndvi_max - ndvi_min)) * 100

    # vci.fill_value = -9999
    output_path = "static/output/VCI.tiff"
    
    with rasterio.open(output_path, 'w', **ndvi_meta) as dst:
        dst.write(vci, 1)
    return("VCI berhasil dibuat.\n")

def calculate_tci(path, row):
    inputs = glob(f'static/database/LST/{path}_{row}/*.tiff')

    # lst_stack = []
    # for input in inputs:
    #     with rasterio.open(input) as src:
    #         lst_stack.append(src.read(1, masked=True).astype("float32"))

    # lst_stack = np.ma.stack(lst_stack)
    # lst_min = lst_stack.min(axis=0)
    # lst_max = lst_stack.max(axis=0)

    lst_min = float('inf')
    lst_max = float('-inf')

    for raster_file in inputs:
        with rasterio.open(raster_file) as src:
            data = src.read(1, masked=True).astype("float32")
            lst_min = min(lst_min, data.min())
            lst_max = max(lst_max, data.max())

    with rasterio.open(f'static/output/LST.tiff') as lst_src:
        lst = lst_src.read(1, masked=True).astype('float32')
        lst_meta = lst_src.meta
    
    tci = ((lst_max - lst) / (lst_max - lst_min)) * 100

    # tci.fill_value = -9999
    output_path = "static/output/TCI.tiff"
    
    with rasterio.open(output_path, 'w', **lst_meta) as dst:
        dst.write(tci, 1)
    return("TCI berhasil dibuat.\n")

def calculate_vhi():
    with rasterio.open('static/output/VCI.tiff') as vci_src:
        vci = vci_src.read(1, masked=True).astype('float32')
        vci_meta = vci_src.meta

    with rasterio.open('static/output/TCI.tiff') as tci_src:
        tci = tci_src.read(1, masked=True).astype('float32')
    
    vhi = 0.5 * vci + 0.5 * tci

    # vhi.fill_value = -9999
    output_path = "static/output/VHI.tiff"
    
    with rasterio.open(output_path, 'w', **vci_meta) as dst:
        dst.write(vhi, 1)
    return("VHI berhasil dibuat.\n")

def calculate_spi(year, month):
    # Baca metadata raster referensi
    with rasterio.open('static/output/VHI.tiff') as ref:
        ref_transform = ref.transform
        ref_crs = ref.crs
        ref_width = ref.width
        ref_height = ref.height
        ref_bounds = ref.bounds

    raster_file = f'static/database/SPI/spi_{year}{month}01.tif'
    output_path = "static/output/SPI.tiff"

    resample_raster(raster_file, output_path, ref_transform, ref_crs, ref_width, ref_height)
    apply_nodata_to_raster(output_path)

    return("SPI berhasil dibuat.\n")

def calculate_dem():
    # Baca metadata raster referensi
    with rasterio.open('static/output/VHI.tiff') as ref:
        ref_transform = ref.transform
        ref_crs = ref.crs
        ref_width = ref.width
        ref_height = ref.height
        ref_bounds = ref.bounds

    raster_file = 'static/database/DEM.tiff'
    output_path = "static/output/DEM.tiff"

    resample_raster(raster_file, output_path, ref_transform, ref_crs, ref_width, ref_height)
    apply_nodata_to_raster(output_path, True)
    return("DEM berhasil dibuat.\n")

def calculate_nst():
    with rasterio.open('static/output/LST.tiff') as lst_src:
        lst = lst_src.read(1, masked=True).astype('float32')
        lst_min = lst.min() # dalam °C
        lst_max = lst.max() # dalam °C
        lst_meta = lst_src.meta
    
    nst = (lst - lst_min) / (lst_max - lst_min)

    # nst.fill_value = -9999
    output_path = "static/output/NST.tiff"
    
    with rasterio.open(output_path, 'w', **lst_meta) as dst:
        dst.write(nst, 1)
    return("NST berhasil dibuat.\n")

def classify_drought(pdsi):
    """Fungsi klasifikasi nilai PDSI menjadi level kekeringan."""
    if pdsi == -9999:  # Jika nilai adalah NaN, tetap NaN
        return np.nan
    elif pdsi <= -4:
        return 0  # Kekeringan Ekstrem
    elif -4 < pdsi <= -3:
        return 1  # Kekeringan Parah
    elif -3 < pdsi <= -2:
        return 2  # Kekeringan Sedang
    elif -2 < pdsi <= -1:
        return 3  # Kekeringan Ringan
    else:
        return 4  # Tidak ada kekeringan (Normal)

def calculate_pdsi(year, month):
    # Baca metadata raster referensi
    with rasterio.open('static/output/VHI.tiff') as ref:
        ref_transform = ref.transform
        ref_crs = ref.crs
        ref_width = ref.width
        ref_height = ref.height
        ref_bounds = ref.bounds

    raster_file = f'static/database/PDSI_raw/scpdsi_{year}-{month}-01.tif'
    output_path = "static/output/PDSI.tiff"
    resample_raster(raster_file, output_path, ref_transform, ref_crs, ref_width, ref_height)
    apply_nodata_to_raster(output_path)

    with rasterio.open('static/output/PDSI.tiff') as src:
        meta = src.meta.copy()  # Salin metadata
        data = src.read(1).astype(np.float32)  # Baca band pertama

        # Masking nilai NaN
        # masked_data = np.ma.masked_invalid(data)  # Masking NaN tanpa nilai nodata
        masked_data = np.where(np.isnan(data), -9999, data).astype(np.float32)

        # Terapkan klasifikasi
        classified = np.vectorize(classify_drought)(masked_data)  # Isi sementara dengan NaN

        # Gantikan semua NaN dalam hasil klasifikasi dengan nilai tetap (misalnya -9999)
        classified = np.where(np.isnan(classified), -9999, classified).astype(np.float32)

    # Perbarui metadata
    meta.update(dtype=rasterio.float32, count=1, nodata=-9999)  # Tetapkan nodata ke -9999

    # Simpan raster hasil
    with rasterio.open('static/output/PDSI_classify.tiff', 'w', **meta) as dst:
        dst.write(classified, 1)
    return("PDSI berhasil dibuat.\n")

def copy_pdsi(year, month, path, row):
    shutil.copy2(f'static/database/PDSI/{path}_{row}/scpdsi_{path}_{row}_{year}-{month}-01.tiff', "static/output/PDSI.tiff")
    return("PDSI berhasil dibuat.\n")

def start_process(year, month, path, row):
    calculate_ndvi()
    calculate_lst()
    calculate_savi()
    calculate_vci(path, row)
    calculate_tci(path, row)
    calculate_vhi()
    calculate_spi(year, month)
    calculate_dem()
    calculate_nst()
    klasifikasi.start_classification()
    if 2013 < int(year) < 2024:
        calculate_pdsi(year, month)
    return 'Klasifikasi Berhasil'

# Get Result Process
def convert_classification(type=True):
    if type:
        raster_path = 'static/Peta_Kekeringan.tiff'  # Path ke file TIFF
        output_image_path = 'static/convert/Peta_Kekeringan.png'
        title = 'Klasifikasi Kekeringan'
    else:
        raster_path = 'static/output/PDSI_classify.tiff'  # Path ke file TIFF
        output_image_path = 'static/convert/PDSI.png'
        title = 'Gambar PDSI'

    with rasterio.open(raster_path) as src:
        data = src.read(1, masked=True)

    # Tentukan warna untuk setiap kategori kekeringan
    colors = ['#8b0000', '#ff4500', '#ffa500', '#ffff00', '#00ff00']  # Warna untuk setiap kategori
    cmap = ListedColormap(colors)

    # Tetapkan batas kelas yang sesuai untuk kategori 0-4
    bounds = [0, 1, 2, 3, 4, 5]  # Setiap batas mencakup satu kategori
    norm = BoundaryNorm(bounds, cmap.N)

    # Plot data raster dengan colormap dan batas kategori yang ditentukan
    plt.figure(figsize=(8, 6))
    img = plt.imshow(data, cmap=cmap, norm=norm)
    # cbar = plt.colorbar(img, ticks=[0, 1, 2, 3, 4], label="Kategori Kekeringan")
    # cbar.ax.set_yticklabels(['Kekeringan Ekstrem', 'Kekeringan Parah', 'Kekeringan Sedang', 
    #                          'Kekeringan Ringan', 'Tidak ada kekeringan (Normal)'])
    plt.title(title)
    plt.axis('off')  # Menyembunyikan sumbu

    # Simpan gambar
    plt.savefig(output_image_path, format='png', bbox_inches='tight', pad_inches=0, dpi=300)
    plt.close()
    
    return output_image_path

def convert_raster(raster_path, output_image_path, cmap='YlOrRd', title='Gambar'):
    with rasterio.open(raster_path) as src:
        # Baca data band pertama (asumsi data kekeringan hanya satu band)
        raster_data = src.read(1)
        
        # Terapkan nilai NoData sebagai NaN agar tidak muncul di gambar
        if src.nodata is not None:
            raster_data = np.where(raster_data == src.nodata, np.nan, raster_data)

    # Atur ukuran gambar berdasarkan ukuran raster
    plt.figure(figsize=(10, 10))
    
    # Tampilkan data raster dengan colormap yang sesuai dan simpan sebagai gambar
    plt.imshow(raster_data, cmap=cmap, vmin=np.nanmin(raster_data), vmax=np.nanmax(raster_data))
    plt.colorbar(label=title)  # Menambahkan colorbar untuk skala
    plt.axis('off')  # Hilangkan sumbu

    # Simpan gambar
    plt.savefig(output_image_path, format='png', bbox_inches='tight', pad_inches=0, dpi=300)
    plt.close()
    
    return output_image_path

def show_uploaded_raster():
    uploaded_url = {}
    uploaded_url["B4"] = convert_raster(f'static/output/Resample/landsat_SR_B4.tiff', 'static/convert/b4.png', 'Greys', "B4")
    uploaded_url["B5"] = convert_raster(f'static/output/Resample/landsat_SR_B5.tiff', 'static/convert/b5.png', 'Greys', "B5")
    uploaded_url["B10"] = convert_raster(f'static/output/Resample/landsat_ST_B10.tiff', 'static/convert/b10.png', 'Greys', "B10")
    return uploaded_url

def show_preprocessing_raster():
    preprocessing_url = {}
    preprocessing_url['NDVI (Normalized Difference Vegetation Index): Indeks yang menggambarkan kesehatan vegetasi berdasarkan perbedaan pantulan cahaya merah dan inframerah-dekat; nilai tinggi menunjukkan vegetasi sehat dan lebat, sedangkan nilai rendah menunjukkan vegetasi jarang atau tidak ada.'] = convert_raster('static/output/NDVI.tiff', 'static/convert/NDVI.png','Greens', 'NDVI')
    preprocessing_url['LST (Land Surface Temperature): Parameter suhu permukaan tanah; nilai tinggi menunjukkan suhu lebih panas yang dapat mengindikasikan kekeringan, sedangkan nilai rendah menunjukkan suhu lebih dingin.'] = convert_raster('static/output/LST.tiff', 'static/convert/LST.png','hot', 'LST')
    preprocessing_url['SAVI (Soil Adjusted Vegetation Index): Indeks vegetasi yang mengoreksi pengaruh pantulan tanah untuk area dengan tutupan vegetasi rendah; nilai tinggi menunjukkan vegetasi sehat, sedangkan nilai rendah menunjukkan lahan gundul atau vegetasi terdegradasi.'] = convert_raster('static/output/SAVI.tiff', 'static/convert/SAVI.png', 'YlGn', 'SAVI')
    preprocessing_url['VCI (Vegetation Condition Index): Indeks yang membandingkan NDVI saat ini dengan nilai maksimum dan minimum historisnya; nilai tinggi menunjukkan kondisi vegetasi baik, sedangkan nilai rendah menunjukkan stres pada vegetasi akibat kekeringan.'] = convert_raster('static/output/VCI.tiff', 'static/convert/VCI.png', 'RdYlGn', 'VCI')
    preprocessing_url['TCI (Temperature Condition Index): Indeks yang membandingkan suhu saat ini dengan nilai maksimum dan minimum historisnya; nilai tinggi menunjukkan suhu lebih dingin dan kondisi lebih baik, sedangkan nilai rendah menunjukkan suhu tinggi dan kondisi yang lebih kering.'] = convert_raster('static/output/TCI.tiff', 'static/convert/TCI.png', 'coolwarm', 'TCI')
    preprocessing_url['VHI (Vegetation Health Index): Indeks gabungan dari VCI dan TCI untuk menilai kesehatan vegetasi; nilai tinggi menunjukkan vegetasi dalam kondisi baik, sedangkan nilai rendah menunjukkan vegetasi tertekan akibat kekeringan.'] = convert_raster('static/output/VHI.tiff', 'static/convert/VHI.png', 'RdYlGn', 'VHI')
    preprocessing_url['Normalize LST: Normalisasi suhu untuk mengurangi variabilitas; nilai tinggi menunjukkan suhu lebih tinggi dari rata-rata normal, sedangkan nilai rendah menunjukkan suhu lebih rendah dari rata-rata normal.'] = convert_raster('static/output/NST.tiff', 'static/convert/NST.png', 'coolwarm', 'NST')
    preprocessing_url['SPI (Standardized Precipitation Index): Indeks yang mengevaluasi tingkat kekeringan berdasarkan curah hujan; nilai tinggi menunjukkan curah hujan melimpah (basah), sedangkan nilai rendah menunjukkan defisit curah hujan (kekeringan).'] = convert_raster('static/output/SPI.tiff', 'static/convert/SPI.png', 'BrBG', 'SPI')
    preprocessing_url['DEM (Digital Elevation Model): Representasi topografi permukaan bumi; nilai tinggi menunjukkan elevasi yang lebih tinggi seperti perbukitan atau pegunungan, sedangkan nilai rendah menunjukkan daerah dataran rendah yang lebih rentan terhadap kekeringan.'] = convert_raster('static/output/DEM.tiff', 'static/convert/DEM.png', 'terrain', 'DEM')
    # preprocessing_url['image'] = convert_classification()
    convert_classification()
    return preprocessing_url

def show_result():
    result_url = {}
    result_url['PDSI'] = convert_classification(False)
    result_url['Klasifikasi'] = convert_classification()
    return result_url

def calculate_drought_category_percentage():
    raster_path = 'static/Peta_Kekeringan.tiff'  # Path ke file TIFF
    # Buka file raster
    with rasterio.open(raster_path) as src:
        data = src.read(1)

    # Daftar kategori kekeringan
    categories = [0, 1, 2, 3, 4]  # 0 = Kekeringan Ekstrem, 4 = Tidak ada kekeringan (Normal)
    category_labels = [0, 1, 2, 3, 4]

    # Hitung total piksel (tidak termasuk NoData jika ada)
    total_pixels = np.count_nonzero(~np.isnan(data)) if src.nodata is not None else data.size

    # Hitung jumlah piksel untuk setiap kategori
    percentages = {}
    for category, label in zip(categories, category_labels):
        count = np.count_nonzero(data == category)
        percentage = (count / total_pixels) * 100
        percentages[label] = round(percentage, 2)

    return percentages

def calculate_average_category():
    # Load categorical raster
    with rasterio.open('static/Peta_Kekeringan.tiff') as category_src:
        category_data = category_src.read(1)

    # Function to calculate average per category for an index raster
    def calculate_average_per_category(index_data, category_data, category_values):
        results = {}
        for category in category_values:
            # Mask the index data with the category
            mask = (category_data == category)
            masked_index_data = np.where(mask, index_data, np.nan)
            
            # Calculate the mean for the masked area and replace NaN with 0
            mean_value = np.nanmean(masked_index_data)
            formatted_mean = round(np.nan_to_num(mean_value, nan=0.0), 2)  # Round to 2 decimal places
            results[category] = formatted_mean
        return results

    # List of category values (0 to 4)
    category_values = [0, 1, 2, 3, 4]

    # Load each index raster and calculate the average per category
    indices = {
        'NDVI': 'static/output/NDVI.tiff',
        'LST': 'static/output/LST.tiff',
        'SAVI': 'static/output/SAVI.tiff',
        'VCI': 'static/output/VCI.tiff',
        'TCI': 'static/output/TCI.tiff',
        'VHI': 'static/output/VHI.tiff',
        'SPI': 'static/output/SPI.tiff',
        'DEM': 'static/output/DEM.tiff',
        'NST': 'static/output/NST.tiff'
    }

    # Dictionary to hold results, structured by category
    results = {category: {} for category in category_values}

    for index_name, path in indices.items():
        with rasterio.open(path) as index_src:
            index_data = index_src.read(1)
        
        # Calculate averages for each category and store in dictionary
        category_averages = calculate_average_per_category(index_data, category_data, category_values)
        
        # Store each index average in the results dictionary under each category
        for category, avg_value in category_averages.items():
            results[category][index_name] = avg_value

    # Display the structured dictionary result
    return results

def get_accuracy_matrix():
    # Fungsi klasifikasi kekeringan berdasarkan nilai PDSI
    def classify_drought(pdsi):
        if pdsi <= -4:
            return 0  # Kekeringan Ekstrem
        elif -4 < pdsi <= -3:
            return 1  # Kekeringan Parah
        elif -3 < pdsi <= -2:
            return 2  # Kekeringan Sedang
        elif -2 < pdsi <= -1:
            return 3  # Kekeringan Ringan
        else:
            return 4  # Tidak ada kekeringan (Normal)

    # Raster prediksi dan ground truth
    pred_raster_path = "static/Peta_Kekeringan.tiff"
    gt_raster_path = "static/output/PDSI.tiff"  # Raster berisi nilai PDSI

    # Baca raster prediksi
    with rasterio.open(pred_raster_path) as pred_src:
        pred = pred_src.read(1)

    # Baca raster ground truth (PDSI)
    with rasterio.open(gt_raster_path) as gt_src:
        gt_pdsi = gt_src.read(1)

    # Pastikan ukuran dan resolusi raster sama
    assert pred.shape == gt_pdsi.shape, "Prediksi dan ground truth harus memiliki dimensi yang sama"

    # Hilangkan piksel nodata (jika ada)
    nodata_value_pred = pred_src.nodata
    nodata_value_gt = gt_src.nodata

    if nodata_value_pred is not None:
        pred[pred == nodata_value_pred] = np.nan
    if nodata_value_gt is not None:
        gt_pdsi[gt_pdsi == nodata_value_gt] = np.nan

    # Filter nilai NaN dari kedua array
    valid_mask = ~np.isnan(pred) & ~np.isnan(gt_pdsi)
    pred = pred[valid_mask]
    gt_pdsi = gt_pdsi[valid_mask]

    # Klasifikasi data ground truth (PDSI) ke dalam kategori kekeringan
    gt_classes = np.array([classify_drought(value) for value in gt_pdsi.flatten()])

    # Pastikan tidak ada NaN dalam prediksi dan ground truth
    valid_mask = ~np.isnan(pred) & ~np.isnan(gt_classes)
    pred = pred[valid_mask]
    gt_classes = gt_classes[valid_mask]

    # Hitung matrix kebingungan dan akurasi
    conf_matrix = confusion_matrix(gt_classes, pred)
    accuracy = accuracy_score(gt_classes, pred)

    # Laporan klasifikasi
    report = classification_report(gt_classes, pred, output_dict=True)

    # Daftar kategori yang diinginkan
    categories = {
        '0': 'Ekstrem',
        '1': 'Parah',
        '2': 'Sedang',
        '3': 'Ringan',
        '4': 'Normal'
    }

    # Buat salinan report baru dengan nama kategori yang diperbarui
    updated_report = {}

    # Iterasi melalui kategori untuk memproses nilai dari report
    for key, label in categories.items():
        if key in report:
            updated_report[label] = report[key]
        elif f"{key}.0" in report:  # Cek key sebagai float string
            updated_report[label] = report[f"{key}.0"]
        else:
            updated_report[label] = {'precision': 0.0, 'recall': 0.0, 'f1-score': 0.0, 'support': 0.0}

    # Tambahkan bagian lainnya (accuracy, macro avg, weighted avg)
    for key in ['accuracy', 'macro avg', 'weighted avg']:
        updated_report[key] = report[key]

    # Buat table akurasi
    report_df = pd.DataFrame(updated_report).transpose()
    report_df.drop(columns=['support'], inplace=True)  # Drop the support column if not needed
    fig, ax = plt.subplots(figsize=(10, len(report_df) // 2))
    sns.heatmap(report_df, annot=True, cmap='Blues', ax=ax, fmt='.2f')
    plt.title('Classification Report')
    plt.savefig('static/convert/Test_Report.png')
    # plt.show()
    plt.close()

    # Buat table matrix
    plt.figure(figsize=(10, 7))
    sns.heatmap(conf_matrix, annot=True, fmt='d', cmap='Blues', xticklabels=['Ekstrem', 'Parah', 'Sedang', 'Ringan', 'Normal'], yticklabels=['Ekstrem', 'Parah', 'Sedang', 'Ringan', 'Normal'])
    plt.title('Confusion Matrix')
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.savefig('static/convert/Test_Matrix.png')
    plt.close()

    return "{:.2f}%".format(accuracy * 100)