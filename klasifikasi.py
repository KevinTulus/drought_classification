# %%
import os
import numpy as np
import rasterio
from xgboost import XGBClassifier

# %%
# Fungsi untuk membaca data raster beserta mask NoData
def read_raster_to_array(raster_path):
    with rasterio.open(raster_path) as src:
        array = src.read(1)  # Membaca band pertama saja
        nodata_value = src.nodata  # Mendapatkan nilai NoData
        mask = array == nodata_value  # Membuat mask untuk NoData
        array = array.astype('float32')  # Ubah tipe array menjadi float
        array[array == nodata_value] = np.nan  # Ganti NoData dengan NaN
    return array, mask

# Fungsi untuk membaca semua data raster dalam urutan tertentu
def stack_rasters_with_nodata(folder_path, expected_order):
    # Buat path ke setiap file berdasarkan urutan yang diharapkan
    raster_files = [os.path.join(folder_path, f"{name}.tiff") for name in expected_order]
    arrays = []
    nodata_masks = []
    
    for raster_file in raster_files:
        array, mask = read_raster_to_array(raster_file)
        arrays.append(array)
        nodata_masks.append(mask)
    
    stacked_arrays = np.dstack(arrays)  # Menggabungkan semua array menjadi stack
    combined_nodata_mask = np.logical_or.reduce(nodata_masks)  # Gabungkan mask NoData
    
    return stacked_arrays, combined_nodata_mask

def start_classification():
    # Tentukan urutan data yang diharapkan
    expected_order = ["NDVI", "LST", "SAVI", "VCI", "TCI", "VHI", "SPI", "DEM", "NST"]

    # Membaca data raster dalam urutan yang ditentukan
    folder = 'static/output'
    raster_stack, nodata_mask = stack_rasters_with_nodata(folder, expected_order)

    # Merapikan data (Flatten untuk mempermudah prediksi)
    raster_shape = raster_stack.shape  # Menyimpan shape asli untuk nanti
    raster_flat = raster_stack.reshape(-1, raster_shape[2])  # Flatten menjadi 2D array

    # Model yang sudah dilatih sebelumnya
    model = XGBClassifier()  # Load model Anda jika belum ada, pastikan model sudah dilatih
    model.load_model('static/xgboost_klasifikasi_kekeringan_VHI.json')  # Load model XGBoost yang sudah dilatih

    # Prediksi menggunakan model
    predictions_flat = model.predict(raster_flat)

    # Mengembalikan hasil prediksi ke bentuk asli (2D)
    predictions = predictions_flat.reshape(raster_shape[0], raster_shape[1])

    # Ubah tipe data predictions menjadi float untuk mendukung NaN
    predictions = predictions.astype('float32')

    # Set kembali nilai NoData pada hasil prediksi
    predictions[nodata_mask] = np.nan  # Atur NoData di area yang sesuai

    # Simpan hasil prediksi sebagai raster baru
    with rasterio.open('static/output/VHI.tiff') as src:
        meta = src.meta

    # Update metadata untuk file output
    meta.update(count=1, dtype=rasterio.float32, nodata=np.nan)  # Set NoData pada metadata

    # Menyimpan raster hasil prediksi dengan NoData
    with rasterio.open('static/Peta_Kekeringan.tiff', 'w', **meta) as dst:
        dst.write(predictions.astype(rasterio.float32), 1)

    return ("Klasifikasi berhasil dan peta prediksi kekeringan telah disimpan.")