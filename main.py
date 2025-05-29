from flask import Flask, request, render_template, send_file, session, flash
import os
import server as ml

app = Flask(__name__)
app.secret_key = 'dev'
app.config['SECRET KEY'] = 'dev'
app.config['OUTPUT_FOLDER'] = 'static/output/'

@app.route('/')
def index():
    return render_template("index.html")

@app.route('/test')
def test():
    return render_template("testing.html")

@app.route('/upload', methods=['POST'])
def upload_file():
    is_error = False

    # Input data
    session['year']  = request.form['year']
    session['month']  = request.form['month']
    session['wrs_path']  = request.form['wrs_path']
    session['wrs_row']  = request.form['wrs_row']  

    # Download raster
    error = ml.start_download(session['year'], session['month'], session['wrs_path'], session['wrs_row'])
    if error:
        is_error = True
        flash(f"Tidak ada citra memenuhi syarat untuk bulan {session['year']}-{session['month']}. Silahkan pilih lokasi atau waktu yang berbeda")
    else:
        # Check raster
        error, messages = ml.check_rasters_alignment()
        if error:
            is_error = True
            for message in messages:
                flash(message)

    # Render template
    if is_error == True:
        return render_template("index.html", scroll=1)
    else:
        uploaded_url = ml.show_uploaded_raster()
        metadata = ml.get_raster_metadata()
        return render_template("index.html",
                               uploaded_url=uploaded_url,
                               metadata=metadata,
                               show_raster=True,
                               scroll=1)

@app.route('/klasifikasi')
def klasifikasi():
    # Get uploaded raster and metadata
    uploaded_url = ml.show_uploaded_raster()
    metadata = ml.get_raster_metadata()

    # Start preprocessing and classification
    ml.start_process(session['year'], session['month'], session['wrs_path'], session['wrs_row'])

    # Get preprocessing raster
    preprocessing_url = ml.show_preprocessing_raster()
    
    # Get percentage classification
    percentage = ml.calculate_drought_category_percentage()

    # Get averages per classification 
    averages = ml.calculate_average_category()

    return render_template("index.html",
                           uploaded_url=uploaded_url,
                           metadata=metadata,
                           preprocessing_url=preprocessing_url,
                           percentage=percentage,
                           averages=averages,
                           show_raster=True, 
                           show_preprocessing=True,
                           scroll=2)

@app.route('/uploadtesting', methods=['POST'])
def uploadtesting():
    is_error = False

    # Input data
    session['year']  = request.form['year']
    session['month']  = request.form['month']
    session['wrs_path']  = request.form['wrs_path']
    session['wrs_row']  = request.form['wrs_row']

    # Download raster
    error = ml.start_download(session['year'], session['month'], session['wrs_path'], session['wrs_row'])
    if error:
        is_error = True
        flash(f"Tidak ada citra memenuhi syarat untuk bulan {session['year']}-{session['month']}. Silahkan pilih lokasi atau waktu yang berbeda")
    else:
        # Check raster
        error, messages = ml.check_rasters_alignment()
        if error:
            is_error = True
            for message in messages:
                flash(message)

    # Render template
    if is_error == True:
        return render_template("testing.html", scroll=1)
    else:
        uploaded_url = ml.show_uploaded_raster()
        metadata = ml.get_raster_metadata()
        return render_template("testing.html",
                               uploaded_url=uploaded_url,
                               metadata=metadata,
                               show_raster=True,
                               scroll=1)

@app.route('/testing')
def testing():
    # Get uploaded raster and metadata
    uploaded_url = ml.show_uploaded_raster()
    metadata = ml.get_raster_metadata()

    # Start preprocessing and classification
    ml.start_process(session['year'], session['month'], session['wrs_path'], session['wrs_row'])

    # Get preprocessing raster
    preprocessing_url = ml.show_preprocessing_raster()

    # Get preprocessing raster
    result_url = ml.show_result()
    
    # Get percentage classification
    percentage = ml.calculate_drought_category_percentage()

    # Get result
    accuracy = ml.get_accuracy_matrix()

    return render_template("testing.html",
                           uploaded_url=uploaded_url,
                           metadata=metadata,
                           preprocessing_url=preprocessing_url,
                           result_url=result_url,
                           percentage=percentage,
                           accuracy=accuracy,
                           show_raster=True, 
                           show_preprocessing=True,
                           scroll=2)

@app.route('/map')
def serve_raster():
    return send_file('static/Peta_Kekeringan.tiff')

@app.route('/download/tiff')
def download_tiff():
    return send_file('static/Peta_Kekeringan.tiff', as_attachment=True)

@app.route('/download/png')
def download_png():
    return send_file('static/convert/Peta_Kekeringan.png', as_attachment=True)

if __name__ == '__main__':
    app.run(debug=True)