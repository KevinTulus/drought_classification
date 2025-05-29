// JavaScript for image zoom functionality
document.addEventListener("DOMContentLoaded", function () {
    const modal = document.getElementById("image-modal");
    const modalImg = document.getElementById("modal-img");
    const closeModal = document.querySelector(".close");

    document.querySelectorAll(".zoomable").forEach((img) => {
        img.addEventListener("click", function () {
            modal.style.display = "block";
            modalImg.src = this.src;
        });
    });

    closeModal.onclick = function () {
        modal.style.display = "none";
    };

    window.onclick = function (event) {
        if (event.target == modal) {
            modal.style.display = "none";
        }
    };
});

// Scroll Function
function scrollFunc(vars) {
    document.addEventListener("DOMContentLoaded", function () {
        // Check for upload-results parameter and scroll to that section
        if (vars == 1) {
            document.getElementById('upload-section').scrollIntoView({
                behavior: 'smooth'
            });
        }

        // Check for classification-results parameter and scroll to that section
        if (vars == 2) {
            document.getElementById('results-section').scrollIntoView({
                behavior: 'smooth'
            });
        }
    });
}

// loading function
document.addEventListener("DOMContentLoaded", function () {
    // Button for Klasifikasi
    const klasifikasiButton = document.getElementById("klasifikasiButton");
    if (klasifikasiButton) {
        klasifikasiButton.addEventListener("click", function (event) {
            const loadingModal = new bootstrap.Modal(document.getElementById("loadingModal"));
            loadingModal.show();
        });
    }

    // Button for Upload
    const uploadButton = document.getElementById("uploadButton");
    if (uploadButton) {
        uploadButton.addEventListener("click", function (event) {
            const uploadLoadingModal = new bootstrap.Modal(document.getElementById("uploadLoadingModal"));
            uploadLoadingModal.show();
        });
    }
});

// Initialize map
var map = L.map('map').setView([-2.5, 118], 5);  // Center over Indonesia

// Add a base map layer (e.g., OpenStreetMap)
L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
    maxZoom: 19,
    attribution: '© OpenStreetMap contributors'
}).addTo(map);

// URL to the GeoTIFF file in the static folder
var urlToGeoTIFF = "/map";
var rasterLayer;  // Variable to store the raster layer

// Function to load the raster layer
async function loadRasterLayer() {
    try {
        const response = await fetch(urlToGeoTIFF);
        if (!response.ok) {
            throw new Error(`Failed to fetch GeoTIFF: ${response.statusText}`);
        }
        const arrayBuffer = await response.arrayBuffer();
        const georaster = await parseGeoraster(arrayBuffer);

        // Create the raster layer using GeoRasterLayer
        rasterLayer = new GeoRasterLayer({
            georaster: georaster,
            opacity: 0.7,
            pixelValuesToColorFn: values => {
                const pixelValue = values[0];
                switch (pixelValue) {
                    case 0:
                        return '#8b0000'; // Kekeringan Ekstrem
                    case 1:
                        return '#ff4500'; // Kekeringan Parah
                    case 2:
                        return '#ffa500'; // Kekeringan Sedang
                    case 3:
                        return '#ffff00'; // Kekeringan Ringan
                    case 4:
                        return '#00ff00'; // Tidak ada kekeringan
                    default:
                        return null; // Transparent for unknown values
                }
            },
            resolution: 256
        });

        // Add layer to the map
        rasterLayer.addTo(map);
        map.fitBounds(rasterLayer.getBounds());
    } catch (error) {
        console.error("Error loading raster layer:", error);
        alert("Error loading raster layer. Please check the console for details.");
    }
}

// Load raster layer initially
loadRasterLayer();

// Variable for raster visibility status
let isRasterVisible = true;

// Toggle function to show/hide raster layer
const toggleRasterButton = document.getElementById('toggle-raster');
if (toggleRasterButton) {
    toggleRasterButton.addEventListener('click', function () {
        if (isRasterVisible) {
            map.removeLayer(rasterLayer);
            this.textContent = "Tampilkan Raster";
        } else {
            rasterLayer.addTo(map);
            this.textContent = "Sembunyikan Raster";
        }
        isRasterVisible = !isRasterVisible;
    });
}

// Overlay layer to show administrative boundaries
var overlayLayer;
var isOverlayVisible = false;  // Visibility status of the overlay

// URL to the GeoJSON file for administrative boundaries
var urlToGeoJSON = "/static/export.geojson";

// Function to add overlay layer to the map
function addOverlayLayer() {
    fetch(urlToGeoJSON)
        .then(response => {
            if (!response.ok) {
                throw new Error(`Failed to fetch GeoJSON: ${response.statusText}`);
            }
            return response.json();
        })
        .then(data => {
            overlayLayer = L.geoJSON(data, {
                filter: function (feature) {
                    return feature.geometry.type === "Polygon" || feature.geometry.type === "MultiPolygon";
                },
                onEachFeature: function (feature, layer) {
                    if (feature.properties && feature.properties.name) {
                        layer.bindPopup(feature.properties.name);
                    }
                },
                style: function () {
                    return { color: "#ff7800", weight: 1 };
                }
            });
            overlayLayer.addTo(map);
        })
        .catch(error => {
            console.error("Error loading GeoJSON:", error);
            alert("Error loading GeoJSON. Please check the console for details.");
        });
}

// Toggle function to show/hide overlay layer
const toggleOverlayButton = document.getElementById('toggle-overlay');
if (toggleOverlayButton) {
    toggleOverlayButton.addEventListener('click', function () {
        if (isOverlayVisible) {
            map.removeLayer(overlayLayer);
            this.textContent = "Tampilkan Batas-batas Administratif";
        } else {
            addOverlayLayer();
            this.textContent = "Sembunyikan Batas-batas Administratif";
        }
        isOverlayVisible = !isOverlayVisible;
    });
}
