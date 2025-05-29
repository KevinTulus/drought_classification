// JavaScript for handling raster layer toggling
document.addEventListener("DOMContentLoaded", function () {
    // Initialize the map
    var map = L.map('map').setView([-2.5, 118], 5);

    // Add a base map layer
    L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
        maxZoom: 19,
        attribution: '© OpenStreetMap contributors'
    }).addTo(map);

    // URLs to each GeoTIFF file for different drought categories
    var urls = {
        extreme: "static/map/category_0.tif",
        severe: "static/map/category_1.tif",
        moderate: "static/map/category_2.tif",
        mild: "static/map/category_3.tif",
        none: "static/map/category_4.tif"
    };

    // Function to load and create a GeoRasterLayer, excluding NoData values
    async function loadRaster(url, color) {
        const response = await fetch(url);
        const arrayBuffer = await response.arrayBuffer();
        const georaster = await parseGeoraster(arrayBuffer);

        return new GeoRasterLayer({
            georaster: georaster,
            opacity: 0.7,
            pixelValuesToColorFn: values => {
                const pixelValue = values[0];
                if (pixelValue === -9999 || pixelValue === null || isNaN(pixelValue)) {  // Handle NoData values
                    return null;  // Return null to make these pixels transparent
                }
                return color;  // Apply color for valid data pixels
            },
            resolution: 256
        });
    }

    // Initialize variables to store layers and their visibility status
    var layers = {};
    var layerVisibility = {
        extreme: false,
        severe: false,
        moderate: false,
        mild: false,
        none: false
    };

    // Load all layers initially
    async function initializeLayers() {
        layers.extreme = await loadRaster(urls.extreme, '#8b0000'); // Kekeringan Ekstrem
        layers.severe = await loadRaster(urls.severe, '#ff4500');  // Kekeringan Parah
        layers.moderate = await loadRaster(urls.moderate, '#ffa500');  // Kekeringan Sedang
        layers.mild = await loadRaster(urls.mild, '#ffff00');  // Kekeringan Ringan
        layers.none = await loadRaster(urls.none, '#00ff00');  // Tidak ada kekeringan
    }

    initializeLayers();

    // Toggle layer function
    function toggleLayer(layerKey) {
        if (layerVisibility[layerKey]) {
            map.removeLayer(layers[layerKey]);
        } else {
            map.addLayer(layers[layerKey]);
        }
        layerVisibility[layerKey] = !layerVisibility[layerKey];
    }

    // Event listeners for toggling layers
    document.getElementById('toggle-extreme').addEventListener('click', function () {
        toggleLayer('extreme');
    });

    document.getElementById('toggle-severe').addEventListener('click', function () {
        toggleLayer('severe');
    });

    document.getElementById('toggle-moderate').addEventListener('click', function () {
        toggleLayer('moderate');
    });

    document.getElementById('toggle-mild').addEventListener('click', function () {
        toggleLayer('mild');
    });

    document.getElementById('toggle-none').addEventListener('click', function () {
        toggleLayer('none');
    });
});