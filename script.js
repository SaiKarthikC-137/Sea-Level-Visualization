// Initialize the map with options
var map = L.map('map', {
    worldCopyJump: false, // Disable map wrapping
    rotate: false, // Disable map rotation
    maxBounds: [ // Restrict map to world bounds
        [-90, -180], // Southwest corner
        [90, 180] // Northeast corner
    ],
    maxBoundsViscosity: 1.0 // Strictly enforce bounds
}).setView([0, 0], 2);

// Add a tile layer (base map)
var baseLayer = L.tileLayer('https://{s}.tile.openstreetmap.org/{z}/{x}/{y}.png', {
    attribution: ''
}).addTo(map);

// Define available years
const YEARS = [2020, 2021, 2022, 2023, 2024];
const DEFAULT_YEAR = YEARS[YEARS.length - 1]; // Default to the latest year

// Initialize the sea level layer with a default year
var seaLevelLayer = L.tileLayer(`tiles_${DEFAULT_YEAR}/{z}/{x}/{y}.png`, {
    opacity: 0.7,
    attribution: `Sea Level Data (${DEFAULT_YEAR})`, // Initial attribution
    tms: true // Assuming your tiles are TMS. Change if not.
}).addTo(map);

// Add layer control
var layerControl = L.control.layers({
    'Base Map': baseLayer,
}, {
    'Sea Level Data': seaLevelLayer // Generic name, content will change
}).addTo(map);

// Add to your existing Leaflet code (createGradientLayer - not used for tiles but kept if needed elsewhere)
function createGradientLayer(data) {
    var heatPoints = data.map(function(row) { /* ... */ });
    var gradient = { /* ... */ };
    var heatLayer = L.heatLayer(heatPoints, { /* ... */ });
    return heatLayer;
}

// Create a marker cluster group
var markers = L.markerClusterGroup({
    maxClusterRadius: 30,
    disableClusteringAtZoom: 8,
    spiderfyOnMaxZoom: true,
    showCoverageOnHover: false,
    zoomToBoundsOnClick: true,
});

// Load CSV data for markers
Papa.parse("data/tide_gauges.csv", {
    download: true,
    header: true,
    dynamicTyping: true,
    complete: function(results) {
        var data = results.data;
        displayTideGauges(data);
    }
});

// Function to display tide gauges on the map
function displayTideGauges(data) {
    data.forEach(function(row) {
        var lat = row['Lat'];
        var lon = row['Lon'];
        var name = row['Station Name'] || 'Unnamed Station';
        var sealevel = row['Sea Level'];
        var circleMarker = L.circleMarker([lat, lon], {
            radius: 5, fillColor: "#0077ff", color: "#ffffff",
            weight: 1, opacity: 1, fillOpacity: 0.8
        }).bindPopup(`<b>${name}</b><br>Sea Level (from CSV): ${sealevel}<br>Lat: ${lat}, Lon: ${lon}`);
        markers.addLayer(circleMarker);
    });
    map.addLayer(markers);
}

// --- Dashboard Elements and Logic ---
const yearSelect = document.getElementById('yearSelect');
const year1CompareSelect = document.getElementById('year1CompareSelect');
const year2CompareSelect = document.getElementById('year2CompareSelect');
const updateCompareBtn = document.getElementById('updateCompareBtn');
const comparisonHint = document.getElementById('comparisonHint');

// Helper function to populate a select dropdown with years
function populateYearDropdown(selectElement, yearsArray, defaultValue) {
    yearsArray.forEach(y => {
        const opt = document.createElement('option');
        opt.value = y;
        opt.textContent = y;
        selectElement.appendChild(opt);
    });
    if (defaultValue !== undefined) {
        selectElement.value = defaultValue;
    }
}

// 1. Populate the year dropdowns
populateYearDropdown(yearSelect, YEARS, DEFAULT_YEAR);
populateYearDropdown(year1CompareSelect, YEARS, YEARS[0]); // Default to earliest year for year1
populateYearDropdown(year2CompareSelect, YEARS, DEFAULT_YEAR); // Default to latest year for year2



// 2. Event listener for single year selection
yearSelect.addEventListener('change', e => {
    const year = +e.target.value;
    updateMapForSingleYear(year);
});

// 3. Event listener for comparison update button
updateCompareBtn.addEventListener('click', () => {
    const year1 = +year1CompareSelect.value;
    const year2 = +year2CompareSelect.value;

    if (year1 === year2) {
        comparisonHint.textContent = "Please select two different years for comparison.";
        comparisonHint.style.color = "red";
        return;
    }
    // Enforce year2 > year1 if your tiles_year2_vs_year1 strictly means year2 is newer
    if (year1 >= year2) {
        comparisonHint.textContent = "'Compare Year (Year 2)' should be later than 'With Year (Year 1)'.";
        comparisonHint.style.color = "orange";
        // You could choose to proceed or return here
    } else {
        comparisonHint.textContent = ""; // Clear previous hint
    }
    comparisonHint.textContent = ""; // Clear previous hint if valid

    updateMapForComparison(year2, year1); // Order according to folder structure: tiles_year2_vs_year1
});


// Function to update map for single year view
function updateMapForSingleYear(year) {
    console.log('User picked single year:', year);
    const newTileUrl = `tiles_${year}/{z}/{x}/{y}.png`;
    seaLevelLayer.setUrl(newTileUrl);
    seaLevelLayer.options.attribution = `Sea Level Data (${year})`; // Update attribution
    // If you update the layer name in the control, you'd need to remove and re-add the control or find a way to update it.
    // For simplicity, keeping 'Sea Level Data' generic.
}

// Function to update map for comparison view
function updateMapForComparison(year2, year1) { // year2 is the primary, year1 is the baseline
    console.log(`User comparing year: ${year2} vs ${year1}`);
    const newTileUrl = `comparison_data/tiles_${year2}_vs_${year1}/{z}/{x}/{y}.png`;
    seaLevelLayer.setUrl(newTileUrl);
    seaLevelLayer.options.attribution = `Sea Level Comparison (${year2} vs ${year1})`; // Update attribution
}

// Initialize map view based on default single year selection
updateMapForSingleYear(DEFAULT_YEAR);