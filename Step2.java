var landsat_sr = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2"),
var mod_sr = ee.ImageCollection("MODIS/061/MOD09GA"),
var mangrove2020 = ee.FeatureCollection("projects/users/assets/Mangrove_China_2020"),
var ZJK = 
    /* color: #bf04c2 */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[117.41226581297528, 23.933912422491076],
          [117.41226581297528, 23.921281282091655],
          [117.41818813048016, 23.921281282091655],
          [117.41818813048016, 23.933912422491076]]], null, false)；
// Add Mangrove2018 layer to the map and set it to yellow
Map.addLayer(mangrove2020,{color: 'yellow'},"mangrove2020");

// Set the region of interest (ROI) to QA and add it to the map
var ROI =   ZJK; 
Map.addLayer(ROI,{color: 'blue'},"ROI");

// Define a function to extract QA bits
var getQABits = function(image, start, end) {
    var pattern = 0;
    for (var i = start; i <= end; i++) {
       pattern += 1 << i;  // Generate bitmask using bitwise operations
    }
    return image.bitwiseAnd(pattern).rightShift(start);  // Apply bitmask and shift bits
};

// Set the start and end dates for data processing
var begin_date = ee.String('2019-10-01');
var end_date = ee.String('2021-05-30');

// Calculate MODIS NDVI time series data
var mod_ndvi = mod_sr.filterDate(begin_date, end_date).filterBounds(ROI)
   .map(function(image) {
   var image0 = image;
   var img_sub = image0;
   var ndvi = img_sub.normalizedDifference(['sur_refl_b02', 'sur_refl_b01']);  // Calculate NDVI
   var Q = img_sub.select('state_1km');  // Select the quality assurance (QA) band
   var masknum = getQABits(Q, 0, 2);  // Extract QA bits for cloud masking
   var mask = masknum.eq(ee.Image(0)).or(masknum.eq(ee.Image(3)));  // Create cloud mask
   var ndvi_masked = ndvi.updateMask(mask);  // Apply the mask to NDVI
   return ndvi_masked.rename(['MOD_NDVI']).copyProperties(img_sub, img_sub.propertyNames());});  // Return the masked NDVI image

// Interpolate the MODIS NDVI time series
var interl_m = require('users/Yang_Chen/GF-SG:Interpolation_v1');
var frame  = 8*4;  // Set the interpolation time frame
var nodata = -9999;  // Define the value for missing data
var mod_ndvi_interp = interl_m.linearInterp(mod_ndvi, frame, nodata);  // Perform linear interpolation
var mod_ndvi_interp0 = mod_ndvi_interp.select(['MOD_NDVI_INTER']);  // Select the interpolated NDVI layer

// Smooth the interpolated MODIS NDVI time series using the Savitzky-Golay filter
var sg_filter = require('users/Yang_Chen/GF-SG:SG_filter_v1');
var list_trend_sgCoeff = ee.List([-0.076923102,-1.4901161e-008,0.062937059,0.11188812,0.14685316,0.16783218,
0.17482519,0.16783218,0.14685316,0.11188812,0.062937059,-1.4901161e-008,-0.076923102]);  // Define trend SG filter coefficients
var list_sg_coeff = ee.List([-0.090909064,0.060606077,0.16883118,0.23376624,0.25541127,0.23376624,
0.16883118,0.060606077,-0.090909064]);  // Define smoothing SG filter coefficients
var mod_ndvi_sg = sg_filter.sg_filter_chen(mod_ndvi_interp0,list_trend_sgCoeff,list_sg_coeff);  // Apply the SG filter

// Scale and resample the smoothed MODIS NDVI data
var mod_ndvi_sg30 = mod_ndvi_sg.map(function(img){ 
	return img.multiply(1.023).subtract(0.013).float().resample('bicubic')  // Scale and resample the image
	.copyProperties(img, img.propertyNames());});

// Define the start date and image count for the time series
var startDate = ee.Date(begin_date);
var m_size = ee.Number(mod_ndvi_sg30.size());

// Output the smoothed and resampled NDVI time series
var syn_series_sg_int = mod_ndvi_sg30.map(function(img){
  var constant = ee.Image.constant(10000);  // Create a constant image
  var newimg = img.multiply(constant).toInt16();  // Multiply NDVI by the constant and convert to 16-bit integer
  return newimg});

// Select the time range and convert the image collection to a multi-band image
var syn_series_sg_int_list = syn_series_sg_int.toList(m_size);
var syn_series_sg_part_list = syn_series_sg_int_list.slice(0,m_size);
var syn_series_sg_part = ee.ImageCollection(syn_series_sg_part_list);
var syn_series_sg_part0 = syn_series_sg_part.toBands();

// Rename the bands of the output image
var bandname = syn_series_sg_part0.bandNames();
var bandname0 = bandname.map(function(i){
  var namei = ee.String(i);
  var namei0 = ee.String('Landsat_MOD_NDVI_').cat(namei);  // Add a prefix to the band names
  return namei0;
});
syn_series_sg_part0 = syn_series_sg_part0.rename(bandname0);

// Export the result to the user's asset
Export.image.toAsset({
  image:syn_series_sg_part0,  // Select the image to export
  description: "ZJK2020_mod",  // Set the export task description
  assetId:"ZJK2020_mod",  // Specify the asset ID for export
  region:ROI,  // Define the region to export
  scale:30,  // Set the spatial resolution
  crs:"EPSG:4326",  // Define the coordinate reference system
  maxPixels:1e13  // Set the maximum number of pixels
});
