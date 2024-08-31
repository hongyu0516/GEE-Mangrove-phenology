var mod_sr = ee.ImageCollection("MODIS/006/MOD09Q1"),
var landsat_sr = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2"),
var mangrove2020 = ee.FeatureCollection("projects/users/assets/Mangrove_China_2020"),
var ZJK = 
    /* color: #98ff00 */
    /* displayProperties: [
      {
        "type": "rectangle"
      }
    ] */
    ee.Geometry.Polygon(
        [[[117.40720638888727, 23.935460610367233],
          [117.40720638888727, 23.91663135512823],
          [117.42480168002497, 23.91663135512823],
          [117.42480168002497, 23.935460610367233]]], null, false),
  
// Add a layer representing 2018 mangrove data to the map, colored yellow.
Map.addLayer(mangrove2020, {color: 'yellow'}, "mangrove2020");

// Define a geometric region of interest (ROI) for the analysis.
var geometry = ZJK;

// Function to extract specific bits from an image's QA band.
var getQABits = function(image, start, end) {
    var pattern = 0;
    for (var i = start; i <= end; i++) {
        pattern += 1 << i; // Create a bitmask for the specified bit range.
    }
    return image.bitwiseAnd(pattern).rightShift(start); // Extract the bits and shift them to the right.
};

// Define the time period for the analysis.
var begin_date = ee.String('2019-10-01'); 
var end_date = ee.String('2021-05-30'); 

// Filter the Landsat Surface Reflectance (SR) collection by date and region, calculate NDVI, and apply cloud masking.
var landsat_ndvi_mask = landsat_sr.filterDate(begin_date, end_date)
   .filterBounds(geometry)
   .map(function(img) {
       var img_sub = img;
       var ndvi = img_sub.normalizedDifference(['SR_B5', 'SR_B4']).rename('L_ndvi'); // Calculate NDVI.
       var Q = img_sub.select('QA_PIXEL'); // Select the QA_PIXEL band for cloud masking.
       var cloud2BitMask = (1 << 2);
       var cloud3BitMask = (1 << 3);
       var cloud4BitMask = (1 << 4);
       var mask = Q.bitwiseAnd(cloud2BitMask).eq(0).and(Q.bitwiseAnd(cloud3BitMask).eq(0))
       .and(Q.bitwiseAnd(cloud4BitMask).eq(0)).rename('L_cmask'); // Combine cloud masks.
       return ndvi.addBands(mask) // Return NDVI with the cloud mask band.
       .copyProperties(img_sub, img_sub.propertyNames()); // Copy properties from the original image.
   });

// Function to mosaic images by date, creating a single image for each date.
function mosaicByDate(imcol){
  var imlist = imcol.toList(imcol.size());
  var unique_dates = imlist.map(function(im){
    return ee.Image(im).date().format("YYYY-MM-dd");
  }).distinct(); // Get unique dates in the image collection.
  
  var mosaic_imlist = unique_dates.map(function(d){
    d = ee.Date(d);
    var im = imcol.filterDate(d, d.advance(1, "day")).mosaic(); // Mosaic images taken on the same day.
    return im.set(
      "system:time_start", d.millis(), 
      "system:id", d.format("YYYY-MM-dd")
      );
  });
  return ee.ImageCollection(mosaic_imlist);
}

// Sort Landsat NDVI images by date and mosaic them by day.
var landsat_ndvi_mask_sort = landsat_ndvi_mask.sort('system:time_start');
var landsat_ndvi_mask = mosaicByDate(landsat_ndvi_mask_sort);

// Extract NDVI and cloud mask bands from the Landsat data.
var landsat_ndvi = landsat_ndvi_mask.select('L_ndvi');
var landsat_cmask = landsat_ndvi_mask.select('L_cmask');

// Filter MODIS Surface Reflectance (SR) collection by date and region, calculate NDVI, and apply cloud masking.
var mod_ndvi = mod_sr.filterDate(begin_date, end_date).filterBounds(geometry)
  .map(function(image) {
      var img_sub = image;
      var ndvi = img_sub.normalizedDifference(['sur_refl_b02', 'sur_refl_b01']); // Calculate NDVI.
      var Q = img_sub.select('State'); // Select the 'State' QA band.
      var masknum = getQABits(Q, 0, 2); // Extract cloud-related bits from the 'State' band.
      var mask = masknum.eq(ee.Image(0)).or(masknum.eq(ee.Image(3))); // Mask clear sky pixels.
      var ndvi_masked = ndvi.updateMask(mask); // Apply the cloud mask to NDVI.
      return ndvi_masked.rename(['MOD_NDVI']).copyProperties(img_sub, img_sub.propertyNames()); // Rename and copy properties.
  });
  
// Import a linear interpolation module and perform interpolation on the MODIS NDVI data.
var interl_m = require('users/Yang_Chen/GF-SG:Interpolation_v1');
var frame  = 8 * 4; 
var nodata = -9999; 
var mod_ndvi_interp = interl_m.linearInterp(mod_ndvi, frame, nodata); // Interpolate the MODIS NDVI data.
var mod_ndvi_interp0 = mod_ndvi_interp.select(['MOD_NDVI_INTER']); // Select the interpolated NDVI band.

// Apply Savitzky-Golay filtering to smooth the MODIS NDVI time series.
var sg_filter = require('users/Yang_Chen/GF-SG:SG_filter_v1');
var list_trend_sgCoeff = ee.List([-0.076923102, -1.4901161e-008, 0.062937059, 0.11188812, 0.14685316, 0.16783218,
0.17482519, 0.16783218, 0.14685316, 0.11188812, 0.062937059, -1.4901161e-008, -0.076923102]); 
var list_sg_coeff = ee.List([-0.090909064, 0.060606077, 0.16883118, 0.23376624, 0.25541127, 0.23376624,
0.16883118, 0.060606077, -0.090909064]);  
var mod_ndvi_sg = sg_filter.sg_filter_chen(mod_ndvi_interp0, list_trend_sgCoeff, list_sg_coeff); // Apply the Savitzky-Golay filter.

// Adjust MODIS NDVI data and resample it to match Landsat resolution.
var L_ref = landsat_ndvi.first(); 
var proj = ee.Image(L_ref).select('L_ndvi').projection().crs();
var mod_ndvi_sg30 = mod_ndvi_sg.map(function(img){ 
  return img.multiply(1.023).subtract(0.013).float().resample('bicubic') 
  .copyProperties(img, img.propertyNames()); // Adjust and resample MODIS NDVI.
});

// Combine the Landsat NDVI and MODIS NDVI collections.
var combine_ndvi = landsat_ndvi.merge(mod_ndvi_sg30);
var combine_ndvi0 = landsat_cmask.merge(mod_ndvi_sg30); // Combine NDVI data and masks.

// Calculate the relative day of year (DOY) for Landsat and MODIS images.
var startDate = ee.Date(begin_date);
var l_size = ee.Number(landsat_ndvi.size());
var m_size = ee.Number(mod_ndvi_sg30.size());
var combinendvi_list = combine_ndvi.toList(l_size.add(m_size));
var combinendvi_list0 = combine_ndvi0.toList(l_size.add(m_size));
var landsatndvi_list = combinendvi_list.slice(0, l_size);
var landsatmask_list = combinendvi_list0.slice(0, l_size);
var modndvi_list = combinendvi_list.slice(l_size, l_size.add(m_size));
var list_process = ee.List.sequence(0, l_size.subtract(1));
var ldate = list_process.map(function(i){
    i = ee.Number(i);
    var one = landsatndvi_list.get(i);
    var date = ee.Image(one).date();
    var relDoy = date.difference(startDate, 'day');
    return relDoy.toInt(); // Calculate relative DOY for Landsat images.
});

// Calculate relative DOY for MODIS images.
var list_process = ee.List.sequence(0, m_size.subtract(1));
var mdate = list_process.map(function(i){
    i = ee.Number(i);
    var one = modndvi_list.get(i);
    var date = ee.Image(one).date();
    var relDoy = date.difference(startDate, 'day');
    return relDoy; // Calculate relative DOY for MODIS images.
});

var mdate0 = ee.Array(mdate); // Convert MODIS DOY to an array.

// Match each Landsat image to the closest MODIS image based on relative DOY.
var list_process = ee.List.sequence(0, l_size.subtract(1));    
var lm_loc = list_process.map(function(i){
    i = ee.Number(i);
    var onenum = ee.Number(ldate.get(i));
    var diff = mdate0.subtract(onenum).abs().toList(); // Calculate the difference between Landsat and MODIS DOY.
    var minloc = diff.indexOf(diff.reduce(ee.Reducer.min())); // Find the closest MODIS image.
    return minloc; // Store the index of the closest MODIS image.
});

// Generate initial masks to track valid and invalid MODIS pixels during image matching.
var list_loc_old = ee.List.repeat(1, m_size);
var replace = function(current, previous){
  previous = ee.List(previous);
  var i = ee.Number(current);
  var oneloc = ee.Number(lm_loc.get(i));
  var list_loc0 = previous.set(oneloc, 0);
  return list_loc0;
};
var list_loc0 = ee.List.sequence(0, l_size.subtract(1)).iterate(replace, list_loc_old);
var list_loc = ee.List(list_loc0);

// Pair each Landsat image with the closest MODIS image based on DOY.
var list_process = ee.List.sequence(0, l_size.subtract(1));
var modndviPair_list = list_process.map(function(i){
    i = ee.Number(i);
    var onenum = ee.Number(lm_loc.get(i));
    var onelist = modndvi_list.get(onenum);
    return onelist;
});
var modndviPair = ee.ImageCollection(modndviPair_list);

// Generate a list of indices from 0 to l_size - 1
var list_process = ee.List.sequence(0, l_size.subtract(1));

// Create a list of masked Landsat NDVI images by masking cloudy pixels
var landsatNdviPair_list = list_process.map(function(i) {
  i = ee.Number(i);  // Convert index to a number
  var landsatone = landsatndvi_list.get(i);  // Get the Landsat NDVI image at index i
  var landsatone_img = ee.Image(landsatone);  // Cast to ee.Image
  var maskone = landsatmask_list.get(i);  // Get the corresponding cloud mask
  var maskone_img = ee.Image(maskone);  // Cast to ee.Image
  var landsatone_masked = landsatone_img.updateMask(maskone_img);  // Apply the cloud mask to the NDVI image
  return landsatone_masked;  // Return the masked NDVI image
});

// Convert the list of masked Landsat NDVI images to an ee.ImageCollection
var landsatNdviPair = ee.ImageCollection(landsatNdviPair_list);

// Create a fixed 31x31 kernel with all elements set to 1 for neighborhood analysis
var W0 = ee.List.repeat(1, 31); 
var W = ee.List.repeat(W0, 31);  
var kW = ee.Kernel.fixed(31, 31, W); 

// Generate neighborhood bands for MODIS NDVI pairs
var modndviPair_nei0 = modndviPair.map(function(img){
  var img_nei = img.neighborhoodToBands(kW);  // Extract neighborhood bands using the kernel
  return img_nei;
});

// Generate neighborhood bands for all MODIS NDVI images
var modndvisg30 = ee.ImageCollection(modndvi_list);
var modndvisg30_nei = modndvisg30.map(function(img){
  var img_nei = img.neighborhoodToBands(kW);  // Extract neighborhood bands using the kernel
  return img_nei;
});

// Convert the neighborhood bands of MODIS pairs to a list
var modndviPair_nei_list = modndviPair_nei0.toList(l_size);

// Apply cloud masks to the neighborhood bands of MODIS pairs
var modndviPair_nei_list = list_process.map(function(i){
  i = ee.Number(i);  // Convert index to a number
  var modisone = modndviPair_nei_list.get(i);  // Get the neighborhood bands at index i
  var modisone_img = ee.Image(modisone);  // Cast to ee.Image
  var maskone = landsatmask_list.get(i);  // Get the corresponding cloud mask
  var maskone_img = ee.Image(maskone);  // Cast to ee.Image
  var modisone_masked = modisone_img.updateMask(maskone_img);  // Apply the cloud mask
  return modisone_masked;
});

// Convert the list of masked MODIS neighborhood images to an ee.ImageCollection
var modndviPair_nei = ee.ImageCollection(modndviPair_nei_list);

// Compute the mean of neighborhood bands for MODIS and Landsat NDVI images
var modndviPair_nei_mean = modndviPair_nei.mean();
var landsatNdviPair_mean = landsatNdviPair.mean();

// Calculate the difference from the mean for each MODIS and Landsat neighborhood image
var modndviPair_nei_diffm = modndviPair_nei.map(function(img){
  var img_diff = img.subtract(modndviPair_nei_mean);  // Subtract mean
  return img_diff;
});
var landsatNdviPair_diffm = landsatNdviPair.map(function(img){
  var img_diff = img.subtract(landsatNdviPair_mean);  // Subtract mean
  return img_diff;
});

// Convert the lists of differences from the mean to lists of images
var modndviPair_nei_diffm_list = modndviPair_nei_diffm.toList(l_size);
var landsatNdviPair_diffm_diffm_list = landsatNdviPair_diffm.toList(l_size);

// Multiply corresponding MODIS and Landsat differences from the mean
var modnei_landsat_diff_mul_list = list_process.map(function(i){
  i = ee.Number(i);  // Convert index to a number
  var modnei = modndviPair_nei_diffm_list.get(i);  // Get the MODIS difference
  var modnei_img = ee.Image(modnei);  // Cast to ee.Image
  var landsat = landsatNdviPair_diffm_diffm_list.get(i);  // Get the Landsat difference
  var landsat_img = ee.Image(landsat);  // Cast to ee.Image
  var landsat_modnei_mul = modnei_img.multiply(landsat_img);  // Multiply the differences
  return landsat_modnei_mul;
});

// Convert the list of multiplied images to an ee.ImageCollection
var modnei_landsat_diff_mul = ee.ImageCollection(modnei_landsat_diff_mul_list);

// Sum the multiplied differences
var modnei_landsat_diffmul_sum = modnei_landsat_diff_mul.sum();

// Square each difference image and sum them for both MODIS and Landsat
var modndviPair_nei_diffm_2 = modndviPair_nei_diffm.map(function(img){
  var img_2 = img.multiply(img);  // Square the difference
  return img_2;
}); 
var landsatNdviPair_diffm_2 = landsatNdviPair_diffm.map(function(img){
  var img_2 = img.multiply(img);  // Square the difference
  return img_2;
}); 
var modndviPair_nei_diffm2_sum = modndviPair_nei_diffm_2.sum();
var landsatNdviPair_diffm2_sum = landsatNdviPair_diffm_2.sum();

// Compute the square root of the product of summed squares
var modnei_landsat_diffmul_sqrt = modndviPair_nei_diffm2_sum.multiply(landsatNdviPair_diffm2_sum)
.sqrt();

// Calculate the correlation between MODIS and Landsat NDVI
var modnei_landsat_corre = modnei_landsat_diffmul_sum.divide(modnei_landsat_diffmul_sqrt);

// Define a threshold image for correlation
var threshold_img = ee.Image.constant(0.8);

// Create a mask for valid correlation values above the threshold
var threshold_mask = modnei_landsat_corre.gte(threshold_img);

// Apply the threshold mask to filter out low correlations
var modnei_landsat_corre_valid = modnei_landsat_corre.multiply(threshold_mask);

// Sum all valid correlations
var modnei_landsat_corre_sum = modnei_landsat_corre_valid.reduce(ee.Reducer.sum());

// Normalize the valid correlations by dividing by the sum
var modnei_landsat_corre_normal = modnei_landsat_corre_valid.divide(modnei_landsat_corre_sum);

// Map function to apply the weighted average of MODIS neighborhood values based on the correlation
var modndviPair_match0 = modndviPair_nei0.map(function(img){
  var onemodis_w = img.multiply(modnei_landsat_corre_normal);  // Apply the correlation weight
  var onemodis_normal = onemodis_w.reduce(ee.Reducer.sum());  // Sum the weighted neighborhood values
  return onemodis_normal;
});

// Same as above but applied to all MODIS NDVI images
var modndvisg30_match0 = modndvisg30_nei.map(function(img){
  var onemodis_w = img.multiply(modnei_landsat_corre_normal);  // Apply the correlation weight
  var onemodis_normal = onemodis_w.reduce(ee.Reducer.sum());  // Sum the weighted neighborhood values
  return onemodis_normal;
});

// Define a fixed 21x21 kernel for neighborhood analysis
var wm0 = ee.List.repeat(1, 21); 
var wm = ee.List.repeat(wm0, 21); 

// Calculate the sum of valid correlation masks
var thresholdmask_sum = threshold_mask.reduce(ee.Reducer.sum());

// Create images to flag valid and invalid locations based on the correlation mask
var thresimg = ee.Image.constant(0);
var valid_flag_mask = thresholdmask_sum.neq(thresimg); 
var invalid_flag_mask = thresholdmask_sum.eq(thresimg); 

// Replace invalid MODIS NDVI pixels with the neighborhood mean
var modndviPair_match = modndviPair_match0.map(function(img){
  var meanimg = img.reduceNeighborhood({reducer: ee.Reducer.mean(),
  kernel:ee.Kernel.fixed(21, 21, wm)});  // Calculate the mean using a 21x21 neighborhood
  var newimg = img.multiply(valid_flag_mask).add(meanimg.multiply(invalid_flag_mask));  // Replace invalid pixels with the mean
  return newimg;
});

// Same as above but applied to all MODIS NDVI images
var modndvisg30_match = modndvisg30_match0.map(function(img){
  var meanimg = img.reduceNeighborhood({reducer: ee.Reducer.mean(),
  kernel:ee.Kernel.fixed(21, 21, wm)});  // Calculate the mean using a 21x21 neighborhood
  var newimg = img.multiply(valid_flag_mask).add(meanimg.multiply(invalid_flag_mask));  // Replace invalid pixels with the mean
  return newimg;
});  

// Apply cloud masks to the matched MODIS NDVI images
var modndviPair_match_list = modndviPair_match.toList(l_size);
var modndviPair_matchmask_list = list_process.map(function(i){
  i = ee.Number(i);  // Convert index to a number
  var modisone = modndviPair_match_list.get(i);  // Get the matched MODIS NDVI image
  var modisone_img = ee.Image(modisone);  // Cast to ee.Image
  var maskone = landsatmask_list.get(i);  // Get the corresponding cloud mask
  var maskone_img = ee.Image(maskone);  // Cast to ee.Image
  var modisone_masked = modisone_img.updateMask(maskone_img);  // Apply the cloud mask
  return modisone_masked;
});

// Combine Landsat and matched MODIS images for linear regression
var mergelist = list_process.map(function(i){
  i = ee.Number(i);  // Convert index to a number
  var one = ee.Image(landsatNdviPair_list.get(i)).rename('y');  // Get the Landsat NDVI image
  var one0 = ee.Image(modndviPair_matchmask_list.get(i)).rename('x');  // Get the matched MODIS NDVI image
  var two = one.addBands(one0).addBands(ee.Image.constant(1));  // Combine the images with a constant band
  return two;
});

// Convert the merged list to an ee.ImageCollection
var merge = ee.ImageCollection(mergelist);

// Perform linear regression to derive coefficients
var independents = ee.List(['constant', 'x']);
var dependent = ee.String('y');
var trend = merge.select(independents.add(dependent))
    .reduce(ee.Reducer.linearRegression(independents.length(), 1));

// Extract and flatten the regression coefficients
var coefficients = trend.select('coefficients')
  .arrayProject([0])
  .arrayFlatten([independents]); 

// Extract the slope (b1) and intercept (b0) coefficients
var b1 = coefficients.select('x');
var b0 = coefficients.select('constant');

// Apply linear regression adjustment to MODIS NDVI using the coefficients b1 and b0 from the Landsat-MODIS regression.
var modndvisg30_adjust = modndvisg30_match.map(function(img){
  var newimg = img.multiply(b1).add(b0);
  return newimg;
});

// Combine the adjusted MODIS NDVI with the original Landsat NDVI data into a single collection.
var combine_lm = modndvisg30_adjust.merge(ee.ImageCollection(landsatndvi_list));
var combine_lm_list = combine_lm.toList(l_size.add(m_size));

// Separate the combined collection back into adjusted MODIS and original Landsat NDVI lists.
var modndvisg30_adjust_list = combine_lm_list.slice(0, m_size);
var landsatndvi_list0 = combine_lm_list.slice(m_size, l_size.add(m_size));

// Generate synthetic series by merging adjusted MODIS NDVI with valid Landsat NDVI pixels.
var list_process = ee.List.sequence(0, m_size.subtract(1));
var syn_series_list = list_process.map(function(i){
  i = ee.Number(i);
  var modimg = ee.Image(modndvisg30_adjust_list.get(i));
  var flag = lm_loc.indexOf(i);
  var imgd = ee.Image(modndvi_list.get(i));
  var landsatimg = ee.Image(landsatndvi_list0.get(flag));
  var landsatmask0 = ee.Image(landsatmask_list.get(flag));
  var validone = ee.Number(1).subtract(list_loc.get(i));
  var validimg = ee.Image.constant(validone);
  var landsatmask = landsatmask0.multiply(validimg);
  var relocimg = landsatimg.multiply(landsatmask).add(modimg.subtract(modimg.multiply(landsatmask)));
  return relocimg.rename('MOD_NDVI_INTER').set('system:time_start', imgd.get('system:time_start'));
});

// Convert the synthetic series list to an ImageCollection.
var syn_series = ee.ImageCollection(syn_series_list);

// Apply Savitzky-Golay filtering to the synthetic series to smooth the data.
var list_trend_sgCoeff = ee.List([-0.070588261, -0.011764720, 0.038009040, 0.078733027, 0.11040724, 0.13303168, 0.14660634,
  0.15113123, 0.14660634, 0.13303168, 0.11040724, 0.078733027, 0.038009040, -0.011764720, -0.070588261]); 
var list_sg_coeff = ee.List([0.034965038, -0.12820521, 0.069930017, 0.31468537, 0.41724950, 0.31468537,
  0.069930017, -0.12820521, 0.034965038]); 
var syn_series_sg = sg_filter.sg_filter_chen(syn_series, list_trend_sgCoeff, list_sg_coeff);

// Match smoothed synthetic series to Landsat acquisition dates.
var syn_series_sg_list = syn_series_sg.toList(l_size.add(m_size));
var list_process = ee.List.sequence(0, l_size.subtract(1));
var syn_series_land = list_process.map(function(i){
    var i = ee.Number(i);
    var onenum = ee.Number(lm_loc.get(i));
    var onelist = ee.Image(syn_series_sg_list.get(onenum));
    return onelist.set('system:time_start', onelist.get('system:time_start'));
});

// Convert the matched synthetic series to an ImageCollection.
var landsatCollection = ee.ImageCollection(syn_series_land);

// Convert the image values to 16-bit integers for export.
var syn_series_sg_int = landsatCollection.map(function(img){
  var constant = ee.Image.constant(10000);
  var newimg = img.multiply(constant).toInt16();
  return newimg;
});

// Extract and prepare the final collection for export.
var l_size = syn_series_sg_int.size();
var syn_series_sg_int_list = syn_series_sg_int.toList(l_size);
var syn_series_sg_part_list = syn_series_sg_int_list.slice(0, l_size);
var syn_series_sg_part = ee.ImageCollection(syn_series_sg_part_list);
var syn_series_sg_part0 = syn_series_sg_part.toBands();
var bandname = syn_series_sg_part0.bandNames();
var bandname0 = bandname.map(function(i){
  var namei = ee.String(i);
  var namei0 = ee.String('Landsat_MOD_NDVI_').cat(namei);
  return namei0;
});
syn_series_sg_part0 = syn_series_sg_part0.rename(bandname0);

// Export the final synthetic NDVI series to an Earth Engine asset.
Export.image.toAsset({
  image: syn_series_sg_part0,
  description: "ZJK2020_land",
  assetId: "ZJK2020_land",
  region: geometry,
  scale: 30,
  crs: "EPSG:4326",
  maxPixels: 1e13
});

