// Combine Step3 generate image bands
var ZJK_step3_1 = ee.Image("projects/users/assets/ZJK_step3_1"),
var ZJK_step3_2 = ee.Image("projects/users/assets/ZJK_step3_2"),
var landsat = ZJK_step3_1.addBands(ZJK_step3_2);
var ZJK = 
    /* Define the region of interest (ROI) as a rectangle */
    ee.Geometry.Polygon(
        [[[117.41226581297528, 23.933912422491076],
          [117.41226581297528, 23.921281282091655],
          [117.41818813048016, 23.921281282091655],
          [117.41818813048016, 23.933912422491076]]], null, false);
// Define the region of interest (ROI) using the mangrove2018 dataset
var ROI = ZJK;
Map.addLayer(mangrove2018, {color: 'blue'}, "Geo_intersection");

// Set parameters for time series analysis
var wSize = 30;  // Window size for temporal analysis
var th = 0.5;  // Threshold for binary classification
var point = geometry2;  // Point for plotting charts
var scale = 30;  // Scale for processing
var scaleOut = 5000;  // Output scale
var year = 2020;  // Year for analysis

// Function to prepare time series data by masking NDVI values
var prepareTimeSeries = function(im){
    var nodata = im.mask().not();
    var mask = nodata.mask();
    var ndvi = im.select(['ndvi']);
    var ndviMask = ndvi.lt(0.2);  // Mask NDVI values less than 0.2
    return im.addBands(ndvi.where(ndviMask, 0.2))
             .updateMask(mask)
             .set('system:time_start', im.get('system:time_start'));
};

// Set the processing time range
var begin_date = ee.String('2019-10-01');  // Start date
var end_date = ee.String('2021-05-30');  // End date

// Filter Landsat NDVI data by date and ROI, then sort by time
var landsat_ndvi_mask = landsat_sr.filterDate(begin_date, end_date)
                                   .filterBounds(ROI)
                                   .sort('system:time_start');

// Function to create a mosaic of images by date
function mosaicByDate(imcol){
  var imlist = imcol.toList(imcol.size());
  var unique_dates = imlist.map(function(im){
    return ee.Image(im).date().format("YYYY-MM-dd");
  }).distinct();

  // Create a mosaic for each unique date
  var mosaic_imlist = unique_dates.map(function(d){
    d = ee.Date(d);
    var im = imcol.filterDate(d, d.advance(1, "day")).mosaic();
    return im.set("system:time_start", d.millis(), "system:id", d.format("YYYY-MM-dd"));
  });
  return ee.ImageCollection(mosaic_imlist);
}
var landsat_ndvi_mask = mosaicByDate(landsat_ndvi_mask);

// Get the number of bands in the Landsat image and generate a sequence for iteration
var length = landsat.bandNames().length();
var length_seq = ee.List.sequence(0, length.subtract(1));
var landsat_ndvi_list = landsat_ndvi_mask.toList(length);

// Function to convert each band to an image with NDVI values
function bandToimage(i){
  var bandname = landsat.bandNames().get(i);
  var landsat_image = ee.Image(landsat_ndvi_list.get(i));
  var image = landsat.select([bandname]).divide(10000).toFloat();
  image = image.rename('ndvi');  // Rename the band to 'ndvi'
  return image.set('system:time_start', landsat_image.get('system:time_start'));
}
var landsat_sg = length_seq.map(bandToimage);

// Combine multiple images to create a large composite image for analysis
var landEstarfm = image.addBands(image2).addBands(image3).addBands(image4).addBands(image5).addBands(image6).addBands(image7).addBands(image8).addBands(image9);

// Get the number of bands in the composite image and generate a sequence for iteration
var length0 = landEstarfm.bandNames().length();
var length_seq0 = ee.List.sequence(0, length0.subtract(1));

// Function to convert each band in the composite image to an NDVI image
function bandToimage0(i){
  var bandname = ee.String(landEstarfm.bandNames().get(i));
  var landsat_image = ee.Image(landEstarfm.select(bandname));
  landsat_image = landsat_image.rename('ndvi');  // Rename the band to 'ndvi'
  var system_time_start = bandname.slice(-13);  // Extract the time from the band name
  return landsat_image.set('system:time_start', ee.Number.parse(system_time_start));
}
var landsat_sg0 = length_seq0.map(bandToimage0);

// Combine two image collections and filter by the ROI and sort by time
var collection1 = ee.ImageCollection(landsat_sg);
var collection2 = ee.ImageCollection(landsat_sg0);
var landsat_collection = collection1.merge(collection2)
                                    .sort('system:time_start')
                                    .filterBounds(ROI);

// Filter the collection to a specific time range and apply the preparation function
var TOC = landsat_collection.filterDate(year + '-01-01', (year + 1) + '-05-30')
                             .map(prepareTimeSeries);

// Function to convert images to binary based on a threshold
var covertToBinary = function(im){
  var im2 = im.gt(thresh);
  return im2.set('system:time_start', im.get('system:time_start'));
};

// Function to apply a binary mask to images
var maskBinary = function(im){
  return im.set('system:time_start', im.get('system:time_start'));
};

// Calculate amplitude and threshold for the time series
var bioBandName = 'ndvi';
var amplitude = TOC.select(bioBandName).reduce(ee.Reducer.percentile([95]))
                   .subtract(TOC.select(bioBandName).reduce(ee.Reducer.percentile([5])));
var min = TOC.select(bioBandName).reduce(ee.Reducer.percentile([5]));
var thresh = amplitude.multiply(th).add(TOC.select(bioBandName).reduce(ee.Reducer.percentile([5]))).rename('thresh');
var threshdict = thresh.reduceRegion(ee.Reducer.mean(), point, scale);
print('thresh:', threshdict.get('thresh'));

// Define parameters for temporal analysis
var lag = 1;
var startDate = year + '-01-01';
var endDate = (year + 1) + '-05-30';
var listDates = ee.List.sequence(ee.Date(startDate).millis(), ee.Date(endDate).millis(), 86400000 * lag);

// Estimate the difference in NDVI between two time periods
var out_diff = ee.ImageCollection(listDates.map(function(dd){
  var targetDay = ee.Date(dd);

  var startDate = targetDay.advance(-wSize, 'day');
  var endDate = targetDay;

  var TOC = landsat_collection
        .filterDate(startDate, endDate)
        .filterBounds(ROI)
        .map(prepareTimeSeries)
        .map(covertToBinary);

  var median_before = TOC.select(bioBandName).reduce(ee.Reducer.mean()).rename(bioBandName);

  startDate = endDate;
  endDate = targetDay.advance(wSize, 'day');

  TOC = landsat_collection
        .filterDate(startDate, endDate)
        .filterBounds(ROI)
        .map(prepareTimeSeries)
        .map(covertToBinary);

  var median_after = TOC.select(bioBandName).reduce(ee.Reducer.mean()).rename(bioBandName);

  var diff_ndvi = median_after.subtract(median_before);

  var doy = ee.Image(targetDay.getRelative('day', 'year')).rename('doy');
  var out = diff_ndvi.multiply(-1).rename('diff_nd')
                     .addBands(diff_ndvi.rename('diff_nd_inv'))
                     .addBands(doy.int());

  return out.set('system:time_start', targetDay.millis());
}));

// Determine the Start of Season (SoS) and End of Season (EoS) based on NDVI differences
if (0){
  var EoSdiff = out_diff.max().select('diff_nd').rename('EoS_maxDiff');
  var SoSdiff = out_diff.max().select('diff_nd_inv').rename('SoS_maxDiff');
  
  var out_diff_maskedSoS = out_diff.map(function(im){
    var binSoS = im.select('diff_nd_inv').eq(SoSdiff);
    return im.updateMask(binSoS).copyProperties(im, ['system:time_start']);
  });

  var out_diff_maskedEoS = out_diff.map(function(im){
    var binEoS = im.select('diff_nd').eq(EoSdiff);
    return im.updateMask(binEoS).copyProperties(im, ['system:time_start']);
  });

  var SoS_multi = out_diff_maskedSoS.select('doy').mean().rename('SoS');
  var EoS_multi = out_diff_maskedEoS.select('doy').mean().rename('EoS');

  var countEoS = out_diff_maskedEoS.select('diff_nd').reduce(ee.Reducer.count());
  var countSoS = out_diff_maskedSoS.select('diff_nd_inv').reduce(ee.Reducer.count());

  var SoS = out_diff.qualityMosaic('diff_nd_inv').select('doy').rename('SoS');
  var EoS = out_diff.qualityMosaic('diff_nd').select('doy').rename('EoS');

  var SoS = SoS.where(countSoS.gte(3), SoS_multi);
  var EoS = EoS.where(countEoS.gte(3), EoS_multi);
} else {
  var SoS = out_diff.qualityMosaic('diff_nd_inv').select('doy').rename('SoS');
  var EoS = out_diff.qualityMosaic('diff_nd').select('doy').rename('EoS');
}

// Export the Start of Season (SoS) image to Google Drive
Export.image.toDrive({
  image: SoS.int16(),
  description: 'SoS_ZJK',
  scale: 30,
  region: ROI
});

// Create and print time series charts for NDVI and binary data
var chart1 = ui.Chart.image.series(TOC.select(bioBandName), point, ee.Reducer.mean(), scale)
  .setOptions({title: 'MOD09 NDVI', 
              lineWidth: 0,
              pointSize: 4});
print(chart1);

var chart2 = ui.Chart.image.series(TOC.select(bioBandName).map(covertToBinary), point, ee.Reducer.mean(), scale)
  .setOptions({title: 'Binary time series', 
              lineWidth: 0,
              pointSize: 4});
print(chart2);

var chart3 = ui.Chart.image.series(out_diff.select('diff_nd'), point, ee.Reducer.mean(), scale)
  .setOptions({title: 'Ratio observations > threshold (after-before)', 
              lineWidth: 0,
              pointSize: 4});
print(chart3);

// Adjust End of Season (EoS) day of year (DOY) values to correct for calendar year boundaries
var adjustedDoy = EoS.expression(
  'doy < 100 ? doy + 366 : doy', {
    'doy': EoS
}).rename('EoS');
print('adjustedDoy', adjustedDoy);

// Initialize a date for calculations
var init = ee.Image(ee.Date(year + '-01-01').millis());

// Calculate percentiles for SoS and EoS, then apply a mask for SoS
var SoS25 = SoS.reduceRegion(ee.Reducer.percentile([25]), point, scale).get('SoS');
var SoS75 = SoS.reduceRegion(ee.Reducer.percentile([75]), point, scale).get('SoS');
print('SoS 25:', ee.Date(SoS.reduceRegion(ee.Reducer.percentile([25]), point, scale).get('SoS')));
print('SoS 75:', ee.Date(SoS.reduceRegion(ee.Reducer.percentile([75]), point, scale).get('SoS')));

var mask = SoS.gte(ee.Number(SoS25)).and(SoS.lte(ee.Number(SoS75)));
var maskedSoS = SoS.updateMask(mask);
var mean = maskedSoS.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: point,
  scale: 30
}).get('SoS');
print('SoS mean:', ee.Date(mean));
print('SoS stdDev:', ee.Date(maskedSoS.reduceRegion(ee.Reducer.stdDev(), point, scale).get('SoS')));
var SoSdate = maskedSoS.multiply(86400000).add(init);
print('Start of Season:', ee.Date(SoSdate.reduceRegion(ee.Reducer.mean(), point, scale).get('SoS')));

// Apply a similar mask for EoS
var EoS25 = adjustedDoy.reduceRegion(ee.Reducer.percentile([25]), point, scale).get('EoS');
var EoS75 = adjustedDoy.reduceRegion(ee.Reducer.percentile([75]), point, scale).get('EoS');
print('EoS 25:', ee.Date(adjustedDoy.reduceRegion(ee.Reducer.percentile([10]), point, scale).get('EoS')));
print('EoS 75:', ee.Date(adjustedDoy.reduceRegion(ee.Reducer.percentile([75]), point, scale).get('EoS')));

var mask = adjustedDoy.gte(ee.Number(EoS25)).and(adjustedDoy.lte(ee.Number(EoS75)));
var maskedEoS = adjustedDoy.updateMask(mask);
var mean = maskedEoS.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: point,
  scale: 30
}).get('EoS');
print('EoS stdDev:', ee.Date(maskedEoS.reduceRegion(ee.Reducer.stdDev(), point, scale).get('EoS')));
print('EoS mean:', ee.Date(mean));

var EoSdate = maskedEoS.multiply(86400000).add(init);
print('End of Season:', ee.Date(EoSdate.reduceRegion(ee.Reducer.mean(), point, scale).get('EoS')));

// Calculate the number of days between Start of Season (SoS) and End of Season (EoS)
var startDate = ee.Date(SoSdate.reduceRegion(ee.Reducer.median(), point, scale).get('SoS'));
var endDate = ee.Date(EoSdate.reduceRegion(ee.Reducer.median(), point, scale).get('EoS'));

var daysBetween = endDate.difference(startDate, 'day');
print('Days between:', daysBetween);
