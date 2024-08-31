// Import necessary datasets: MODIS surface reflectance, Landsat surface reflectance, and mangrove data from 2018
var mod_sr = ee.ImageCollection("MODIS/006/MOD09GA"),
var landsat_sr = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2"),
var mangrove2020 = ee.FeatureCollection("projects/users/assets/Mangrove_China_2020"),
var ZJK_land = ee.Image("projects/users/assets/ZJK2020_land"),
var ZJK_mod = ee.Image("projects/users/assets/ZJK2020_mod"),
var ZJK = 
    /* Define the region of interest (ROI) as a rectangle */
    ee.Geometry.Polygon(
        [[[117.41226581297528, 23.933912422491076],
          [117.41226581297528, 23.921281282091655],
          [117.41818813048016, 23.921281282091655],
          [117.41818813048016, 23.933912422491076]]], null, false);

// Define the time range for the analysis
var begin_date = ee.String('2019-10-01');  
var end_date = ee.String('2021-05-30'); 
var ROI = ZJK;  // Set the region of interest


// Add the mangrove 2018 layer and ROI to the map for visualization
Map.addLayer(mangrove2020, {color: 'blue'}, "mangrove2020");


// Combine multiple QA (Quality Assurance) bands for Landsat and MODIS
var landsat = ZJK2020_land;
var modis = ZJK2020_mod;

// Filter Landsat NDVI data by date and region
var landsat_ndvi_mask = landsat_sr.filterDate(begin_date, end_date)
  .filterBounds(Geo_intersection)
  .sort('system:time_start');  // Sort by acquisition time

// Function to create a mosaic of images by date
function mosaicByDate(imcol){
  var imlist = imcol.toList(imcol.size());  // Convert ImageCollection to a list
  var unique_dates = imlist.map(function(im){
    return ee.Image(im).date().format("YYYY-MM-dd");
  }).distinct();  // Get distinct dates
  
  // Create a mosaic for each unique date
  var mosaic_imlist = unique_dates.map(function(d){
    d = ee.Date(d);
    var im = imcol.filterDate(d, d.advance(1, "day")).mosaic();
    return im.set("system:time_start", d.millis(), "system:id", d.format("YYYY-MM-dd"));
  });
  return ee.ImageCollection(mosaic_imlist);  // Return as ImageCollection
}
var landsat_ndvi_mask = mosaicByDate(landsat_ndvi_mask);  // Apply the mosaic function

// Get the length of the band names list
var length = landsat.bandNames().length();
var length_seq = ee.List.sequence(0, length.subtract(1));
var landsat_ndvi_list = landsat_ndvi_mask.toList(length);

// Function to convert each band to an image with scaled values
function bandToimage(i){
  var bandname = landsat.bandNames().get(i);
  var landsat_image = ee.Image(landsat_ndvi_list.get(i));
  var image = landsat.select([bandname]).divide(10000).toFloat();  // Scale and convert to float
  return image.set('system:time_start', landsat_image.get('system:time_start'));
}
var landsat_sg = length_seq.map(bandToimage);  // Apply the band conversion function

// Define new time range for MODIS data
var begin_date = ee.String('2019-10-01'); 
var end_date = ee.String('2021-05-30');

// Filter MODIS NDVI data by date and region
var modis_ndvi_mask = mod_sr.filterDate(begin_date, end_date)
  .filterBounds(Geo_intersection);

// Get the length of the MODIS band names list
var length_modis = modis.bandNames().length();
var length_seq_modis = ee.List.sequence(0, length_modis.subtract(1));
var modis_ndvi_list = modis_ndvi_mask.toList(length_modis);

// Function to convert MODIS bands to images with scaled values
function bandToimage_modis(i){
  var bandname = modis.bandNames().get(i);
  var modis_image = ee.Image(modis_ndvi_list.get(i));
  var image = modis.select([bandname]).divide(10000).toFloat();  // Scale and convert to float
  return image.set('system:time_start', modis_image.get('system:time_start'));
}
var modis_sg = length_seq_modis.map(bandToimage_modis);  // Apply the band conversion function for MODIS

// Define the analysis period and region for both Landsat and MODIS
// Since the amount of daily interpolation calculation is too large, the total time can be divided into n segments (generally n equals 5).
// interpolation is performed in each time period, and finally the running results are merged.
var startDate = '2021-01-01';
var endDate = '2021-05-30';
var region = Geo_intersection;  // Set the region for analysis
var commonBandNames = ee.List(['ndvi']);  // Common band names for both datasets
var kernelRadius = ee.Number(10);  // Radius for the square kernel
var kernel = ee.Kernel.square(kernelRadius);  // Define the square kernel
var numPixels = kernelRadius.add(kernelRadius.add(1)).pow(2);  // Number of pixels in the kernel
var coverClasses = 7;  // Number of cover classes

// Function to rename NDVI bands
function reName(image){
    var ndvi = image.rename('ndvi');  // Rename band to 'ndvi'
    return ndvi;
}

// Function to set dates and day of the year (DOY) for images
function DateImage(image){
    return image.setMulti({
        'system:time_start': ee.Date(image.date().format('y-M-d')).millis(),
        'DOY': image.date().format('D')
    });
}

// Process Landsat data by filtering date and region, and apply renaming and date setting functions
var landsat = ee.ImageCollection(landsat_sg)
                    .filterDate(startDate, endDate)
                    .filterBounds(region)
                    .map(reName)
                    .map(DateImage);

// Process MODIS data similarly
var modis = ee.ImageCollection(modis_sg)
              .filterDate(startDate, endDate)
              .map(reName)
              .map(function(image){ 
                  return image.set('DOY', image.date().format('D'));
              });

// Create filters to join Landsat and MODIS images based on date
var dayfilter = ee.Filter.equals({leftField: 'system:time_start', rightField: 'system:time_start'});
var pairedJoin = ee.Join.simple();  // Simple join for paired images
var invertedJoin = ee.Join.inverted();  // Inverted join for unpaired images

// Apply joins to create paired and unpaired image collections
var landsatPaired = pairedJoin.apply(landsat, modis, dayfilter);
var modisPaired = pairedJoin.apply(modis, landsat, dayfilter);
var modisUnpaired = invertedJoin.apply(modis, landsat, dayfilter);
var paired = [landsatPaired, modisPaired, modisUnpaired];

// Function to get a list of dates from images
function getDates(image, empty_list){
    var date = ee.Image(image).date().format('yyyy-MM-dd');
    var updatelist = ee.List(empty_list).add(date);
    return updatelist;  
}

// Function to create sub-collections of paired images based on dates
function makeSubcollections(paired){
  function getSub(ind){
        var lan_01 = paired[0]
            .filterDate(ee.List(dateList).get(ind),
                        ee.Date(ee.List(dateList).get(ee.Number(ind).add(1)))
                            .advance(1, 'day'))
            .toList(2);
            
        var mod_01 = paired[1]
            .filterDate(ee.List(dateList).get(ind),
                        ee.Date(ee.List(dateList).get(ee.Number(ind).add(1)))
                            .advance(1, 'day'))
            .toList(2);
            
        var mod_p = paired[2]
            .filterDate(ee.List(dateList).get(ind),
                        ee.Date(ee.List(dateList).get(ee.Number(ind).add(1)))
                            .advance(1, 'day'));

        var mod_p = mod_p.toList(mod_p.size());

        var subcollection = ee.List([lan_01, mod_01, mod_p]);

        return subcollection;
      }

    // Generate the list of dates for sub-collections
    var empty_list = ee.List([]);
    var dateList = paired[0].iterate(getDates, empty_list);
    var subcols = ee.List.sequence(0, ee.List(dateList).length().subtract(2))
        .map(getSub);
    return subcols;
}
var subs = makeSubcollections(paired);  // Generate sub-collections
var num_lists = ee.List.sequence(0, subs.length().subtract(1));

// Function to process sub-collections into numerical lists for further processing
function tonumlist(i){
  var num_imgs = ee.List(ee.List(subs.get(i)).get(2)).length();
  var remaining = num_imgs.mod(10);
  var index_seq = ee.List.sequence(0, num_imgs.subtract(remaining), 10);
  var subList = ee.List(ee.List(subs.get(i)).get(2)); 
  var list_index_seq = ee.List.sequence(0, index_seq.length().subtract(1))
  .map(function(x){
    var start = ee.Number(index_seq.get(x));
    var end = ee.Algorithms.If(start.add(10).gt(num_imgs), num_imgs, start.add(10));
    var pred_group = subList.slice(start, end);
    var landsat_t01 = ee.List(ee.List(subs.get(i)).get(0));
    var modis_t01 = ee.List(ee.List(subs.get(i)).get(1));
    var modis_tp = pred_group;
    var startDay = ee.Number.parse(ee.ImageCollection(pred_group).first().get('DOY'));
    var endDay = ee.Number.parse(ee.ImageCollection(pred_group).sort('system:time_start', false).first().get('DOY'));
    var year = ee.Date(ee.ImageCollection(pred_group).sort('system:time_start', false).first().get('system:time_start')).format('Y');
    var doys = landsat_t01.map(function(img){ return ee.String(ee.Image(img).get('DOY')).cat('_')});
    
    // Function to spatially align images
    function registerImages(landsat_t01, modis_t01, modis_tp){
        var landsat_t01_1 = landsat_t01.map(function(image){ return ee.Image(image).resample('bicubic');});
        var modis_t01_1 = modis_t01.map(function(image){ return ee.Image(image).resample('bicubic');});
        var modis_tp_1 = modis_tp.map(function(image){ return ee.Image(image).resample('bicubic');});

        // Register MODIS images to Landsat reference images
        var modis_t01_2 = modis_t01_1.map(function(image){ 
          return ee.Image(image).register({
            referenceImage: ee.Image(landsat_t01_1.get(0)),
            maxOffset: ee.Number(150.0),
            stiffness: ee.Number(7.0)
          });
        });
    
        var modis_tp_2 = modis_tp_1.map(function(image){
          return ee.Image(image).register({
            referenceImage: ee.Image(landsat_t01_1.get(0)),
            maxOffset: ee.Number(150.0),
            stiffness: ee.Number(7.0)
          });
        });
        return [landsat_t01_1, modis_t01_2, modis_tp];
    }
        
    var registerImages_1 = registerImages(landsat_t01, modis_t01, modis_tp);
    var landsat_t01 = registerImages_1[0];
    var modis_t01 = registerImages_1[1];
    var modis_tp = registerImages_1[2];

    // Function to prepare Landsat data for further processing
    function prepLandsat(landsat_t01, kernel, numPixels, commonBandNames, doys, coverClasses){
        var neighLandsat_t01 = landsat_t01.map(function(image){ 
          return ee.Image(image).neighborhoodToBands(kernel);  // Apply kernel to image
        });

        var pixPositions = ee.Image(neighLandsat_t01.get(0)).bandNames()
            .map(function(bn){ return ee.String(bn).replace('[a-z]+_', '_');})
            .slice(0, numPixels);

        var pixelBandNames = pixPositions.map(function(position){ 
            return doys.map(function(doy){ 
            return commonBandNames.map(function(bn){ 
            return ee.String(doy).cat('_').cat(ee.String(bn)).cat(ee.String(position));
          });
          });
          }).map(function(l){ return ee.List(l).flatten();});    

        var lanArr_t01 = ee.Image(neighLandsat_t01.get(0)).toArray()
            .arrayCat(ee.Image(neighLandsat_t01.get(1)).toArray(), 0);

        var pixArrays = ee.List([]);
        var pixArrays = ee.List.sequence(0, numPixels.subtract(1)).map(function(i){ 
          return lanArr_t01.arraySlice(0, ee.Number(i).int(),
                                    numPixels.multiply(commonBandNames.length().multiply(2)).int(),
                                    numPixels);
        });

        var lanSorted = ee.List.sequence(0, numPixels.subtract(1))
            .map(function(i){ return ee.Image(pixArrays.get(i)).arrayFlatten([pixelBandNames.get(i)]);});

        // Function to calculate threshold values for pixel masking
        function threshold(landsat, coverClasses){
            function getThresh(image){
                var stddev = ee.Image(image).reduceRegion({
                    reducer: ee.Reducer.stdDev(),
                    bestEffort: true,
                    maxPixels: ee.Number(1e6)
                });
                var stddev = stddev.toImage();
                var names = stddev.bandNames().map(function(bn){return ee.String(bn).cat('_thresh');});
                var thresh = stddev.multiply(ee.Image.constant(ee.Number(2)).divide(coverClasses));
                var thresh = thresh.rename(names);
                return thresh;
            }

            var threshs = ee.List(landsat).map(getThresh);
            return threshs;
        }

        var thresh = threshold(landsat_t01, coverClasses);

        // Apply the calculated thresholds to mask pixels
        function threshMask(neighLandsat_t01, thresh, commonBandNames){
          var masks = ee.List([0, 1]).map(function(i){
            return commonBandNames.map(function(name){
              return ee.Image(neighLandsat_t01.get(i))
                .select([ee.String(name).cat('_(.+)')])
                .select([ee.String(name).cat('_0_0')])
                .subtract(ee.Image(neighLandsat_t01.get(i)).select([ee.String(name).cat('_(.+)')]))
                .abs()
                .lte(ee.Image(thresh.get(i)).select([ee.String(name).cat('_(.+)')])); 
            });
          });
          return masks;
        }

        var mask_t01 = threshMask(neighLandsat_t01, thresh, commonBandNames);
        var maskArr_t01 = ee.ImageCollection(mask_t01.flatten()).toBands().toArray();

        function maskArrays_f(i){return maskArr_t01.arraySlice(0, ee.Number(i).int(),
                                numPixels.multiply(commonBandNames.length().multiply(2)).int(), numPixels);}

        var maskArrays = ee.List.sequence(0, numPixels.subtract(1)).map(maskArrays_f);
        var masksSorted = ee.List.sequence(0, numPixels.subtract(1)).map(function(i){ return ee.Image(maskArrays.get(i)).arrayFlatten([pixelBandNames.get(i)]);});

        var maskedLandsat = ee.List.sequence(0, numPixels.subtract(1)).map(function(index){
          return ee.Image(lanSorted.get(index)).updateMask(ee.Image(masksSorted.get(index)));
        });
        return [maskedLandsat, pixPositions, pixelBandNames];  
    }    
    
    var prepLandsat_1 = prepLandsat(landsat_t01, kernel, numPixels, commonBandNames, doys, coverClasses);
    var maskedLandsat = prepLandsat_1[0];
    var pixPositions = prepLandsat_1[1];
    var pixBN = prepLandsat_1[2];

    // Function to prepare MODIS data for further processing, similar to Landsat
    function prepMODIS(modis_t01, modis_tp, kernel, numPixels, commonBandNames, pixelBandNames){
      var neighMod_t01 = modis_t01.map(function(image){return ee.Image(image).neighborhoodToBands(kernel);});
      var neighMod_tp = modis_tp.map(function(image){return ee.Image(image).neighborhoodToBands(kernel);});

      var modArr = ee.Image(neighMod_t01.get(0)).toArray()
        .arrayCat(ee.Image(neighMod_t01.get(1)).toArray(), 0);

      var modPixArrays_t01 = ee.List.sequence(0, numPixels.subtract(1)).map(function(i){ 
      return modArr.arraySlice(0, ee.Number(i).int(),
                              numPixels.multiply(commonBandNames.length().multiply(2)).int(),
                              numPixels);
      });

      var modPixArrays_tp = neighMod_tp.map(function(image){ 
        return ee.List.sequence(0, numPixels.subtract(1)).map(function(i){
          return ee.Image(image).toArray().arraySlice(0, ee.Number(i).int(),
                                  numPixels.multiply(commonBandNames.length()).int(),
                                  numPixels);
        });
      });

      var modSorted_t01 = ee.List.sequence(0, numPixels.subtract(1))
        .map(function(i){ return ee.Image(modPixArrays_t01.get(i)).arrayFlatten([pixelBandNames.get(i)]);});

      var modSorted_tp = ee.List.sequence(0, modPixArrays_tp.length().subtract(1))
        .map(function(i){ return ee.List.sequence(0, numPixels.subtract(1))
                    .map(function (x){ return ee.Image(ee.List(modPixArrays_tp.get(i)).get(x))
                                              .arrayFlatten([commonBandNames])
                                              .set('DOY', ee.Image(modis_tp.get(i)).get('DOY'));});
      });
      return [modSorted_t01, modSorted_tp];
    }

    var prepMODIS_1 = prepMODIS(modis_t01, modis_tp, kernel, numPixels, commonBandNames, pixBN);
    var modSorted_t01 = prepMODIS_1[0];
    var modSorted_tp = prepMODIS_1[1];

    // Function to calculate the spectral distance between Landsat and MODIS pixels
    function calcSpecDist(maskedLandsat, modSorted_t01, numPixels, pixPositions){
        var sDist = ee.List.sequence(0, numPixels.subtract(1)).map(function(index){ 
          return ee.Image(maskedLandsat.get(index)).subtract(ee.Image(modSorted_t01.get(index)))
                          .abs().reduce(ee.Reducer.mean());
        });
        return ee.ImageCollection(sDist)
              .toBands()
              .rename(pixPositions.map(function(name){ return ee.String('sDist').cat(name);}));
    }

    var specDist = calcSpecDist(maskedLandsat, modSorted_t01, numPixels, pixPositions);      

    // Function to calculate the spatial distance between pixels
    function calcSpatDist(positions){
        var w2 = positions.length().sqrt().subtract(1).divide(2);

        var dist = positions.map(function(position){
            return ee.Image.constant(ee.Number(1).add(
                                      ee.Number.parse(ee.String(position)
                                      .match('(-?[0-9]+)', 'g')
                                      .get(0)).pow(2).add(ee.Number.parse(
                                      ee.String(position).match('(-?[0-9]+)', 'g')
                                                          .get(1)).pow(2)).sqrt().divide(w2)));
        });

        return ee.ImageCollection(dist).toBands()
              .rename(positions.map(function(bn){return ee.String('corr').cat(bn);}));
    }

    var spatDist = calcSpatDist(pixPositions);  

    // Function to calculate weights based on spectral and spatial distances
    function calcWeight(spatDist, specDist){
      var disIndex = specDist.multiply(spatDist);
      var num = ee.Image.constant(1).divide(disIndex);
      var sumNum = ee.ImageCollection(num.bandNames().map(function(bn){
                                  return num.select([bn]).rename('num');})).sum();
      var W = num.divide(sumNum);
      return W.unmask().toArray().toArray(1).matrixToDiag();
    }

    var weights = calcWeight(spatDist, specDist);

    // Function to calculate conversion coefficients between Landsat and MODIS
    function calcConversionCoeff(maskedLandsat, modSorted_t01, doys, numPixels, commonBandNames){
      var lanMod = doys.map(function(doy){ 
                      return ee.List.sequence(0, numPixels.subtract(1))
                      .map(function(index){
                      return ee.Image.constant(1).rename(['intercept'])
                                            .addBands(ee.Image(modSorted_t01.get(index))
                                            .select(ee.String(doy).cat('.+'))
                                            .rename(commonBandNames
                                            .map(function(bn){ return ee.String(bn).cat('_modis');})))
                                            .addBands(ee.Image(maskedLandsat.get(index))
                                                      .select(ee.String(doy).cat('.+'))
                                                      .rename(commonBandNames.map(function (bn){ return ee.String(bn).cat('_landsat');})));
                      });
      });

      var lanMod = ee.ImageCollection(lanMod.flatten());
      var coeffs = lanMod.reduce(ee.Reducer.linearRegression(commonBandNames.length().add(1),
                                            commonBandNames.length()))
        .select([0], ['coefficients'])
        .arraySlice(0, 1, commonBandNames.length().add(1));
      return coeffs;
    }

    var coeffs = calcConversionCoeff(maskedLandsat, modSorted_t01, doys, numPixels, commonBandNames);

    // Function to predict Landsat data from MODIS using calculated coefficients and weights
    function predictLandsat(landsat_t01, modSorted_t01, doys, modSorted_tp, weights, coeffs, commonBandNames, numPixels){
      var modSplit_t01 = doys.map(function(doy){
          return modSorted_t01.map(function(image){
            return ee.Image(image).select(ee.String(doy).cat('.+')).rename(commonBandNames);
          });
      });

      var diffMod = modSplit_t01.map(function(imagelist){
        return ee.ImageCollection(ee.List.sequence(0, numPixels.subtract(1))
                  .map(function(index){
                    return ee.Image(modSorted_tp.get(index))
                                          .subtract(ee.Image(ee.List(imagelist).get(index)));
                  }))
                  .toArray();
      });

      var sumPixel = diffMod.map(function(arrayImage){
        return weights.matrixMultiply(ee.Image(arrayImage).matrixMultiply(coeffs))
                  .arrayReduce(ee.Reducer.sum(), [0])
                  .arrayProject([1])
                  .arrayFlatten([commonBandNames]);
      });

      var sumDiff_t0 = ee.Image.constant(1).divide(ee.Image(diffMod.get(0))
                .arrayReduce(ee.Reducer.sum(), [0]).abs());
      var sumDiff_t1 = ee.Image.constant(1).divide(ee.Image(diffMod.get(1))
                .arrayReduce(ee.Reducer.sum(), [0]).abs());
      var weight_t0 = sumDiff_t0.divide(sumDiff_t0.add(sumDiff_t1)) 
                .arrayProject([1]).arrayFlatten([commonBandNames]);
      var weight_t1 = sumDiff_t1.divide(sumDiff_t0.add(sumDiff_t1)).arrayProject([1])
                .arrayFlatten([commonBandNames]);
      var temporalWeight = ee.List([weight_t0, weight_t1]);

      var predictions = ee.List([0, 1]).map(function(time){
        return ee.Image(landsat_t01.get(time))
                  .add(ee.Image(sumPixel.get(time)))
                  .multiply(ee.Image(temporalWeight.get(time)));
      });

      var mergedPrediction = ee.Image(predictions.get(0)).add(ee.Image(predictions.get(1)));

      return mergedPrediction.setMulti({
            'DOY': ee.Image(modSorted_tp.get(0)).get('DOY'),
            'system:time_start': ee.Image(modSorted_tp.get(0)).get('system:time_start')
      });
    }

    var prediction = modSorted_tp.map(function(image){
      return predictLandsat(landsat_t01, modSorted_t01, doys, ee.List(image), weights, coeffs, commonBandNames, numPixels);
    });

    var preds = ee.ImageCollection(prediction).toBands();
    var dates = modis_tp.map(function(img){ return ee.Image(img).get('system:time_start');});
     
    var predNames = ee.List.sequence(0, prediction.length().subtract(1)).map(function(i){
      return commonBandNames.map(function(name){
        return ee.String(name).cat(ee.String(ee.Number(dates.get(i)).format()));
      });
    }).flatten();

    var preds = preds.rename(predNames);
    return preds;
  });
  return list_index_seq;
}

// Generate results for all sub-collections
var results = num_lists.map(tonumlist);
var results_list = results.flatten();
var results_tobands = ee.ImageCollection(results_list).toBands();
var bandname = results_tobands.bandNames();

// Rename bands for consistent output
var bandname0 = bandname.map(function(i){
  var namei0 = ee.String(i).slice(-17);
  return namei0;
});
results_tobands = results_tobands.rename(bandname0);

// Export the final result as an asset in Google Earth Engine
Export.image.toAsset({
  image:results_tobands,
  description: 'ZJK_step3_1',
  assetId:'ZJK_step3_1',
  region:ROI,
  scale:30,
  crs:"EPSG:4326",
  maxPixels:1e13
});