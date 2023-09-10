################################################################################
###
### Code illustrating the analyses in:
###
### Heat-related mortality in Europe during the summer of 2022
### Nature Medicine 29, 1857–1866 (2023)
###
### https://www.nature.com/articles/s41591-023-02419-z
### https://doi.org/10.1038/s41591-023-02419-z
###
### Contact: Joan Ballester (joan.ballester@isglobal.org)
###
### Note: This is a simplified code with sample data, results are expected to
###       differ from those published in the article
### 
################################################################################

rm( list = ls() );
cat("\014");

# Required Libraries
suppressMessages( library(lubridate) ); # wday
suppressMessages( library(ISOweek) ); # ISOweek2date
suppressMessages( library(dlnm) ); # crossbasis
suppressMessages( library(splines) ); # ns, bs
suppressMessages( library(mixmeta) ); # mixmeta
suppressMessages( library(tsModel) ); # Lag
suppressMessages( library(MASS) ); # mvrnorm



################################################################################
### Parameter Definition
################################################################################

print("");
print("= Parameter Definition =");
print("");

# Period for the Calibration: 2015-2019
DATE1_CALI = as.Date( "2015-01-01" );
DATE2_CALI = as.Date( "2019-12-26" );

# Period for the Predictions: 2015-2022
DATE1_PRED = as.Date( "2015-01-01" );
DATE2_PRED = as.Date( "2022-11-03" );

# Period for Summer 2022
DATE1_SU22 = as.Date( "2022-06-02" );
DATE2_SU22 = as.Date( "2022-09-01" );

# For Convenience, We Identify the Date of Each Week with Its Thursday
if( wday( DATE1_CALI, week_start = 1 ) != 4 | wday( DATE2_CALI, week_start = 1 ) != 4 ){ stop("ERROR: Invalid Calibration Period !!!"); }
if( wday( DATE1_PRED, week_start = 1 ) != 4 | wday( DATE2_PRED, week_start = 1 ) != 4 ){ stop("ERROR: Invalid Prediction Period !!!"); }
if( wday( DATE1_SU22, week_start = 1 ) != 4 | wday( DATE2_SU22, week_start = 1 ) != 4 ){ stop("ERROR: Invalid Summer 2022 Period !!!"); }

# Exposure-Response: Spline Type and Degree
VAR_FUN = "ns"; VAR_DEG = NA;
# VAR_FUN = "bs"; VAR_DEG = 2;

# Exposure-Response: Temperature Knot Percentiles
VAR_PRC = c(10,50,90) / 100;

# Lag-Response: Minimum and Maximum Lags (in Weeks)
MIN_LAG = 0; if( MIN_LAG <     0    ){ stop( "ERROR: Invalid MIN_LAG !!!" ); }
MAX_LAG = 3; if( MAX_LAG <= MIN_LAG ){ stop( "ERROR: Invalid MAX_LAG !!!" ); }

# Degrees of Freedom per Year for the Seasonal and Long-Term Trends
DF_SEAS = 8; if( DF_SEAS <= 0 ){ stop( "ERROR: Invalid DF_SEAS !!!" ); }

# Cumulative Exposure-Response: Predicted Temperature Percentiles
PRED_PRC = sort( unique( c( seq(  0.0,   1.0, 0.1 ),
                            seq(  1.5,   5.0, 0.5 ),
                            seq(  6.0,  94.0, 1.0 ),
                            seq( 95.0,  98.5, 0.5 ),
                            seq( 99.0, 100.0, 0.1 ) ) / 100 ) );
if( any( 0 > PRED_PRC | PRED_PRC > 1 ) ){ stop("ERROR: Invalid Predicted Temperature Percentiles !!!"); }

# Temperature Percentile Range for the Minimum Mortality Temperature
MIN_PMMT =   5 / 100; if( any( 0 > MIN_PMMT | MIN_PMMT > 1 ) ){ stop("ERROR: Invalid Lower Temperature Percentile Range for the Minimum Mortality Temperature !!!"); }
MAX_PMMT = 100 / 100; if( any( 0 > MAX_PMMT | MAX_PMMT > 1 ) ){ stop("ERROR: Invalid Upper Temperature Percentile Range for the Minimum Mortality Temperature !!!"); }



################################################################################
### Data Preparation
################################################################################

print("");
print("= Data Preparation =");
print("");

# Reading the Datatable
DATATABLE = read.csv( "./data.csv" );

# For Convenience, We Identify the Date of Each Week with Its Thursday
DATATABLE$date = ISOweek2date( paste0( DATATABLE$year, "-W", sprintf( "%02d", DATATABLE$woy ), "-", 4 ) );

# Reading the Metatable
METATABLE = read.csv( "./metadata.csv" );

# Creating the Vector of NUTS Region Codes
vREG = METATABLE$location;
nREG = length(vREG);

# Restricting the Datatable to the Calibration and Prediction Periods
# WARNING: The Datatables Include MAX_LAG Additional Weeks to Calculate All the Attributable Mortality Numbers
DATATABLE_CALI = DATATABLE[ which( DATE1_CALI <= DATATABLE$date & DATATABLE$date <= DATE2_CALI + 7 * MAX_LAG ), ];
DATATABLE_PRED = DATATABLE[ which( DATE1_PRED <= DATATABLE$date & DATATABLE$date <= DATE2_PRED + 7 * MAX_LAG ), ];
rm(DATATABLE);

# Creating the Datalists of the Calibration and Prediction Periods
DATALIST_CALI = lapply( vREG, function(x) DATATABLE_CALI[ DATATABLE_CALI$location == x, ] ); names(DATALIST_CALI) = vREG;
DATALIST_PRED = lapply( vREG, function(x) DATATABLE_PRED[ DATATABLE_PRED$location == x, ] ); names(DATALIST_PRED) = vREG;
rm(DATATABLE_CALI,DATATABLE_PRED, METATABLE);

# Adding the Variable Week-of-Period to the Calibration Period
for( iREG in 1:nREG ){
  DATALIST_CALI[[iREG]]$wop = 1:length(DATALIST_CALI[[iREG]]$woy);
}
rm(iREG);



################################################################################
### Location-Specific Associations
################################################################################

print("");
print("= Location-Specific Associations =");
print("");

# Reduced Coefficients
if      ( VAR_FUN == "ns" ){ COEF_MODEL = matrix( data = NA, nREG, length(VAR_PRC) +    1   , dimnames = list(vREG) );
}else if( VAR_FUN == "bs" ){ COEF_MODEL = matrix( data = NA, nREG, length(VAR_PRC) + VAR_DEG, dimnames = list(vREG) );
}else                      { stop("ERROR: Invalid VAR_FUN !!!"); }

# Reduced Covariance Matrices
VCOV_MODEL = vector( "list", nREG );
names(VCOV_MODEL) = vREG;

for( iREG in 1:nREG ){
  
  # Cross-Basis of Temperature
  if( VAR_FUN == "ns" ){
    CROSS_BASIS = crossbasis( DATALIST_CALI[[iREG]]$temp,
                              c( MIN_LAG, MAX_LAG ),
                              argvar = list( fun = VAR_FUN,
                                             knots = quantile( DATALIST_CALI[[iREG]]$temp, VAR_PRC, na.rm = TRUE ),
                                             Boundary.knots = range( DATALIST_CALI[[iREG]]$temp, na.rm = TRUE ) ),
                              arglag = list( fun = "integer" ) );
  }else if( VAR_FUN == "bs" ){
    CROSS_BASIS = crossbasis( DATALIST_CALI[[iREG]]$temp,
                              c( MIN_LAG, MAX_LAG ),
                              argvar = list( fun = VAR_FUN,
                                             degree = VAR_DEG,
                                             knots = quantile( DATALIST_CALI[[iREG]]$temp, VAR_PRC, na.rm = TRUE ) ),
                              arglag = list( fun = "integer" ) );
  }else{
    stop("ERROR: Invalid VAR_FUN !!!");
  }
  
  # Fitting of the Model
  GLM_MODEL = glm( formula = mort ~ ns( wop, df = round( DF_SEAS * length(wop) * 7 / 365.25 ) ) + CROSS_BASIS,
                   DATALIST_CALI[[iREG]],
                   family = quasipoisson,
                   na.action = "na.exclude" );
  
  # Cumulative Exposure-Response without Centring (Before the Meta-Analysis)
  suppressMessages( CROSS_PRED_REG_NOMETA <- crosspred( CROSS_BASIS,
                                                        GLM_MODEL,
                                                        at = quantile( DATALIST_CALI[[iREG]]$temp, PRED_PRC, na.rm = TRUE ) ) );
  
  # The Minimum Mortality Temperature Is the Lowest Relative Risk Value within the Predefined Temperature Percentile Range [MIN_PMMT,MAX_PMMT]
  MMT = CROSS_PRED_REG_NOMETA$predvar[ which( PRED_PRC == MIN_PMMT ) - 1 +
                                       which.min( CROSS_PRED_REG_NOMETA$allRRfit[ which( PRED_PRC == MIN_PMMT ) : which( PRED_PRC == MAX_PMMT ) ] ) ];
  rm(CROSS_PRED_REG_NOMETA);
  
  # Reduced Coefficients and Covariance Matrices
  REDUCED = crossreduce( CROSS_BASIS, GLM_MODEL, cen = MMT );
  COEF_MODEL[iREG,] = coef( REDUCED );
  VCOV_MODEL[[iREG]] = vcov( REDUCED );
  rm(REDUCED);
  
  rm(CROSS_BASIS,GLM_MODEL, MMT);
  
}
rm(iREG);



################################################################################
### Best Linear Unbiased Predictions
################################################################################

print("");
print("= Best Linear Unbiased Predictions =");
print("");

# Regional Meta-Predictors: Temperature Average (TEMP_AVG), Temperature Inter-Quartile Range (TEMP_IQR)
TEMP_AVG = sapply( DATALIST_CALI, function(x) mean( x$temp, na.rm = TRUE ) );
TEMP_IQR = sapply( DATALIST_CALI, function(x)  IQR( x$temp, na.rm = TRUE ) );

# Regional Minimum Mortality Temperature
MMT_REG = array( NA, dim = c( nREG ), dimnames = list( vREG ) );

# Cumulative Exposure-Response (After the Meta-Analysis)
CROSS_PRED_REG_META = vector( "list", nREG );
names(CROSS_PRED_REG_META) = vREG;

# Meta-Analysis of the Reduced Coefficients
# WARNING: For Simplicity, this Includes Two Meta-Predictors without Country Random Effects
#          The Article Includes Three Meta-Predictors with Country Random Effects
MULTIVAR = mixmeta( COEF_MODEL ~ TEMP_AVG + TEMP_IQR,
                    VCOV_MODEL,
                    data = data.frame( vREG = vREG ),
                    control = list( showiter = TRUE, igls.inititer = 10 ),
                    method = "reml" );
rm(COEF_MODEL,VCOV_MODEL, TEMP_AVG,TEMP_IQR);
print("");
print( summary( MULTIVAR ) );

# Best Linear Unbiased Predictions
BLUP = blup( MULTIVAR, vcov = TRUE );
rm(MULTIVAR);

for( iREG in 1:nREG ){
  
  # Temperature Basis
  if( VAR_FUN == "ns" ){
    BASIS_VAR = onebasis( quantile( DATALIST_CALI[[iREG]]$temp, PRED_PRC, na.rm = TRUE ),
                          fun = VAR_FUN,
                          knots = quantile( DATALIST_CALI[[iREG]]$temp, VAR_PRC, na.rm = TRUE ),
                          Boundary.knots = range( DATALIST_CALI[[iREG]]$temp, na.rm = TRUE ) );
  }else if( VAR_FUN == "bs" ){
    BASIS_VAR = onebasis( quantile( DATALIST_CALI[[iREG]]$temp, PRED_PRC, na.rm = TRUE ),
                          fun = VAR_FUN,
                          degree = VAR_DEG,
                          knots = quantile( DATALIST_CALI[[iREG]]$temp, VAR_PRC, na.rm = TRUE ) );
  }else{
    stop("ERROR: Invalid Temperature Basis !!!");
  }
  
  # Vector of Percentile Temperatures for the Predictions
  PRC_VAR = quantile( DATALIST_CALI[[iREG]]$temp, PRED_PRC, na.rm = TRUE );
  
  # Cumulative Exposure-Response without Centring (After the Meta-Analysis)
  suppressMessages( PRED_MORT <- crosspred( BASIS_VAR,
                                            coef = BLUP[[iREG]]$blup,
                                            vcov = BLUP[[iREG]]$vcov,
                                            model.link = "log",
                                            at = PRC_VAR ) );
  
  # The Minimum Mortality Temperature Is the Lowest Relative Risk Value within the Predefined Temperature Percentile Range [MIN_PMMT,MAX_PMMT]
  MMT_REG[iREG] = PRC_VAR[ which( PRED_PRC == MIN_PMMT ) - 1 +
                           which.min( PRED_MORT$allRRfit[ which( PRED_PRC == MIN_PMMT ) : which( PRED_PRC == MAX_PMMT ) ] ) ];
  rm(PRED_MORT);
  
  # Cumulative Exposure-Response (After the Meta-Analysis)
  CROSS_PRED_REG_META[[iREG]] = crosspred( BASIS_VAR,
                                           coef = BLUP[[iREG]]$blup,
                                           vcov = BLUP[[iREG]]$vcov,
                                           model.link = "log",
                                           at = PRC_VAR,
                                           cen = MMT_REG[iREG] );
  rm(BASIS_VAR,PRC_VAR);
  
}
rm(iREG);



################################################################################
### Attributable Mortality
################################################################################

print("");
print("= Attributable Mortality =");
print("");

# Vector of Predicted Time Periods: Individual Weeks
# WARNING: We Exclude the MAX_LAG Weeks that Were Initially Added to Calculate All the Attributable Mortality Numbers
vPER_PRED = as.character( DATALIST_PRED[[1]]$date[ 1 : ( dim(DATALIST_PRED[[1]])[1] - MAX_LAG ) ] );
nPER_PRED = length(vPER_PRED);

# Vector of Predicted Time Periods: Individual Weeks, Summer 2022
# WARNING: Other Multi-Week Periods Can Be Added Here
vPER_PRED = c( vPER_PRED, "Summer 2022" );

# Vector of Temperature Ranges
vRNG = c( "Total", "Total Cold", "Total Heat", "Moderate Cold", "Moderate Heat", "Extreme Cold", "Extreme Heat" );
nRNG = length(vRNG);

# Number of Simulations
nSIM = 10^3;

# Values and Simulations of the Regional Attributable Numbers
ATT_VAL = array( NA, dim = c( nREG, length(vPER_PRED), nRNG       ), dimnames = list( vREG, vPER_PRED, vRNG         ) );
ATT_SIM = array( NA, dim = c( nREG, length(vPER_PRED), nRNG, nSIM ), dimnames = list( vREG, vPER_PRED, vRNG, 1:nSIM ) );

# Attributable Number (AN, in Deaths), Fraction (AF, in %) and Incidence (AI, in Deaths per Million) for Regions (REG) and Overall (TOT)
AN_REG = AF_REG = AI_REG = array( NA, dim = c( nREG, length(vPER_PRED), nRNG, 3 ), dimnames = list(    vREG  , vPER_PRED, vRNG, c( "att_val", "att_low", "att_upp" ) ) );
AN_TOT = AF_TOT = AI_TOT = array( NA, dim = c(   1 , length(vPER_PRED), nRNG, 3 ), dimnames = list( "Overall", vPER_PRED, vRNG, c( "att_val", "att_low", "att_upp" ) ) );

# Mortality (MORT) and Population (POPU) for Regions (REG) and Overall (TOT)
MORT_REG = POPU_REG = array( NA, dim = c( nREG, length(vPER_PRED) ), dimnames = list(    vREG  , vPER_PRED ) );
MORT_TOT = POPU_TOT = array( NA, dim = c(   1 , length(vPER_PRED) ), dimnames = list( "Overall", vPER_PRED ) );

# Regional Thresholds for the Temperature Ranges
TMIN = P025 = P975 = TMAX = array( NA, dim = nREG, dimnames = list(vREG) );
for( iREG in 1:nREG ){
  TMIN[iREG] =      min( DATALIST_PRED[[iREG]]$temp,        na.rm = TRUE );
  P025[iREG] = quantile( DATALIST_PRED[[iREG]]$temp, 0.025, na.rm = TRUE );
  P975[iREG] = quantile( DATALIST_PRED[[iREG]]$temp, 0.975, na.rm = TRUE );
  TMAX[iREG] =      max( DATALIST_PRED[[iREG]]$temp,        na.rm = TRUE );
}
rm(iREG);

### Loop for Regions ###
for( iREG in 1:nREG ){
  print( paste0( "   Region ", iREG, ": ", vREG[iREG] ) );
  
  # Temperature Basis
  if( VAR_FUN == "ns" ){
    BASIS_VAR = onebasis( DATALIST_PRED[[iREG]]$temp,
                          fun = VAR_FUN,
                          knots = quantile( DATALIST_CALI[[iREG]]$temp, VAR_PRC, na.rm = TRUE ),
                          Boundary.knots = range( DATALIST_CALI[[iREG]]$temp, na.rm = TRUE ) );
    
  }else if( VAR_FUN == "bs" ){
    BASIS_VAR = onebasis( DATALIST_PRED[[iREG]]$temp,
                          fun = VAR_FUN,
                          degree = VAR_DEG,
                          knots = quantile( DATALIST_CALI[[iREG]]$temp, VAR_PRC, na.rm = TRUE ) );
  }else{
    stop("ERROR: Invalid Temperature Basis !!!");
  }
  
  # Minimum Mortality Temperature for the Centring of the Temperature Basis
  if( VAR_FUN == "ns" ){
    BASIS_MMT = onebasis( MMT_REG[iREG],
                          fun = VAR_FUN,
                          knots = quantile( DATALIST_CALI[[iREG]]$temp, VAR_PRC, na.rm = TRUE ),
                          Boundary.knots = range( DATALIST_CALI[[iREG]]$temp, na.rm = TRUE ) );
  }else if( VAR_FUN == "bs" ){
    BASIS_MMT = onebasis( MMT_REG[iREG],
                          fun = VAR_FUN,
                          degree = VAR_DEG,
                          knots = quantile( DATALIST_CALI[[iREG]]$temp, VAR_PRC, na.rm = TRUE ) );
  }else{
    stop("ERROR: Invalid Minimum Mortality Temperature for the Centring of the Temperature Basis !!!");
  }
  
  # Temperature Basis After Centring at the Minimum Mortality Temperature
  BASIS_CEN = scale( BASIS_VAR, center = BASIS_MMT, scale = FALSE );
  rm(BASIS_VAR,BASIS_MMT);
  
  # Matrix of Lagged Deaths
  LAGGED_MORT_MATRIX = Lag( DATALIST_PRED[[iREG]]$mort, -MIN_LAG:-MAX_LAG );
  
  # Average of Deaths Across Lags
  LAGGED_MORT_VECTOR = rowMeans( LAGGED_MORT_MATRIX );
  
  # Weekly Time Series Vector of Attributable Numbers (Dimensions: Time x 1)
  ATT_NUM_TS_REF = ( 1 - exp( -BASIS_CEN %*% BLUP[[iREG]]$blup ) ) * LAGGED_MORT_VECTOR;
  
  # Perturbed BLUPs (Dimensions: Coefficients x Simulations)
  set.seed( 5634654 );
  COEF_SIM = t( mvrnorm( nSIM, BLUP[[iREG]]$blup, BLUP[[iREG]]$vcov ) );
  
  # Weekly Time Series Matrix of Attributable Numbers (Dimensions: Time x Simulations)
  ATT_NUM_TS_SIM = ( 1 - exp( -BASIS_CEN %*% COEF_SIM ) ) * LAGGED_MORT_VECTOR;
  rm(COEF_SIM);
  
  # Loop for Predicted Time Periods
  for( iPER in 1:length(vPER_PRED) ){
    
    # Selection of Time Period Data
    if     (           iPER  <=   nPER_PRED   ){ vTIM = which(               DATALIST_PRED[[iREG]]$date == as.Date( vPER_PRED[iPER] )              ); } # Individual Weeks
    else if( vPER_PRED[iPER] == "Summer 2022" ){ vTIM = which( DATE1_SU22 <= DATALIST_PRED[[iREG]]$date & DATALIST_PRED[[iREG]]$date <= DATE2_SU22 ); } # Summer 2022
    else                                       { stop("ERROR: Invalid Time Period Data !!!"); }
    if( length(vTIM) <= 0 ){ stop("ERROR: Invalid vTIM !!!"); }
    
    # Mortality and Population for the Region and Predicted Time Period
    MORT_REG[iREG,iPER] =  sum( DATALIST_PRED[[iREG]]$mort[vTIM], na.rm = TRUE );
    POPU_REG[iREG,iPER] = mean( DATALIST_PRED[[iREG]]$popu[vTIM], na.rm = TRUE );
    
    # Loop for Temperature Ranges
    for( iRNG in 1:nRNG ){
      
      # Temperature Range Data
      if     ( vRNG[iRNG] == "Total"         ){ vRNG_THRES =                                                    1:length(vTIM)                                             ; } # Total
      else if( vRNG[iRNG] == "Total Cold"    ){ vRNG_THRES = which(                                                     DATALIST_PRED[[iREG]]$temp[vTIM] <  MMT_REG[iREG] ); } # Total Cold
      else if( vRNG[iRNG] == "Total Heat"    ){ vRNG_THRES = which( MMT_REG[iREG] <  DATALIST_PRED[[iREG]]$temp[vTIM]                                                     ); } # Total Heat
      else if( vRNG[iRNG] == "Moderate Cold" ){ vRNG_THRES = which(    P025[iREG] <= DATALIST_PRED[[iREG]]$temp[vTIM] & DATALIST_PRED[[iREG]]$temp[vTIM] <  MMT_REG[iREG] ); } # Moderate Cold
      else if( vRNG[iRNG] == "Moderate Heat" ){ vRNG_THRES = which( MMT_REG[iREG] <  DATALIST_PRED[[iREG]]$temp[vTIM] & DATALIST_PRED[[iREG]]$temp[vTIM] <=    P975[iREG] ); } # Moderate Heat
      else if( vRNG[iRNG] == "Extreme Cold"  ){ vRNG_THRES = which(                                                     DATALIST_PRED[[iREG]]$temp[vTIM] <     P025[iREG] ); } # Extreme Cold
      else if( vRNG[iRNG] == "Extreme Heat"  ){ vRNG_THRES = which(    P975[iREG] <  DATALIST_PRED[[iREG]]$temp[vTIM]                                                     ); } # Extreme Heat
      else                                    { stop("ERROR: Invalid Temperature Range !!!"); }
      if( length(vRNG_THRES) > 0 & sum( LAGGED_MORT_VECTOR[vTIM] ) > 0 ){
        
        # Attributable Number: Value
        ATT_VAL[iREG,iPER,iRNG ] = sum( ATT_NUM_TS_REF[ vTIM[vRNG_THRES] ], na.rm = TRUE ) *
                                   sum( rowMeans( LAGGED_MORT_MATRIX[ vTIM, , drop=FALSE ], na.rm = TRUE ), na.rm = TRUE ) /
                                   sum(           LAGGED_MORT_VECTOR[ vTIM               ]                , na.rm = TRUE );
        if(      !is.finite( ATT_VAL[iREG,iPER,iRNG ] )   ){ stop("ERROR: Invalid ATT_VAL !!!"); }
        
        # Attributable Number: Simulations
        ATT_SIM[iREG,iPER,iRNG,] = colSums( ATT_NUM_TS_SIM[ vTIM[vRNG_THRES], , drop=FALSE ], na.rm = TRUE ) *
                                   sum( rowMeans( LAGGED_MORT_MATRIX[ vTIM, , drop=FALSE ], na.rm = TRUE ), na.rm = TRUE ) /
                                   sum(           LAGGED_MORT_VECTOR[ vTIM               ]                , na.rm = TRUE );
        if( any( !is.finite( ATT_SIM[iREG,iPER,iRNG,] ) ) ){ stop("ERROR: Invalid ATT_SIM !!!"); }
        
        # Attributable Number: Value and Confidence Interval
        AN_REG[iREG,iPER,iRNG, 1 ] =           ATT_VAL[iREG,iPER,iRNG ]                  ;
        AN_REG[iREG,iPER,iRNG,2:3] = quantile( ATT_SIM[iREG,iPER,iRNG,], c(0.025,0.975) );
        
        # Attributable Fraction (in %)
        if( MORT_REG[iREG,iPER] > 0 ){ AF_REG[iREG,iPER,iRNG,] = 10^2 * AN_REG[iREG,iPER,iRNG,] / MORT_REG[iREG,iPER]; }
        
        # Attributable Incidence (in Deaths per Million)
        if( POPU_REG[iREG,iPER] > 0 ){ AI_REG[iREG,iPER,iRNG,] = 10^6 * AN_REG[iREG,iPER,iRNG,] / POPU_REG[iREG,iPER]; }
        else                         { stop("ERROR: Invalid Regional Population !!!"); }
        
      }else{
        ATT_VAL[iREG,iPER,iRNG ] = 0;
        ATT_SIM[iREG,iPER,iRNG,] = 0;
         AN_REG[iREG,iPER,iRNG,] = 0;
         AF_REG[iREG,iPER,iRNG,] = 0;
         AI_REG[iREG,iPER,iRNG,] = 0;
      }
      rm(vRNG_THRES);
      
    }
    rm(vTIM, iRNG);
    
  }
  rm(BASIS_CEN, LAGGED_MORT_MATRIX,LAGGED_MORT_VECTOR, ATT_NUM_TS_REF,ATT_NUM_TS_SIM, iPER);
  
}
rm(BLUP, TMIN,P025,P975,TMAX, iREG);

### Overall ###
print( "   Overall" );

# Mortality and Population of Overall and Predicted Time Period
MORT_TOT[1,] = colSums( MORT_REG );
POPU_TOT[1,] = colSums( POPU_REG );

# Attributable Number: Value and Confidence Interval
AN_TOT[1,,,1] =        apply( ATT_VAL, c(2,3  ), sum )                           ;
AN_TOT[1,,,2] = apply( apply( ATT_SIM, c(2,3,4), sum ), c(1,2), quantile, 0.025 );
AN_TOT[1,,,3] = apply( apply( ATT_SIM, c(2,3,4), sum ), c(1,2), quantile, 0.975 );

# Loop for Predicted Time Periods
for( iPER in 1:length(vPER_PRED) ){
  
  # Attributable Fraction (in %)
  if( MORT_TOT[1,iPER] > 0 ){ AF_TOT[1,iPER,,] = 10^2 * AN_TOT[1,iPER,,] / MORT_TOT[1,iPER]; }
  
  # Attributable Incidence (in Deaths per Million)
  if( POPU_TOT[1,iPER] > 0 ){ AI_TOT[1,iPER,,] = 10^6 * AN_TOT[1,iPER,,] / POPU_TOT[1,iPER]; }
  else                      { stop("ERROR: Invalid Overall Population !!!"); }
  
}
rm(iPER);

rm(ATT_VAL,ATT_SIM, MORT_REG,POPU_REG, MORT_TOT,POPU_TOT);

rm(DATE1_CALI,DATE2_CALI, DATE1_PRED,DATE2_PRED, DATE1_SU22,DATE2_SU22, VAR_FUN,VAR_DEG,VAR_PRC, MIN_LAG,MAX_LAG, DF_SEAS, PRED_PRC, MIN_PMMT,MAX_PMMT);
rm(DATALIST_CALI,DATALIST_PRED);
rm(vRNG,nRNG, nSIM);



################################################################################
### Plots
################################################################################

print("");
print("= Plots =");
print("");


# Heat Related Mortality for Summer 2022
# WARNING: Numbers Do Not Need to Coincide with Those in the Manuscript
print( "Heat Related Mortality for Summer 2022:" );
print( round( rbind( AN_REG[ , "Summer 2022", "Total Heat", ],
                     AN_TOT[ , "Summer 2022", "Total Heat", ] ), digits = 0 ) );

# Attributable Number

pdf( "./plot_AN.pdf", width = 8, height = 4 );
layout( matrix( seq(1*1), nrow = 1, byrow = TRUE ) );
par( mex = 0.8, mgp = c(2.5,1,0), las = 0 );

plot   (    as.Date( vPER_PRED[1:nPER_PRED] )                                            ,    AN_TOT[1,1:nPER_PRED,"Total Cold","att_val"]                                                       , col = rgb(0.00,0.00,1.00    ), lwd = 2, lty = 1, type = "l", main = "Temperature Related Mortality", xlab = "Time (week)", ylab = "Mortality (weekly deaths)", xaxt = "n" );
polygon( c( as.Date( vPER_PRED[1:nPER_PRED] ), rev( as.Date( vPER_PRED[1:nPER_PRED] ) ) ), c( AN_TOT[1,1:nPER_PRED,"Total Cold","att_low"], rev( AN_TOT[1,1:nPER_PRED,"Total Cold","att_upp"] ) ), col = rgb(0.00,0.00,1.00,1/3), border = FALSE );
polygon( c( as.Date( vPER_PRED[1:nPER_PRED] ), rev( as.Date( vPER_PRED[1:nPER_PRED] ) ) ), c( AN_TOT[1,1:nPER_PRED,"Total Heat","att_low"], rev( AN_TOT[1,1:nPER_PRED,"Total Heat","att_upp"] ) ), col = rgb(1.00,0.00,0.00,1/3), border = FALSE );
lines  (    as.Date( vPER_PRED[1:nPER_PRED] )                                            ,    AN_TOT[1,1:nPER_PRED,"Total Cold","att_val"]                                                       , col = rgb(0.00,0.00,1.00    ), lwd = 2, lty = 1 );
lines  (    as.Date( vPER_PRED[1:nPER_PRED] )                                            ,    AN_TOT[1,1:nPER_PRED,"Total Heat","att_val"]                                                       , col = rgb(1.00,0.00,0.00    ), lwd = 2, lty = 1 );
abline ( col = rgb(0.00,0.00,0.00), lwd = 1, lty = 1, h = 0 );
legend ( "top", c( "Total Cold", "Total Heat" ), col = rgb( c(0.00,1.00), c(0.00,0.00), c(1.00,0.00) ), lwd = 2, lty = 1, box.lty = 0, horiz = TRUE, bg = "transparent" );
axis   ( 1, at = as.Date( c( "2015-07-01", "2016-07-01", "2017-07-01", "2018-07-01", "2019-07-01", "2020-07-01", "2021-07-01", "2022-07-01" ) ),
            labels = c( "2015", "2016", "2017", "2018", "2019", "2020", "2021", "2022" ) );

dev.off();


# Cumulative Exposure-Response

pdf( "./plot_RR.pdf", width = 18, height = 9 );
layout( matrix( seq(6*3), nrow = 3, byrow = TRUE ) );
par( mex = 0.8, mgp = c(2.5,1,0), las = 0 );

for( iREG in 1:nREG ){
  plot   (    CROSS_PRED_REG_META[[iREG]]$predvar                                              ,    CROSS_PRED_REG_META[[iREG]]$allRRfit                                                , col = rgb(0.00,0.00,0.00    ), lwd = 4, lty = 1, type = "l", main = vREG[iREG], xlab = expression( paste( "Temperature (", degree, "C)" ) ), ylab = "Relative Risk", axes = T );
  polygon( c( CROSS_PRED_REG_META[[iREG]]$predvar, rev( CROSS_PRED_REG_META[[iREG]]$predvar ) ), c( CROSS_PRED_REG_META[[iREG]]$allRRlow, rev( CROSS_PRED_REG_META[[iREG]]$allRRhigh ) ), col = rgb(0.00,0.00,0.00,1/3), border = FALSE );
  lines  (    CROSS_PRED_REG_META[[iREG]]$predvar                                              ,    CROSS_PRED_REG_META[[iREG]]$allRRfit                                                , col = rgb(0.00,0.00,0.00    ), lwd = 4, lty = 1 );
  abline( v = MMT_REG[iREG], lwd = 1, lty = 1, col = rgb(0.00,0.00,0.00) );
  abline( h =       1      , lwd = 1, lty = 1, col = rgb(0.00,0.00,0.00) );
}
rm(iREG);

dev.off();


# "We are a way for the cosmos to know itself" (Carl Sagan)
# "We've tried nothing, and we're all out of ideas" (The Simpsons)
