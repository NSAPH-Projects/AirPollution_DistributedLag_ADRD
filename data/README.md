# Data Dictionary

## pm25 exposure 

Location: `/nfs/home/D/dam9096/shared_space/ci3_exposure/pm25/whole_us/annual/zipcode/qd_predictions_ensemble/ywei_aggregation/`

```
"ZIP","pm25","year"
"00001",6.21,2000
"00002",8.67,2000
"00003",5.98,2000
"00004",9.30,2000
```

## pm25_components

Location: `/nfs/home/D/dam9096/shared_space/ci3_exposure/pm25_components/`

```
"ZIP","br","ca","cu","ec","fe","k","nh4","ni","no3","oc","pb","si","so4","v","z"
"00012",1.13,29.478,0.322,0.197,33.10,32.14,0.18,0.015,0.21,1.046,0.60,106.32,0.53,0.69,1.46
"00014",2.46,122.219,1.031,0.189,88.15,77.50,0.24,0.07,0.42,0.87,1.34,309.78,1.1,1.09,2.87
```

For ec, oc, nh4, no3, and so4 the units are microgram per cubic meter; for br, ca, cu, fe, k, ni, pb, si, v, and z the units are nanogram per cubic meter.

## no2 exposure 

Location: `/nfs/home/D/dam9096/shared_space/ci3_exposure/no2/whole_us/annual/zipcode/qd_predictions_ensemble/ywei_aggregations/`

```
"ZIP","no2","year"
"00001",6.006,2000
"00002",4.547,2000
"00003",13.28,2000
"00004",3.420,2000
```

## ozone exposure

Location: `/nfs/home/D/dam9096/shared_space/ci3_exposure/ozone/whole_us/seasonal/zipcode/`

```
ZIP,year,ozone_summer,ozone_winter,ozone_fall,ozone_spring
1,2000,41.77,36.11,35.88,42.11
1,2001,43.53,38.53,37.727,46.83
1,2002,45.22,38.20,38.863,46.04
1,2003,45.45,36.12,39.23,46.17
```

## temperature

Location: `/nfs/home/D/dam9096/shared_space/ci3_confounders/data_for_analysis/prepped_temperature/annual/`

| Column | Description |
|--------|-------------|
| ZIP  | zip code |
| tmmx | Maxmimum Daily Temperture in degrees Kelvin, annual average of daily data |
| rmax | Maximum Daily Relative Humidity (0-100%), annual average of daily data |
| pr   | Daily Precipitation (mm), annual average of daily data |
| year | year |

```
ZIP,year,tmmx,rmax,pr
4742,1999,284.41,85.24,2.79
4740,1999,284.58,84.90,2.89
4758,1999,284.80,85.58,3.00
4769,1999,284.74,84.80,2.81
4787,1999,284.65,85.57,2.99
4735,1999,285.13,86.28,3.04
4760,1999,285.30,86.69,3.05
4734,1999,285.17,86.08,3.03
```

## denominator 

Location: `/nfs/home/D/dam9096/shared_space/ci3_health_data/medicare/mortality/1999_2016/wu/cache_data/merged_by_year_v2/`

| Column | Description |
|--------|-------------|
| qid    | ID          |
| year   |             |
| zip    |             |
| hmo_mo |             |
| age    |             |
| dead   |             |
