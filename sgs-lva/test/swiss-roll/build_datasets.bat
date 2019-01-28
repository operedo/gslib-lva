:: first run octave script run.m
:: second, build lva field from swiss roll
perl script.pl swiss_roll_v4.csv > lva_field_swiww_roll_v4.csv
:: third, build gridded version of swiss roll
perl makeGrid.pl swiss_roll_v4.csv > swiss_roll_v4_gridded.csv
:: fourth, build lva from gridded swiss roll
perl script.pl swiss_roll_v4_gridded.csv > lva_field_swiss_roll_v4_gridded.csv
:: fifth, build the drillholes for testing
perl makeDrillholes.pl swiss_roll_v4_labels.csv > swiss_roll_v4_gridded_drillholes.csv
