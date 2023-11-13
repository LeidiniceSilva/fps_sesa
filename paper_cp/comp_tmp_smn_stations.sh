#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = 'Nov 13, 2023'
#__description__ = 'Compute mean temperature'

echo
echo "--------------- THE INIT ----------------"

# Code list
code_list=('87007' '87078' '87047' '87121' '10038' '87222' '87129' '87211' '10067' '87217' '87257' '87305'	'87393' '87344' 
'87345' '87349' '87374' '10114' '87328' '87418' '87420' '87480' '87497' '87453' '87436' '10145' '87534' '87569' '87585' '87571' 
'87576' '87593' '87270' '87509' '87022' '87582' '87448' '87016' '87289' '87097' '10358' '87178' '87467' '87416' '10423' '10444' 
'10445' '87371' '87548' '10460' '10464' '87166' '87395' '87162' '87155' '10492' '87311' '87046' '10502' '10550' '10901' '10902' 
'10903' '10904' '10905' '10906' '10907' '10908' '10909' '10910' '10911' '10912' '86330' '86580' '86560' '86530' '11004' '86440' 
'86490' '86460' '86430' '86585' '86350' '86565' '86360' '86500' '11014' '11015' '11016' '11017' '11018' '11019' '11020' '86068' 
'86086' '86097' '86128' '86134' '86185' '86192' '86210' '86218' '86233' '86234' '86255' '86260' '86268' '86285' '86297')

dt="D_1979-01-01_2021-12-31"
path_tmax="/afs/ictp.it/home/m/mda_silv/Documents/FPS_SESA/database/obs/smn_ii/smn_nc/tmax/"
path_tmin="/afs/ictp.it/home/m/mda_silv/Documents/FPS_SESA/database/obs/smn_ii/smn_nc/tmin/"
path_tmp="/afs/ictp.it/home/m/mda_silv/Documents/FPS_SESA/database/obs/smn_ii/smn_nc/tmp/"

for code in ${code_list[@]}; do

    echo "Compute the mean temperature:" ${code}
    cdo add ${path_tmax}tmax_${code}_${dt}.nc ${path_tmin}tmin_${code}_${dt}.nc ${path_tmp}tmax_tmin_${code}_${dt}.nc 
    cdo divc,2 ${path_tmp}tmax_tmin_${code}_${dt}.nc ${path_tmp}tmp_${code}_${dt}.nc 
    
    echo "Delete files"
    rm ${path_tmp}tmax_tmin_${code}_${dt}.nc
    
done

echo
echo "--------------- THE END ----------------"


