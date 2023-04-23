#!/bin/bash

#__author__      = 'Leidinice Silva'
#__email__       = 'leidinicesilva@gmail.com'
#__date__        = '04/17/2023'
#__description__ = 'Select stations from the code list'

echo
echo "--------------- INIT ----------------"

# Variable
var='pre'

# Code list
code_list=( 'A899' 'A836' 'A802' 'A811' 'A827' 'A881' 'A804' 'A838' 'A812' 'A831' 'A832' 'A801' 'A886' 'A813' 'A809' 'A826' 'A803' 'A889' 'A884' 'A879' 'A808' 'A833' 'A840' 'A897' 'A867' 'A837' 'A829' 'A894' 'A883' 'A830' 'A866' 'A853' 'A814' 'A880' 'A852' 'A815' 'A839' 'A844' 'A845' 'A856' 'A810' 'A805' 'A865' 'A870' 'A828' 'A806' 'A863' 'A854' 'A860' 'A841' 'A868' 'A858' 'A817' 'A859' 'A876' 'A875' 'A848' 'A862' 'A874' 'A823' 'A873' 'A807' 'B804' 'B806' 'A818' 'A820' 'A714' 'A715' 'A871' 'A821' 'A752' 'A835' 'A869' 'A850' 'A718' 'A705' 'A763' 'A707' 'A768' 'A743' 'A727' 'A762' 'A747' 'A734' 'A748' 'A702' 'A756' 'A729' 'A520' 'A519' 'A724' 'A507' 'A616' 'A035' 'A016' 'A003' 'A026' 'A909' 'A029' 'A033' )

for code in ${code_list[@]}; do

    echo "Starting to move:" ${code}
    mv ${var}_${code}_H_* /home/nice/Documentos/FPS_SESA/weather_stations/inmet_nc/${var}
    
done

echo
echo "------------- THE END --------------"

