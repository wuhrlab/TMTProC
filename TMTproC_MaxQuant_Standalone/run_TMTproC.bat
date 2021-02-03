ECHO OFF
monocle_converter\Monocle.CLI.exe -f data\%1.raw -t exmzxml -x
python3 reporter_quant_python\main.py --peptides data\%2.csv --scans data\%1.mzXML --quant parameters\quant.ini
matlab_functions\TMTproC_maxquant_compatible.exe data\%2_intermediate.csv