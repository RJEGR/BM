# Procesando datos CTD

0) 

```bash
head -n31 lan001b.cnv
# Example
# name 0 = scan: Scan Count
# name 1 = prDM: Pressure, Digiquartz [db]
# name 2 = t090C: Temperature [ITS-90, deg C]
# name 3 = c0mS/cm: Conductivity [mS/cm]
# name 4 = sal00: Salinity, Practical [PSU]
# name 5 = sbeox0ML/L: Oxygen, SBE 43 [ml/l]
# name 6 = flSP: Fluorescence, Seapoint
# name 7 = wetCDOM: Fluorescence, WET Labs CDOM [mg/m^3]
# name 8 = CStarTr0: Beam Transmission, WET Labs C-Star [%]
# name 9 = haardtY: Fluorescence, Dr. Haardt Yellow Sub
# name 10 = sigma-t00: Density [sigma-t, kg/m^3 ]
# name 11 = sbeox0Mm/Kg: Oxygen, SBE 43 [umol/kg], WS = 2
# name 12 = latitude: Latitude [deg]
# name 13 = longitude: Longitude [deg]
# name 14 = timeJ: Julian Days
# name 15 = altM: Altimeter [m]
# name 16 = nbin: number of scans per bin
# name 17 = user1: cast number
# interval = decibars: 1
# file_type = ascii
# bad_flag = NaN
# *END*
```



1) Relabel files by Cruice_Station grep

```bash
#!bin/bash

for i in *.cnv
do
   S=$(grep Station $i | sed 's/** Station://'g)
   C=$(grep Cruise $i | sed 's/** Cruise://g')
   echo 'labeling' $i 'by' $C 'with' $S
   cp $i ${i%.cnv}_${C}_${S}.tmp
done

```

2) Clean headers from tmp lanefiles and reformat to csv

```bash
for i in *.tmp
do
   awk 'NR>31' $i > ${i%.tmp}.csv
done && rm *.tmp

```

3) Lets sort and uniq stations in order to define repetition mesurements:

```bash
ls *.csv | cut -d'_' -f3 | sort |uniq -c | sort -n

# Ejemplo:

#   3 G44.csv
#   2 Y6.csv
#   3 Y7.csv
# And check unique samples
```

4) And remove repetition based on the less deep info

```bash
# 1)
# Example
for i in G44; do wc -l *${i}.csv ; done | sort -k2,2

#      51 lan022b_XIX06_G44.csv
#    2373 lan023b_XIX06_G44.csv <--- SELECT IT!
#    1010 lan024b_XIX06_G44.csv

# 2)
ls *.csv | cut -d'_' -f3 | sort |uniq -c | sort -n

   1 G44.csv
   1 Y6.csv
   1 Y7.csv
```

5) Subsetting data <= 200m and select four intentioned variables 

- name 12 = latitude: **Latitude** [deg] - ($13)
- name 13 = longitude: **Longitude** [deg] - ($14)
- name 1 = prDM: **Pressure**, Digiquartz [db] - ($2)
- name 2 = t090C: **Temperature** [ITS-90, deg C] - ($3)
- name 4 = sal00: **Salinity**, Practical [PSU] - ($5)
- name 5 = sbeox0ML/L: **Oxygen**, SBE 43 [ml/l] - ($6) 
- name 6 = flSP: **Fluorescence**, Seapoint -  ($7)
- name 10 = sigma-t00: **Density** [sigma-t, kg/m^3 ] - ($11)

```bash
for i in *.csv
do
   awk '$2 <= 200.0 {print $13, $14, $2, $3, $5, $6, $7, $11}' $i > ${i%.csv}_200m.csv
done
```



6) Calculate the average of data:

> saving the follow script as `mean.awk`

```bash
{
    for(i=1; i<=NF; i++) {
        a[i]+=$i
        if($i!="") 
            b[i]++}
    } 
END {
    for(i=1; i<=NF; i++) 
        printf "%s%s", a[i]/b[i], (i==NF?ORS:OFS)
}
```

And run

```bash
# awk -f mean.awk lan017b_XIX06_J49_200m.csv

# Loop

for i in *_200m.csv
do
	f=$(ls $i | cut -d"_" -f3)
	m=$(awk -f mean.awk $i)
	echo "$f $m"
done >> lances_medias.csv


```

Other approach to calculate mean:

```bash
awk '{x+=$4; next} END{print x/NR}' lan010b_XIX06_B11_200m.csv


for i in *_200m.csv
do
	f=$(ls $i | cut -d"_" -f3)
	m=$(awk '{x+=$2; next} END{print x/NR}' $i)
	echo "$f $m"
done
	
```



