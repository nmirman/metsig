#!/bin/bash

fileout="file_catalog.txt"
basedir="/mnt/xrootd/user/nmirman/Ntuples/METsig"

rm ${fileout}

for channel in `ls ${basedir}`
do
   for date in `ls ${basedir}/${channel}`
   do

      for process in `ls ${basedir}/${channel}/${date}`
      do
         filelist="${channel} ${date} ${process} "
         for file in `ls ${basedir}/${channel}/${date}/${process}`
         do
            filelist=${filelist}"${file} "
         done
         echo ${filelist} >> ${fileout}

      done

      echo >> ${fileout}
   done
done

