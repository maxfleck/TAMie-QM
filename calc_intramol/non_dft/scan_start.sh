#!/bin/bash
dfne() { 
     if [ ! -f thisdef2 ]; then 
       sdg intdef > thisdef2 
     fi 
     cp coord coord.back 
     kdg cartesianstep 
     kdg intdef 
     kdg end 
     cat control > thiscontrol 
     cat thisdef2 >> thiscontrol 
     echo '$end' >> thiscontrol 
     mv thiscontrol control 
     define < define.in > define2.out 2>&1 
     grep 'CAUTION, within' define2.out > /dev/null 2>&1 
     if [ $? -eq 0 ]; then
      echo 'WARNING: define went wrong! check '$i'/define2.out ' 
      echo 'trying again with thisdef ' 
      cp coord.back coord 
      kdg intdef 
      kdg end 
      cat control > thiscontrol 
      cat thisdef >> thiscontrol 
      echo '$end' >> thiscontrol 
      mv thiscontrol control 
      define < define.in > define3.out 2>&1 
      grep 'CAUTION, within' define3.out > /dev/null 2>&1 
      if [ $? -eq 0 ]; then
          echo 'ERROR: define went wrong! check '$i'/define3.out ' 
          cp coord.back coord 
          cd ..  
          continue 
      fi 
     fi 
} 
for i in `cat scan_list`
do
 cd $i
 if [ -f GEO_OPT_CONVERGED ]; then 
   echo 'found GEO_OPT_CONVERGED in '$i'. skipping...' 
   cd .. 
   continue 
 fi 
 if [ -f GEO_OPT_RUNNING ]; then 
   echo 'found GEO_OPT_RUNNING in '$i'. something is strange. skipping...' 
   cd .. 
   continue 
 fi 
 if [ -f stop -o -f STOP ]; then 
   echo 'ERROR: found stop file in '$i'. something is strange. skipping...' 
   echo 'try the concatenated scan on this structure. '
   cd .. 
   continue 
 fi 
 echo 'running in '$i
 if [ ! -f GEO_OPT_FAILED ]; then 
     if [ ! -f thisdef2 ]; then 
       sdg intdef > thisdef2 
     fi 
     dfne 
 fi 
 if [ ! -f GEO_OPT_FAILED ]; then 
   actual -r
   finit
   echo ' initial run 1' >> jobex.out 
   rm -f hessapprox nextstep optinfo tmp.input 
   jobex -level cc2 -c 10 >> jobex.out 
   dfne 
   actual -r
   finit
   echo ' initial run 2' >> jobex.out 
   rm -f hessapprox nextstep optinfo tmp.input 
   jobex -level cc2 -c 10 >> jobex.out 
   dfne 
 fi 
 finit
 rm -f hessapprox nextstep optinfo tmp.input 
     echo ' normal run ' >> jobex.out 
 jobex -level cc2 >> jobex.out 
  if [ -f GEO_OPT_FAILED ]; then 
   echo '      did not yet converge!' 
   echo '      check for the reason in '$i 
   echo '      and restart the script to continue.' 
  fi
 t2x -c > $i.xyz
 cd ..
done
rm -rf scan_movie.xyz 
rm -rf scan_energy.dat 
for i in `cat scan_list`
do
 cat $i/$i.xyz >> scan_movie.xyz 
 echo -n $i' : ' >> scan_energy.dat 
 grep cyc $i/gradient |awk '{print $7}'|tail -n 1 >> scan_energy.dat 
done
sed 's/_/ /g' scan_energy.dat > scan_tmp
mv scan_tmp scan_energy.dat 
