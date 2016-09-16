#!/bin/tcsh

set WDir = `pwd`
set normal = ( "benzene" "bodipy" "T3" )
set prj = ( "benzene_prj" "bodipy_prj" "T3_prj" )
set logfile = "$WDir/logfile"
set failed = "0"
set passed = "0"

if ( -e $logfile ) then
    rm $logfile
endif

foreach test ( $normal )

  cd $test

  echo "Running test:  $test" >> $logfile

  ../../stable/JACoPO.py --chg1 mon1.chg --chg2 mon2.chg --cub1 mon1.cub --cub2 mon2.cub -o test.out
  sed -i '/Calculation time/d' test.out >& /dev/null
  diff test.out correct.out >& /dev/null

  if ( $status != "0" ) then
    echo "Test $test FAILED!" >> $logfile
    echo >> $logfile
    @ failed ++ 
  else
    echo "Test $test PASSED!" >> $logfile
    echo >> $logfile
    @ passed ++ 
  endif

  cd ..

end


foreach test ( $prj )

  cd $test

  echo "Running test:  $test" >> $logfile

  if ( $test == "bodipy_prj" ) then

    ../../stable/JACoPO.py --cub1 temp.cub --geo1 1.inc --cub2 temp.cub --selcub2 cub2sel.txt --geo2 2.inc --selgeo2 geo2sel.txt -o test.out
    sed -i '/Calculation time/d' test.out >& /dev/null
    diff test.out correct.out >& /dev/null

    if ( $status != "0" ) then
      echo "Test $test FAILED!" >> $logfile
      echo >> $logfile
      @ failed ++ 
    else
      echo "Test $test PASSED!" >> $logfile
      echo >> $logfile
      @ passed ++ 
    endif

    cd ..

  else
  
    ../../stable/JACoPO.py --cub1 temp.cub --geo1 1.inc --cub2 temp.cub --geo2 2.inc -o test.out
    sed -i '/Calculation time/d' test.out >& /dev/null
    diff test.out correct.out >& /dev/null

    if ( $status != "0" ) then
      echo "Test $test FAILED!" >> $logfile
      echo >> $logfile
      @ failed ++ 
    else
      echo "Test $test PASSED!" >> $logfile
      echo >> $logfile
      @ passed ++ 
    endif

    cd ..

  endif

end

echo "SUMMARY" >> $logfile
echo "TOTAL Tests:   6" >> $logfile
echo "PASSED Tests:  $passed" >> $logfile
echo "FAILED Tests:  $failed" >> $logfile
