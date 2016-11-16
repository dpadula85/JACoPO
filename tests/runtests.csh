#!/bin/tcsh

set WDir = `pwd`
set normal = ( "benzene" "bodipy" "T3" )
set prj = ( "benzene_prj" "bodipy_prj" "T3_prj" )
set logfile = "$WDir/logfile"
set failed = "0"
set passed = "0"
set jac = "../../dev/JACoPO.py"
set total = "0"

if ( -e $logfile ) then
    rm $logfile
endif

foreach test ( $normal )

  cd $test

  echo "Running test:  $test" >> $logfile

  $jac --chg1 mon.chg --chg2 mon.chg --geo1 mon1.xyz --geo2 mon2.xyz --coup chgs -vvv > chgs.out
  $jac --cub1 mon1.cub --cub2 mon2.cub --coup tdc -vvv > tdc.out
  sed -i '/Calculation Time/d' chgs.out >& /dev/null
  sed -i '/Calculation Time/d' tdc.out >& /dev/null

  foreach out ( chgs.out tdc.out )

    diff $out ref.$out >& /dev/null

    if ( $status != "0" ) then
      echo "Test $test ${out:r} FAILED!" >> $logfile
      echo >> $logfile
      @ failed ++ 
    else
      echo "Test $test ${out:r} PASSED!" >> $logfile
      echo >> $logfile
      @ passed ++ 
    endif

    @ total ++

  end

  cd ..

end

foreach test ( $prj )

  cd $test

  echo "Running test:  $test" >> $logfile

  if ( $test == "bodipy_prj" ) then

    $jac --cub1 temp.cub --geo1 1.inc --cub2 temp.cub --selcub2 cub2sel.txt --geo2 2.inc --selgeo2 geo2sel.txt --coup tdc -vvv > tdc.out
    sed -i '/Calculation Time/d' tdc.out >& /dev/null
    diff tdc.out ref.tdc.out >& /dev/null

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
  
    $jac --cub1 temp.cub --geo1 1.inc --cub2 temp.cub --geo2 2.inc --coup tdc -vvv > tdc.out
    sed -i '/Calculation Time/d' tdc.out >& /dev/null
    diff tdc.out ref.tdc.out >& /dev/null

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

  @ total ++

end

echo "SUMMARY" >> $logfile
echo "TOTAL Tests: $total" >> $logfile
echo "PASSED Tests:  $passed" >> $logfile
echo "FAILED Tests:  $failed" >> $logfile
