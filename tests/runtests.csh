#!/bin/tcsh

set normal = ( "benzene" "bodipy" "T3" )
set prj = ( "benzene_prj" "bodipy_prj" "T3_prj" )


foreach test ( $normal )

  cd $test
  ../../stable/JACoPO.py --chg1 mon1.chg --chg2 mon2.chg --cub1 mon1.cub --cub2 mon2.cub -o test.out
  sed -i '/Calculation time/d' test.out
  diff test.out correct.out >& /dev/null

  if ( $status != "0" ) then
      echo
      echo "Test $test FAILED!"
      echo
  endif

  cd ..

end


foreach test ( $prj )

  cd $test

  if ( $test == "bodipy_prj" ) then

    ../../stable/JACoPO.py --cub1 temp.cub --geo1 1.inc --cub2 temp.cub --selcub2 cub2sel.txt --geo2 2.inc --selgeo2 geo2sel.txt -o test.out
    sed -i '/Calculation time/d' test.out
    diff test.out correct.out >& /dev/null

    if ( $status != "0" ) then
        echo
        echo "Test $test FAILED!"
        echo
    endif

    cd ..

  else
  
    ../../stable/JACoPO.py --cub1 temp.cub --geo1 1.inc --cub2 temp.cub --geo2 2.inc -o test.out
    sed -i '/Calculation time/d' test.out
    diff test.out correct.out >& /dev/null

    if ( $status != "0" ) then
        echo
        echo "Test $test FAILED!"
        echo
    endif

    cd ..

  endif

end
