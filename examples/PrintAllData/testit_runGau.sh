#! /bin/bash -x

cd ../PrintAllData
if [ -f outfile ]; then
    rm outfile
fi

# set up the Gaussian environment
export GauBINARY="${GAU_BINARY}"

mkdir test2
cd test2
echo "------------------------" >> ../outfile 2>&1
echo "Pass through MatFile test">> ../outfile
echo "------------------------" >> ../outfile 2>&1
cp ../../data/MatrixFile/rhf_h2-sto3g.mat . >> ../outfile
../DataPrint rhf_h2-sto3g.mat ${GauBINARY} >> ../outfile
echo "-------------------------------" >> ../outfile 2>&1
echo "test.com Create MatrixFile test">> ../outfile
echo "-------------------------------" >> ../outfile 2>&1
cp ../../data/Gaussian_input/test.com . >> ../outfile
../DataPrint test.com ${GauBINARY} >> ../outfile
cd ..
rm -r test2

mkdir test1071
cd test1071
echo "-------------------------------" >> ../outfile 2>&1
echo "test1071 Create MatrixFile test" >> ../outfile 2>&1
echo "-------------------------------" >> ../outfile 2>&1
#
# If mat file from DataSummary run is around, use it
#
if [ -d ../../DataSummary/test1071 ]; then
    cp ../../DataSummary/test1071/test1071.mat .
    ../DataPrint test1071.mat >> ../outfile 2>&1
else
    cp ../../data/Gaussian_input/test1071.com . >> ../outfile
    ../DataPrint test1071.com ${GauBINARY} >> ../outfile 2>&1
    mkdir ../../DataSummary/test1071
    cp test1071.mat ../../DataSummary/test1071/test1071.mat
    rm -f test1071.com
fi
cd ..
rm -r test1071

mkdir test1129
cd test1129
echo "-------------------------------" >> ../outfile 2>&1
echo "test1129 Create MatrixFile test" >> ../outfile 2>&1
echo "-------------------------------" >> ../outfile 2>&1
#
# If mat file from DataSummary run is around, use it
#
if [ -d ../../DataSummary/test1129 ]; then
    cp ../../DataSummary/test1129/test1129.dat .
    ../DataPrint test1129.dat >> ../outfile 2>&1
else
    cp ../../data/Gaussian_input/test1129.com . >> ../outfile
    ../DataPrint test1129.com ${GauBINARY} >> ../outfile 2>&1
    mkdir ../../DataSummary/test1129
    cp test1129.dat ../../DataSummary/test1129/test1129.dat
    rm -f test1129.com
fi
cd ..
rm -r test1129

mkdir test1132
cd test1132
echo "--------------------------------" >> ../outfile 2>&1
echo "test1132 Create MatrixFiles test" >> ../outfile 2>&1
echo "--------------------------------" >> ../outfile 2>&1
#
# If mat file from DataSummary run is around, use it
#
if [ -d ../../DataSummary/test1132 ]; then
    cp ../../DataSummary/test1132/test1132.mat .
    cp ../../DataSummary/test1132/test1132r.dat .
    cp ../../DataSummary/test1132/test1132gr.dat .
    cp ../../DataSummary/test1132/test1132u.dat .
    ../DataPrint test1132.mat >> ../outfile 2>&1
else
    cp ../../data/Gaussian_input/test1132.com . >> ../outfile
    ../DataPrint test1132.com ${GauBINARY} >> ../outfile 2>&1
    mkdir ../../DataSummary/test1132
    cp test1132.mat ../../DataSummary/test1132/test1132.mat
    cp test1132r.dat ../../DataSummary/test1132/test1132r.dat
    cp test1132gr.dat ../../DataSummary/test1132/test1132gr.dat
    cp test1132u.dat ../../DataSummary/test1132/test1132u.dat
    rm -f test1132.com
fi
echo "-----------------------------------------" >> ../outfile 2>&1
echo "test1132 Pass one MatrixFile through test" >> ../outfile 2>&1
echo "-----------------------------------------" >> ../outfile 2>&1
../DataPrint test1132r.dat >> ../outfile 2>&1
echo "------------------------------------------------------------------" >> ../outfile 2>&1
echo "test1132 Create One MartixFile and Pass through 3 MatrixFiles test" >> ../outfile 2>&1
echo "------------------------------------------------------------------" >> ../outfile 2>&1
../DataPrint test1132.mat test1132gr.dat test1132r.dat test1132u.dat ${GauBINARY} >> ../outfile 2>&1
cd ..
rm -r test1132

grep -v "echo argv" outfile > tmpfile 2>&1
grep -v "cp test1132" tmpfile > outfile 2>&1
grep -v "Gaussian Version:" outfile > tmpfile 2>&1
sed -e 'sZ-0.000000Z 0.000000Zg' < tmpfile > outfile
rm tmpfile

cd OUTPUT
rm -f out_runGau
var1=$(stat -c%s out_runGau_v2)
var2=$(stat -c%s ../outfile)
if [[ "$var1" -eq "$var2" ]]; then
    ln -s out_runGau_v2 out_runGau
else
    ln -s out_runGau_v1 out_runGau
fi
cd ..

 diff -b -B outfile OUTPUT/out_runGau

exit
# tests 1130, 1135 and 1136 require R.  So a pain to run.
cd test1130
echo "test1130-Should be no problems" >> ../outfile 2>&1
../FullWavefunctionRW test1130.com ${GauBINARY} >> ../outfile 2>&1
cd ..

cd test1135
echo "test1135" >> ../outfile
../FullWavefunctionRW test1135.com g16 >> ../outfile
cd ..

cd test1136
echo "test1136" >> ../outfile
../FullWavefunctionRW test1136.com g16 >> ../outfile
cd ..


