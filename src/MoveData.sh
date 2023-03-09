cp *java ../gradle_directory/app/src/main/java/singlecellqc
cd ../gradle_directory
./gradlew build
cd -
cp ../gradle_directory/app/build/libs/singlecellqc.jar ../Jar/SingleCellQC.jar 
