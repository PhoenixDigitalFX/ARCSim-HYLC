#called from within build dir
testcase=${1:-0} # get argument, defaulting to "0"

if [[ $testcase = *.json ]]
then # select by filename ending with json
  configfile=$testcase
else # select by id
  if [ $testcase == 1 ] ;
  then
    configfile=conf/hylc_sphere.json
  else
    configfile=conf/hylc_sphere_noremesh.json
  fi
fi

cmake .. -DCMAKE_BUILD_TYPE=Release \
  && make -j arcsim_0.2.1 \
  && pushd .. \
  && ./build-Release/bin/arcsim_0.2.1 simulate $configfile \
  && popd
