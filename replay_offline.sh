#called from within build dir

cmake .. -DCMAKE_BUILD_TYPE=Release \
  && make -j arcsim_0.2.1 \
  && pushd .. \
  && ./build-Release/bin/arcsim_0.2.1 replay ${1} \
  && popd
