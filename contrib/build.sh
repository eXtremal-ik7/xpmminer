#!/bin/bash
set -e
VERSION="10.5-beta2"

# Linux static build

# gmp
cd /home/user/build/deps-linux
tar --lzip -xvf ../gmp-6.1.2.tar.lz
cd gmp-6.1.2
./configure --build corei7 --prefix=/home/user/install/x86_64-Linux --enable-cxx --enable-static --disable-shared 
make -j`nproc`
make install

# openssl
cd /home/user/build/deps-linux
tar -xzf ../openssl-1.0.2.tar.gz
cd openssl-1.0.2
./config --prefix=/home/user/install/x86_64-Linux no-shared
make
make install

# curl
cd /home/user/build/deps-linux
tar -xzf ../curl-7.68.0.tar.gz
cd curl-7.68.0
./configure --prefix=/home/user/install/x86_64-Linux --enable-static --disable-shared 
make -j`nproc`
make install

# jansson
cd /home/user/build/deps-linux
tar -xzf ../jansson-2.11.tar.gz
cd jansson-2.11
./configure --prefix=/home/user/install/x86_64-Linux --enable-static --disable-shared
make -j`nproc`
make install

# CLRX
mkdir $HOME/build/deps-linux/CLRX
cd $HOME/build/deps-linux/CLRX
cmake $HOME/build/CLRX-mirror -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$HOME/install/x86_64-Linux
make -j`nproc`
make install
rm $HOME/install/x86_64-Linux/lib64/libCLRX*.so*

# xpmminer
mkdir /home/user/build/xpmminer/x86_64-Linux
cd /home/user/build/xpmminer/x86_64-Linux
cmake ../src -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/home/user/install/x86_64-Linux \
  -DSTATIC_BUILD=ON \
  -DCUDA_DRIVER_LIBRARY=/usr/local/cuda-11.2/compat/libcuda.so \
  -DBUILDOPENCLMINER=OFF
make -j`nproc`

# make CPU distr
mkdir xpmminer-cpu-$VERSION-linux
cd xpmminer-cpu-$VERSION-linux
cp ../CPU/xpmcpuminer ./miner
echo "./miner --url primea.primecoin.org:9912 --user primecoinrpc --pass PASSWORD   --wallet AbMXJscCQUP28rySq8NYorqHouEo5zF8eH" >> xpmminercpu
chmod +x xpmminercpu
cd ..
tar -czf xpmminer-cpu-$VERSION-linux.tar.gz xpmminer-cpu-$VERSION-linux

# make NVidia distr
mkdir xpmminer-cuda-$VERSION-linux
cd xpmminer-cuda-$VERSION-linux
cp ../Cuda/xpmcuda ./miner
echo "#/bin/bash" > xpmminernv
echo "DIR=\$(dirname \"\$0\")" >> xpmminernv
echo "LD_LIBRARY_PATH=\$DIR/. ./miner --url primea.primecoin.org:9912 --user primecoinrpc --pass PASSWORD   --wallet AbMXJscCQUP28rySq8NYorqHouEo5zF8eH" >> xpmminernv
chmod +x xpmminernv
mkdir -p xpm/cuda
cp ../../src/Cuda/*.cu xpm/cuda
cp /usr/local/cuda-11.2/lib64/libnvrtc.so.11.2 .
cp /usr/local/cuda-11.2/lib64/libnvrtc-builtins.so.11.2 .
cd ..
tar -czf xpmminer-cuda-$VERSION-linux.tar.gz xpmminer-cuda-$VERSION-linux

# win64 static build

# gmp
cd /home/user/build/deps-win32
tar --lzip -xvf ../gmp-6.1.2.tar.lz
cd gmp-6.1.2
./configure --host=x86_64-w64-mingw32 --build corei7 --prefix=/home/user/install/x86_64-w64-mingw32 --enable-cxx --enable-static --disable-shared
make -j`nproc`
make install

# openssl
cd /home/user/build/deps-win32
tar -xzf ../openssl-1.0.2.tar.gz
cd openssl-1.0.2
sed -i 's/:.dll.a/ -Wl,--export-all -shared:.dll.a/g' Configure
sed -i 's,.*target already defined.*,$target=$_;,g' Configure
./Configure --cross-compile-prefix="x86_64-w64-mingw32-" mingw64  --prefix=/home/user/install/x86_64-w64-mingw32 shared
make
make install

# curl
cd /home/user/build/deps-win32
tar -xzf ../curl-7.68.0.tar.gz
cd curl-7.68.0
./configure  --host=x86_64-w64-mingw32 --prefix=/home/user/install/x86_64-w64-mingw32
make -j`nproc`
make install

# jansson
cd /home/user/build/deps-win32
tar -xzf ../jansson-2.11.tar.gz
cd jansson-2.11
./configure --host=x86_64-w64-mingw32 --prefix=/home/user/install/x86_64-w64-mingw32 --enable-static --disable-shared
make -j`nproc`
make install

# CLRX
mkdir $HOME/build/deps-win32/CLRX
cd $HOME/build/deps-win32/CLRX
cmake $HOME/build/CLRX-mirror -DCMAKE_BUILD_TYPE=Release -DCMAKE_TOOLCHAIN_FILE=$HOME/build/xpmminer/src/cmake/Toolchain-cross-mingw32-linux.cmake -DCMAKE_INSTALL_PREFIX=$HOME/install/x86_64-w64-mingw32
make -j`nproc`
make install
rm -f $HOME/install/x86_64-w64-mingw32/lib/libCLRX*.dll*

# xpmminer
mkdir /home/user/build/xpmminer/x86_64-w64-mingw32
cd /home/user/build/xpmminer/x86_64-w64-mingw32
cmake ../src -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_TOOLCHAIN_FILE=../src/cmake/Toolchain-cross-mingw32-linux.cmake \
  -DCMAKE_INSTALL_PREFIX=/home/user/install/x86_64-w64-mingw32 \
  -DSTATIC_BUILD=ON \
  -DCUDA_TOOLKIT_ROOT_DIR=/usr/local/cuda-win32 \
  -DOpenCL_INCLUDE_DIR=/usr/local/cuda-win32/include \
  -DOpenCL_LIBRARY=/usr/local/cuda-win32/lib/x64/OpenCL.lib \
  -DBUILDOPENCLMINER=OFF -DBUILDCPUMINER=OFF \
  -DCUDA_INCLUDE_DIRS=/usr/local/cuda-win32/include  \
  -DCUDA_DRIVER_LIBRARY=/usr/local/cuda-win32/lib/x64/cuda.lib \
  -DCUDA_nvrtc_LIBRARY=/usr/local/cuda-win32/lib/x64/nvrtc.lib \
  -DOPENSSL_CRYPTO_LIBRARY=/home/user/install/x86_64-w64-mingw32/lib/libcrypto.dll.a \
  -DOPENSSL_SSL_LIBRARY=/home/user/install/x86_64-w64-mingw32/lib/libssl.dll.a
make -j`nproc`
cd Cuda
x86_64-w64-mingw32-strip xpmcuda.exe

# make Nvidia distr
mkdir xpmminer-cuda-$VERSION-win64
cd xpmminer-cuda-$VERSION-win64
cp ../xpmcuda.exe miner.exe
mkdir -p xpm/cuda
cp ../../../src/Cuda/*.cu xpm/cuda/
cp /usr/lib/gcc/x86_64-w64-mingw32/7.3-posix/libgcc_s_seh-1.dll .
cp /usr/lib/gcc/x86_64-w64-mingw32/7.3-posix/libstdc++-6.dll .
cp /usr/x86_64-w64-mingw32/lib/libwinpthread-1.dll .
cp /usr/local/cuda-win32/bin/nvrtc64_112_0.dll ./
cp /usr/local/cuda-win32/bin/nvrtc-builtins64_112.dll .
cp /home/user/install/x86_64-w64-mingw32/bin/libcurl-4.dll ./
cp /home/user/install/x86_64-w64-mingw32/bin/libeay32.dll ./
echo ".\miner.exe --url primea.primecoin.org:9912 --user primecoinrpc --pass PASSWORD   --wallet AbMXJscCQUP28rySq8NYorqHouEo5zF8eH" > xpmminernv.bat

cd ..

zip -9 -r xpmminer-cuda-$VERSION-win64.zip xpmminer-cuda-$VERSION-win64
# Calculate SHA256 checksum
cd /home/user/build/xpmminer
if [ -f /home/user/build/xpmminer/x86_64-Linux/xpmminer-cpu-$VERSION-linux.tar.gz ]; then
sha256sum /home/user/build/xpmminer/x86_64-Linux/xpmminer-cpu-$VERSION-linux.tar.gz >> xpmminer-$VERSION-sha256.txt
fi
if [ -f /home/user/build/xpmminer/x86_64-Linux/xpmminer-cuda-$VERSION-linux.tar.gz ]; then
sha256sum /home/user/build/xpmminer/x86_64-Linux/xpmminer-cuda-$VERSION-linux.tar.gz >> xpmminer-$VERSION-sha256.txt
fi
if [ -f /home/user/build/xpmminer/x86_64-w64-mingw32/xpmminer-cuda-$VERSION-win64.zip ]; then
sha256sum /home/user/build/xpmminer/x86_64-w64-mingw32/xpmminer-cuda-$VERSION-win64.zip >> xpmminer-$VERSION-sha256.txt
fi