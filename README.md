xpmminer
========

Open-source primecoin(XPM) GPU & CPU miner (http://primecoin.io/). Only solo mining now available.

1 Requirements

- OS supported: Linux (Windows support coming soon), but you can try build on other Unix-like systems
- GPU manifacturer supported: AMD (NVidia is untested now)

2 Performance

Now miner can be slow than concurrents.. With this miner Radeon R9 290X @ 1120/1500 16.2 times faster than Core i7 920 @ 4.2GHz. It means, you can reach 1.05 chains/day on overclocked Radeon 290X.

3 Building

3.1 Prepare to build

But before this, you must fetch all dependencies:
 - gmp
 - OpenSSL
 - curl
 - jansson
 - ncurses
 - OpenCL SDK
 - GPU drivers (Catalyst 14.1 or later)

AMD APP SDK contains broken OpenCL compiler, it shows segmentation fault when kernel.cl compiling. Workaround is using compiler (libamdocl64.so) from Catalyst package. On Ubuntu you can run this shell command:

sudp cp /usr/lib/fglrx/libamdocl64.so /opt/AMDAPP/lib/x86_64
 - /usr/lib/fglrx - driver installation path
 - /opt/AMDAPP - AMD APP SDK installation path


3.2 Build

Miner uses CMake build system.

cd 'project root directory'<BR>
mkdir build<BR>
cd build<BR>
cmake -DCMAKE_BUILD_TYPE=Release ../src<BR>
make<BR>
make install

CMake may be failed found OpenCL SDK. In this case, setup path to SDK manually:

cmake ../src -DOPENCL_INCLUDE_DIRS=/opt/AMDAPP/include -DOPENCL_LIBRARIES=/opt/AMDAPP/lib/x86_64/libOpenCL.so -DCMAKE_BUILD_TYPE=Release

4 Usage

- Run official primecoin client with RPC support. Examples:<BR>
  primecoind -rpcuser=userName -rpcpassword=password<BR>
  primecoin-qt -rpcuser=userName -rpcpassword=password -server<BR>
  
- Run CPU or GPU miner:

./xpmclminer -u <userName> -p <password> -w Ac9ycgpEL4vzXRndS93Q7A7VGBHof1Jqzy (paste you own wallet)

for other options, see --help

5 Donations

BTC 17TQurzvatmsqZzGEe8jXnMDoiC4ACZjH7<BR>
LTC LdZcF4WejhC46DHfdyjvXomZzqQRy5xCj2<BR>
XPM Ac9ycgpEL4vzXRndS93Q7A7VGBHof1Jqzy
