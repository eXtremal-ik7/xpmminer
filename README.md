xpmminer
========

Open-source primecoin(XPM) GPU & CPU miner (http://primecoin.io/). Only 

1. Requirements

- OS supported: Linux (Windows support coming soon), but you can try build on other Unix-like systems
- GPU manifacturer supported: AMD (NVidia is untested now)

2. Performance

Now miner can be slow than concurrents.. With this miner Radeon R9 290X @ 1120/1500 16.2 times faster than Core i7 920 @ 4.2GHz. It means, you can reach 1.05 chains/day on overclocked Radeon 290X.

3. Building

Miner uses CMake build system:

cd <project root directory>
mkdir build
cd build
cmake ../src
make
make install

But before this, you must fetch all dependencies:
 - gmp
 - OpenSSL
 - curl
 - jansson
 - OpenCL implementation

4. Usage

- Run official primecoin client with RPC support. Examples:<BR>
  primecoind -rpcuser=userName -rpcpassword=password<BR>
  primecoin-qt -rpcuser=userName -rpcpassword=password -server<BR>
  
- Run CPU or GPU miner:

./xpmclminer -u <userName> -p <password> -w Ac9ycgpEL4vzXRndS93Q7A7VGBHof1Jqzy (paste you own wallet)

for other options, see --help
