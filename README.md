xpmminer
========

Open-source primecoin(XPM) GPU & CPU miner (http://primecoin.io/). Only solo mining now available.

1 Requirements

- OS supported: Linux, Windows
- GPU manifacturer supported: NVidia (AMD is untested now)

2 Performance

Now miner can be slow than concurrents.. With this miner Radeon R9 290X @ 1120/1500 16.2 times faster than Core i7 920 @ 4.2GHz. It means, you can reach 1.05 chains/day on overclocked Radeon 290X.

3 Building

in linux see contrib/README.md

4 Usage

- Run official primecoin client with RPC support. Examples:<BR>
  primecoind -rpcuser=userName -rpcpassword=password<BR>
  primecoin-qt -rpcuser=userName -rpcpassword=password -server<BR>
  
- Run CPU or GPU miner,
method 1
./miner --url RPCaddress --user primecoinrpc --pass PASSWORD   --wallet youaddress
for other options, see --help

method 2
modify xpmminercpu (in linux cpu) / xpmminernv (in linux NVidia) / xpmminernv.bat (in windows NVidia)
use you RPC address, rpcuser, rpcpassword and you wallet address to replace the corresponding contents. Then run these files. For example, in windows, run `.\xpmminernv.bat` in power shell.


5 Donations

BTC 17TQurzvatmsqZzGEe8jXnMDoiC4ACZjH7<BR>
LTC LdZcF4WejhC46DHfdyjvXomZzqQRy5xCj2<BR>
XPM Ac9ycgpEL4vzXRndS93Q7A7VGBHof1Jqzy
