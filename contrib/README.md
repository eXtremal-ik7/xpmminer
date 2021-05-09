How to deploy xpmclient:

1. Install docker (on Windows 10 also needed msys2 or different bash terminal)
2. Download latest CUDA 10 for Windows: https://developer.nvidia.com/cuda-downloads and put it to contrib directory
3. Setup CUDA10_INSTALLER variable in deploy.sh file:

  CUDA10_INSTALLER="cuda_10.0.130_411.31_win10.exe"
  
  In contrib 
  chmod +111 build.sh

4. run deploy.sh
