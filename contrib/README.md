How to deploy xpmclient:

1. Install docker (on Windows 10 also needed msys2 or different bash terminal)
2. Download latest CUDA 10 for Windows: https://developer.download.nvidia.com/compute/cuda/11.2.1/local_installers/cuda_11.2.1_461.09_win10.exe and put it to contrib directory
3. Setup CUDA10_INSTALLER variable in deploy.sh file:

  CUDA10_INSTALLER="cuda_11.2.1_461.09_win10.exe"
  
  In contrib 
  chmod +111 build.sh

4. run deploy.sh
