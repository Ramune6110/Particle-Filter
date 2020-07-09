# Particle-Filter
Particle-Filter
AtushiさんのMATLABcodeを参考にC++にてPFを実装しました.  
https://myenigma.hatenablog.com/entry/20140628/1403956852
C++で演算した結果をtxtファイルに保存して, 描画はMATLABで行っています.
![PF](https://github.com/Ramune6110/Particle-Filter/blob/master/Particle_Filter.png)
## Environment
Ubuntu18.04
## Procedure
```bash
g++ main.cpp particle_filter.cpp -I /usr/include/eigen3
```
```bash
./a.out
```
