Step 1:compile
```
mkdir build

cd build

g++ -I../include -O3 $(find ../src -name "*.cc") -lgflags -o my_program
```
Step 2:running BiT-DT

```./my_program --filepath="../dataset/yourdataset" --mode="standard"```

Step 3:running Batch-DT

```./my_program --filepath="../dataset/yourdataset" --mode="batch"```

All dataset files need to be compressed into a .bin file named graph-sort.bin and placed in the /dataset/yourdataset directory, where yourdataset is the name of the dataset. All datasets can be downloaded from http://konect.cc/networks/.
